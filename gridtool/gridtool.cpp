#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "grid.hpp"
#include "gridutil.hpp"
#include "smoothgrid.hpp"
#include "bltmatrx.hpp"
#include "get_image_path.h"

using namespace std;

static void run_commands( list<string> &commands );

typedef list<string> commandlist;

void help()
{
    string helpfile(get_image_path());
    helpfile += ".help";
    ifstream hf(helpfile.c_str());
    if( hf)
    {
        string line;
        while(hf)
        {
            getline(hf,line);
            cout << line << endl;
        }
        hf.close();
        exit(0);
    }
}


static string &lowercase( string &s )
{
    transform(s.begin(),s.end(),s.begin(),::tolower);
    return s;
}

template<class T> bool parse_string( const string &str, T &value )
{
    istringstream s(str);
    char c;
    bool result = true;
    if( !(s >> value) || s.get(c)  ) result = false;
    // if( !(s >> value) ) result = false;
    return result;
}


void print_syntax()
{
    help();
    cout << "Arguments: input_grid_file [actions] output_grid_file" << endl;
    cout << "       or: command_file" << endl;
    exit(0);

}

int main( int argc, char *argv[] )
{
    commandlist commands;

    if( argc < 2 )
    {
        print_syntax();
    }

    for( int i = 1; i < argc-1; i++ ) commands.push_back(argv[i]);
    if( commands.size() == 0 ) commands.push_front("run");
    if( commands.front() != "run" && commands.front() != "read" && commands.front() != "create" )
    {
        commands.push_front("read");
        commands.push_back("write");
    }
    commands.push_back(argv[argc-1]);
    run_commands(commands);
}

vector<string> split(const string &s, char delim ) 
{
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) 
    {
         elems.push_back(item);
    }
    return elems;
}

static string next_command( commandlist &commands, const string reason="" )
{
    string result("");
    if( commands.size() )
    {
        result = commands.front();
        commands.pop_front();
    }
    else if( reason != "")
    {
        throw runtime_error(reason+" missing");
    }
    return result;
}

bool next_command_is( commandlist &commands, const string &value )
{
    bool result = false;
    if( commands.size() > 0 && commands.front() == value )
    {
        commands.pop_front();
        result = true;
    }
    return result;
}

template<class T> bool next_command_value( commandlist &commands, T &value, const string reason="" )
{
    bool result = false;
    if( commands.size() )
    {
        string svalue = commands.front();
        result = parse_string(svalue,value);
        if( result )
        {
            commands.pop_front();
        }
        else if( reason != "" )
        {
            throw runtime_error(reason+": Invalid value " + svalue);
        }
    }
    else if( reason != "")
    {
        throw runtime_error(reason+" missing");
    }
    return result;
}

static string run_read_grid( grid &g, commandlist &commands, const string &operation="read" )
{
    int maxcols = 99;
    char delimiter=' ';
    string filename;
    bool output=operation == "read";

    while( true )
    {
        filename = next_command(commands,string("Filename for ")+operation+" operation");
        if( filename == "maxcols" )
        {
           next_command_value(commands,maxcols,"Maximum number of columns for read");
        }
        else if( filename == "csv" )
        {
            delimiter=',';
        }
        else
        {
            break;
        }
    }
    if( output ) cout << "Reading file from " << filename << endl;
    try
    {
        g.readfile(filename.c_str(),delimiter,maxcols);
    }
    catch( runtime_error e )
    {
        string message = filename + " " + e.what();
        throw(runtime_error( message ));
    }
    if( output ) cout << "Grid has " << g.nrow() << " rows and " << g.ncol() << " columns" << endl;
    if( output ) cout << "Each point has " << g.nvalue() << " data values" << endl;
    return filename;
}

static void run_create_grid( grid &g, commandlist &commands )
{
    double minx, maxx, miny, maxy, incx, incy;
    vector<string> columns;
    if( next_command_is(commands,"extents"))
    {
        if( next_command_is(commands,"wkt"))
        {
            string filename = next_command( commands, "Extents wkt_file file for ");
            wkt_extents(filename.c_str(),minx,maxx,miny,maxy);
        }
        else
        {
            grid gext;
            run_read_grid( gext, commands, "create" );
            gext.extents(minx,maxx,miny,maxy);
        }
        next_command_value(commands,incx,"Grid longitude spacing");
        next_command_value(commands,incy,"Grid latitude spacing");
        minx = incx*floor(minx/incx);
        miny = incy*floor(miny/incy);
    }
    else
    {
        next_command_value(commands,minx,"Grid minimum longitude");
        next_command_value(commands,maxx,"Grid maximum longitude");
        next_command_value(commands,incx,"Grid longitude spacing");

        next_command_value(commands,miny,"Grid minimum latitude");
        next_command_value(commands,maxy,"Grid maximum latitude");
        next_command_value(commands,incy,"Grid latitude spacing");
    }

    if( next_command_is(commands,"columns") )
    {
        string columndef=next_command(commands,"Grid columns to create");
        columns=split(columndef,'+');
    }
    g.create( minx, maxx, incx, miny, maxy, incy, columns );
}

static void mark_grid( grid &g, commandlist &commands, string markcommand, grid *refgrid=0 )
{
    grid::markaction action = grid::on;
    g.clearMarked();


    while( true )
    {
        string command = next_command( commands, string("Selection action for ") + markcommand);
        if( command == "nearest_to" )
        {
            string filename = next_command( commands, string("Points file for ") + markcommand);
            mark_xy_file( g, filename.c_str(), action );
        }
        else if( command == "on_grid" )
        {
            grid *rgrid;
            grid ongrid;
            if( refgrid && next_command_is(commands,refgrid->filename()))
            {
                rgrid=refgrid;
            }
            else
            {
                string filename = run_read_grid( ongrid, commands, markcommand+" on_grid" );
                rgrid=&ongrid;
            }
            grid::point xy;
            for( int row = 0; row < rgrid->nrow(); row++ ) 
                for( int col = 0; col < rgrid->ncol(); col++ )
                {
                    rgrid->nodexy(row,col,xy.x,xy.y);
                    g.markNearest(xy,action,0.0001);
                }
        }
        else if( command == "inside" )
        {
            string filename = next_command( commands, string("Inside wkt_file file for ") + markcommand);
            mark_wkt_file( g, false, filename.c_str(), action );
        }
        else if( command == "outside" )
        {
            string filename = next_command( commands, string("Outside wkt_file file for ") + markcommand);
            mark_wkt_file( g, true, filename.c_str(), action );
        }
        else if( command == "edge" )
        {
            int nedge = 1;
            next_command_value(commands,nedge,string("Width of edge for ")+markcommand);
            g.markEdge(nedge,false,action);
        }
        else if( command == "inside_edge" )
        {
            int nedge = 1;
            next_command_value(commands,nedge,string("Width of inside_edge for ")+markcommand);
            g.markEdge(nedge,true,action);
        }
        else if( command == "expand" && action == grid::on )
        {
            double distance = 0.0;
            next_command_value(commands,distance,string("Distance for expand for ")+markcommand);
            expand_marked(g,distance,true);
        }
        else 
        {
           string field = command;
           if( field == "where" ) field = next_command( commands, string("Field name for ")+markcommand);
           string op = next_command( commands, string("Operation type for ")+markcommand);
           double value;
           next_command_value(commands,value,string("Field value for ")+markcommand);
           g.markWhere(field,op,value,action);
        }
        if( next_command_is(commands, "not") )
        {
            action = grid::off;
        }
        else if( next_command_is(commands,"and"))
        {
            action = grid::on;
        }
        else
        {
            break;
        }
    }
}

static void mark_grid( grid &g, commandlist &commands, string markcommand, grid &refgrid )
{
    mark_grid( g, commands, markcommand, &refgrid );
}

static void run_write_grid( grid &g, commandlist &commands )
{
    bool csv = false;
    const char *crlf="\n";
    string filename;
    string columns = "*";
    int ndp=-1;
    while( 1 )
    {
        filename = next_command(commands,"Filename for write operation");
        if( filename == "csv" ) 
        {
            csv=true;
        }
        else if( filename == "dos" ) 
        {
            crlf="\r\n";
        }
        else if( filename == "columns" )
        {
            columns=next_command(commands,"Grid columns to write");
            if( columns == "none" ) columns = "";
        }
        else if( filename == "ndp" )
        {
            next_command_value(commands,ndp,"Precision (ndp) of writing grid");
        }
        else
        {
            break;
        }
    }
    bool marked=false;
    if( next_command_is(commands,"where") )
    {
        mark_grid( g, commands, "write" );
        marked=true;
    }
    vector<int> colids;
    grid_columns( g, columns, colids );
    cout << "Writing grid" << (marked ? " data" : "") << " to file " << filename << endl;
    g.writefile(filename.c_str(),csv ? "," : "\t",crlf,&colids,marked,ndp);
}

static void run_precision( grid &g, commandlist &commands )
{
    int prec;
    next_command_value(commands,prec,"Precision for writing grid data");
    g.setprecision(prec);
}


static void run_write_linzgrid( grid &g, commandlist &commands )
{
    string crdsys=next_command(commands,"Coordinate system for write_linzgrid operation");
    string header1=next_command(commands,"First header line for write_linzgrid operation");
    string header2;
    string header3;
    if( header1 == "file" )
    {
        string headerfile=next_command(commands,"File name for linzgrid headers");
        ifstream is(headerfile.c_str());
        {
            if( is )
            {
                if( is ) getline(is,header1);
                if( is ) getline(is,header2);
                if( is ) getline(is,header3);
                is.close();
            }
            else
            {
                throw runtime_error(string("Cannot open linzgrid header file ")+header1);
            }
        }
    }
    else
    {
        header2=next_command(commands,"Second header line for write_linzgrid operation");
        header3=next_command(commands,"Third header line for write_linzgrid operation");
    }
    string linzgridfile=next_command(commands,"File name for write_linzgrid operation");
    string vres = "AUTO";
    string columns = "*";
    while( true )
    {
        if( linzgridfile == "resolution" )
        {
            vres = next_command(commands,"Resolution of grid data for write_linzgrid operation");
        }
        else if( linzgridfile == "columns" )
        {
            columns = next_command(commands,"Columns of grid data to write to linzgrid file");
        }
        else 
        {
            break;
        }
        linzgridfile=next_command(commands,"File name for write_linzgrid operation");
    }
    cout << "Writing grid to LINZ ASCII format in " << linzgridfile << endl;
    write_linz_grid_file(g,crdsys.c_str(),header1.c_str(),header2.c_str(),header3.c_str(),vres,columns,linzgridfile.c_str());
}

static void run_stats( grid &g )
{
    print_grid_stats( g );
}

static void apply_zero( grid &g, grid::node &n, void *data )
{
    vector<double>::pointer p = g.values(n);
    for( int i = 0; i < g.nvalue(); i++, p++ )
    {
        *p = 0.0;
    }
}

static void run_zero_grid( grid &g, commandlist &commands )
{
    mark_grid(g,commands,"zero operation");
    cout << "Zeroing grid at " << g.markCount() << " points" << endl;
    g.processMarked( apply_zero );
}

static void run_affected_area( grid &g, commandlist &commands )
{
    mark_grid(g,commands,"affected area operation");
    string wktfile = next_command(commands,"WKT output file name for affected area operation");

    cout << "Writing outline of cells affected by " << g.markCount() << " select nodes to " << wktfile << endl;
    outline_affected_cells(g,wktfile.c_str());
}

static void run_smoothgrid( grid &g, commandlist &commands ) 
{
    string mode = next_command(commands,"Mode for smooth operation");
    if( mode != "linear" && mode != "quadratic" && mode != "cubic" )
    {
        throw runtime_error("Invalid smoothing mode selected on command line: "+mode);
    }
    mark_grid(g,commands,"smooth operation");
    cout << "Applying " << mode << " smoothing to grid to fill " << g.markCount() << " selected nodes" << endl;
    if( mode == "linear" ) smoothgrid(g,1);
    else if( mode == "quadratic" ) smoothgrid(g,2);
    else smoothgrid(g,3);
}

static void run_addgrid( grid &g, commandlist &commands, string &command )
{
    grid gadd;
    string filename = run_read_grid( gadd, commands, command );
    bool marked = false;
    if( next_command_is(commands,"where"))
    {
        mark_grid(g,commands,command);
        marked = true;
    }
    
    string action = 
        command == "add" ? "Adding" : 
        command == "subtract" ? "Subtracting" : "Replacing";
    cout << action << (marked ? " selected" : "")
        << " grid values from grid " << filename << endl;
    double factor0 = command == "replace" ? 0 : 1;
    double factor1 = command == "subtract" ? -1 : 1;
    g.add(gadd,factor0,factor1,marked);
}

static void run_aligntogrid( grid &g, commandlist &commands )
{
    grid galign;
    string filename = run_read_grid( galign, commands, "align" );
    cout <<  "Aligning grid to " << filename << endl;
    g.alignto(galign);
}

static void run_trimto( grid &g, commandlist &commands )
{
    int buffer=0;
    if( next_command_is(commands,"buffer"))
    {
        next_command_value(commands,buffer,"Buffer size for trimto");
    }
    if( next_command_is(commands,"wkt"))
    {
        string wktfile= next_command( commands, "File name for trimto wkt");
        cout <<  "Trimming to wkt " << wktfile;
        if( buffer ) cout << " with buffer " << buffer;
        cout << endl;
        double minx, maxx, miny, maxy;
        wkt_extents(wktfile.c_str(),minx,maxx,miny,maxy);
        g.trimto(minx,maxx,miny,maxy,buffer);
    }
    else
    {
        grid galign;
        string filename = run_read_grid( galign, commands, "align" );
        cout <<  "Trimming to grid " << filename;
        if( buffer ) cout << " with buffer " << buffer;
        cout << endl;
        g.trimto(galign,buffer);
    }
}

static void run_multiply( grid &g, commandlist &commands )
{
    string sfactor = next_command(commands,"Scale factor in multiply command");
    if( sfactor == "by" ) sfactor = next_command(commands,"Scale factor in multiply command");
    double factor;
    if( ! parse_string(sfactor, factor) )
    {
        ostringstream msg;
        msg << "Invalid multiplication factor " << sfactor << " in multiply command";
        throw runtime_error(msg.str());
    }
    cout << "Multiplying grid by " << factor << endl;
    g.multiplyBy( factor );
}

static void run_resize( grid &g, commandlist &commands )
{
    int rowmin, colmin, rowmax, colmax;
    bool relative=next_command_is(commands,"relative");
    next_command_value(commands,rowmin,"Resize minimum row");
    next_command_value(commands,colmin,"Resize minimum column");
    next_command_value(commands,rowmax,"Resize maximum row");
    next_command_value(commands,colmax,"Resize maximum column");
    if( relative ) { rowmax += g.nrow()-1; colmax += g.ncol()-1; }
    cout << "Resizing to include rows " << rowmin << " to " << rowmax
        << " and columns " << colmin << " to " << colmax << endl;
    g.resize(rowmin,colmin,rowmax,colmax);
}

static void run_trim( grid &g, commandlist &commands )
{
    int borderSize;
    next_command_value(commands,borderSize,"Trim border size");
    cout << "Trimming to leave " << borderSize << " zero rows/cols around border" << endl;
    g.trim(borderSize);
}

static void run_evaluate( grid &g, commandlist &commands )
{
    bool haveids=false;
    bool csv=false;
    string infile = next_command(commands,"Input file for evaluate");
    if( infile == "with_ids" ){ 
        haveids=true; 
        infile = next_command(commands,"Input file for evaluate");
        }
    if( infile == "csv" ){ 
        csv=true; 
        infile = next_command(commands,"Input file for evaluate");
        }
    if( infile == "at" ) infile = next_command(commands,"Input file for evaluate");
    string outfile = next_command(commands,"Output file for evaluate");
    if( outfile == "to" )  outfile = next_command(commands,"Output file for evaluate");

    cout << "Evaluating grid at points in " << infile << " - results in " << outfile << endl;
    ifstream fin(infile.c_str());
    if( ! fin ) throw runtime_error(string("Cannot open evaluation point input file ").append(infile));
    ofstream fout(outfile.c_str());
    if( !fout ) throw runtime_error(string("Cannot open evaluation point output file ").append(outfile)); 
    string input;
    string id;
    char delim = csv ? ',' : '\t';
    vector<double> v;
    fout << "lon" << delim << "lat";
    for( int i = 0; i < g.nvalue(); i++ )
    {
        fout << delim << g.fieldName(i);
    }
    fout << endl;

    fout << setprecision(12);
    while( getline(fin,input) )
    {
        if( csv ) std::replace(input.begin(),input.end(),',',' ');
        istringstream s(input);
        grid::point p;
        if( haveids ) s >> id;
        if( ! (s >> p.x >> p.y) ) continue;
        g.valueAt(p,v);
        if( haveids ) fout << id << "\t";

        fout << p.x << delim << p.y;
        for( vector<double>::iterator vi = v.begin(); vi != v.end(); vi++ )
        {
            fout << delim << (*vi);
        }
        fout << endl;
    }
}

static void run_extents_wkt( grid &g, commandlist &commands )
{
    string filename = next_command( commands, "File name for extents_wkt operation");
    cout << "Writing grid extents to " << filename << endl;
    write_grid_extents_wkt(g, filename.c_str() );
}

static void run_colstats( grid &g )
{
    double mean, max, min;
    for( int i = 0; i < g.nvalue(); i++ )
    {
        g.colstats(i,&mean,&min,&max);
        cout << "Column " << g.fieldName(i) 
            << " Mean: " << mean << "  Min: " << min << "  Max: " << max <<endl;
    }
}

static void run_command_file( grid &g, commandlist &commands )
{
    string filename = next_command(commands,"File name for run command");
    istream *f;
    if( filename == "-")
    {
        f = &cin;
    }
    else
    {
        f = new ifstream(filename.c_str());
        if( ! *f )
        {
            throw runtime_error(string("Cannot open command file ").append(filename));
        }
    }
    list<string> filecommands;
    string cmdline;
    while( getline(*f,cmdline) )
    {
        int p = cmdline.find_first_not_of(' ');
        if( p == string::npos || cmdline[p] == '#' ) continue;
        stringstream s(cmdline);
        string cmd;
        while( s >> cmd )
        {
            if( cmd[0] == '"' )
            {
                cmd = cmd.substr(1);
                char c;
                bool escape = false;
                while( s.get(c) )
                {
                    if( escape )
                    {
                        escape = false;
                    }
                    else
                    {
                        if( c == '\\' ) { escape = true; continue; }
                        if( c == '"' ) break;
                    }
                    cmd.append(1,c);
                }
            }
            filecommands.push_back(cmd);
        }
    }
    if( filename != "-" )
    {
        delete f;
    }
    commands.insert(commands.begin(),filecommands.begin(),filecommands.end());
    cout << "Running commands from " << filename << endl;
}
 
static void run_commands( commandlist &commands )
{
    grid g;
    try
    {

        while( commands.size())
        {
            string op = next_command(commands);
            lowercase(op);
            if( op == "zero" ) run_zero_grid(g,commands);
            else if( op == "smooth" ) run_smoothgrid(g,commands);
            else if( op == "precision" ) run_precision(g, commands );
            else if( op == "write" ) run_write_grid(g,commands);
            else if( op == "write_linzgrid" ) run_write_linzgrid(g,commands);
            else if( op == "read" ) run_read_grid(g,commands);
            else if( op == "create" ) run_create_grid(g,commands);
            else if( op == "run" ) run_command_file(g,commands);
            else if( op == "add" ) run_addgrid( g, commands, op );
            else if( op == "subtract" ) run_addgrid(g, commands, op );
            else if( op == "replace" ) run_addgrid(g, commands, op );
            else if( op == "alignto" ) run_aligntogrid( g, commands );
            else if( op == "trimto" ) run_trimto( g, commands );
            else if( op == "multiply" ) run_multiply( g, commands );
            else if( op == "evaluate" ) run_evaluate( g, commands );
            else if( op == "resize" ) run_resize( g, commands );
            else if( op == "stats" ) run_stats( g );
            else if( op == "trim" ) run_trim( g, commands );
            else if( op == "affected_area" ) run_affected_area( g, commands );
            else if( op == "extents_wkt" ) run_extents_wkt( g, commands );
            else if( op == "colstats" ) run_colstats( g );
            else
            {
                throw runtime_error("Invalid grid operation " + op );
            }
        }
    }
    catch( runtime_error e )
    {
        cout << "Error processing grid:" << endl << e.what() << endl;
    }
}
