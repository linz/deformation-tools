#include "grid.hpp"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdexcept>

using namespace std;

grid::grid( const char * filename, char delim ) 
{
    initiallize();
    if( filename ) readfile( filename, delim );
}

grid::~grid()
{
}

void grid::initiallize()
{
    m_header="";
    m_lon="lon";
    m_lat="lat";
    m_nrow=0;
    m_ncol=0;
    m_nvalue=0;
    m_dataprec=-1;
    m_x0=0.0;
    m_y0=0.0;
    m_rowdx=0.0;
    m_rowdy=1.0;
    m_rowlen=1.0;
    m_coldx=1.0;
    m_coldy=0.0;
    m_collen=1.0;
    m_zerotol=0.0;
    m_values.resize(0);
    m_marked.resize(0);
    m_fields.resize(0);
}

void grid::create( int nrow, int ncol, int nvalue )
{
    if( nrow < 2 || ncol < 2 || nvalue < 0 )
    {
        stringstream msg;
        msg << "Invalid parameters to grid::create - nrow: " 
            << nrow << " ncol: " << ncol << " nvalue: " << nvalue;
        throw runtime_error(msg.str()) ;
    }

    m_nrow = nrow;
    m_ncol = ncol;
    m_nvalue = nvalue;
    if( nvalue == 0 ) 
    {
        m_values.resize(0);
    }
    else
    {
        m_values.assign(nrow,ncol*m_nvalue,0.0);
    }
    m_marked.resize(nrow);
    clearMarked();
}

void grid::create( double minx, double maxx, double incx,
            double miny, double maxy, double incy,
            std::vector<std::string> columns )
{
    int ngrdx=(int) floor((maxx-minx)/incx+0.99999);
    int ngrdy=(int) floor((maxy-miny)/incy+0.99999);

    ngrdx++;
    ngrdy++;
    create( ngrdy, ngrdx, columns.size());
    m_fields=columns;
    m_x0=minx;
    m_y0=miny;
    m_rowdx=0.0;
    m_rowdy=incy;
    m_coldx=incx;
    m_coldy=0.0;
    m_rowlen = m_rowdx*m_rowdx + m_rowdy*m_rowdy;
    m_collen = m_coldx*m_coldx + m_coldy*m_coldy;
    m_header="header"; // Kludge to ensure headers output
}

bool grid::readfile( const char * filename, char delim, int maxcols )
{
    initiallize();
    ifstream f(filename);

    if( ! f )
    {
        throw runtime_error(string("Invalid grid filename specified ") + filename );
    }
    m_filename=std::string(filename);
    delim;
    string buffer;
    string spacechar(" \n");
    int nvalue = 0;
    int npt = 0;
    int nrow = 0;
    int ncol = 0;
    double dx = 0, dy = 0;
    double x = 0, y = 0;
    double xr, yr;
    double xc, yc;

    while(getline(f,buffer))
    {
        if( ! f.good() ) break;
	std::replace(buffer.begin(), buffer.end(),delim,' ');
        stringstream s(buffer);
	

        double x0 = x, y0 = y, v1;

        if( ! (s >> x >> y) )
        {
            if( npt > 0 ) continue;
            if( buffer.find_first_not_of(spacechar) == string::npos ) continue;
            buffer.erase(0,buffer.find_first_not_of(spacechar));
            buffer.erase(buffer.find_last_not_of(spacechar)+1,buffer.length());
            m_header = buffer;
            continue;
        }
        npt++;

        if( npt == 1 )
        {
            nrow = 1;
            nvalue = 0;
            m_x0 = x;
            m_y0 = y;
            while( s >> v1 ) nvalue++;
            if( nvalue > maxcols ) nvalue=maxcols;
        }
        else if( npt == 2 )
        {
            dx = x-x0;
            dy = y-y0;
	    xc = x;
	    yc = y;
        }
        else if( (x-x0)*dx + (y-y0)*dy < 0.0 )
        {
            nrow++;
            xr = x;
            yr = y;
            if( ncol == 0 )
            {
                ncol = npt-1;
            }
            else if( (npt-1) % ncol != 0 )
            {
                stringstream message;
                message << "Grid data  does not form a regular grid ("
                        << npt << "," << ncol << ")" << endl << buffer;
                throw runtime_error( message.str());
            }
        }
        else if( nrow == 1 )
        {
            xc = x;
            yc = y;
        }
    }
    if( nrow < 2 || ncol < 2 || nrow * ncol != npt )
    {
        stringstream message;
        message << "Grid data does not form a regular grid ("
                << nrow << "," << ncol << "," << npt << ")";
        throw runtime_error(message.str());
    }

    m_rowdx = (xr - m_x0)/(nrow-1);
    m_rowdy = (yr - m_y0)/(nrow-1);
    m_coldx = (xc - m_x0)/(ncol-1);
    m_coldy = (yc - m_y0)/(ncol-1);

    m_rowlen = m_rowdx*m_rowdx + m_rowdy*m_rowdy;
    m_collen = m_coldx*m_coldx + m_coldy*m_coldy;

    if( fabs(m_rowdx*m_coldx + m_rowdy*m_coldy)/m_rowlen > 0.0000001 )
    {
        throw runtime_error("Grid axes are not orthogonal");
    }

    create( nrow, ncol, nvalue );

    int icol = -1;
    int irow = 0;
    f.clear();
    f.seekg(0);

    while(getline(f,buffer))
    {
        if( ! f.good() ) break;
	std::replace(buffer.begin(), buffer.end(),delim,' ');
        stringstream s(buffer);

        double x0 = x, y0 = y, v1;

        s >> x >> y >> v1;
        // Note: Currently don't test that points are actually on regular grid...
        // ... this is where it should happen.
        if( ! s ) continue;
        icol++;
        if( icol >= m_ncol ){ irow++; icol=0; }
        vector<double>::pointer pv = values( node(irow,icol) );

        *pv++ = v1;
        for( int iv = 1; iv < m_nvalue; iv++ )
        {
            s >> v1;
            if( ! s )
            {
                throw runtime_error("Missing data values in grid");
            }
            *pv++ = v1;
        }
    }
    if( irow != m_nrow-1 || icol != m_ncol-1 )
    {
        throw runtime_error("Failed to read all values from grid");
    }
    setupFields();
    return true;
}

void grid::extents( double &minx, double &maxx, double &miny, double &maxy )
{
    nodexy(0,0,minx,miny);
    nodexy(m_nrow-1,m_ncol-1,maxx,maxy);
    if( minx > maxx ){ double tmp=minx; minx=maxx; maxx=tmp;}
    if( miny > maxy ){ double tmp=miny; miny=maxy; maxy=tmp;}
}

void grid::setupFields()
{
    m_fields.assign(m_nvalue,string(""));
    istringstream s(m_header);
    string name;
    // Skip lon/lat fields
    s >> name;
    if( name != "" ) m_lon=name;
    s >> name;
    if( name != "" ) m_lat=name;
    for( int iv = 0; iv < m_nvalue; iv++ )
    {
        if( ! (s >> name) )
        {
            ostringstream ns;
            ns << "field_" << (iv+1);
            name = ns.str();
        }
        m_fields[iv] = name;
    }
}

string grid::fieldName( int iv )
{
    string name("");
    if( iv >= 0 && iv < nvalue())
    {
        if( iv >= m_fields.size() )
        {
            ostringstream ns;
            ns << "field_" << (iv+1);
            name = ns.str();
        }
        else
        {
            name = m_fields[iv];
        }
    }
    return name;
}

string grid::format_number(double value, int ndp )
{
    std::stringstream number;
    number.setf(std::ios_base::fixed,std::ios_base::floatfield);
    number << std::setprecision(ndp);
    number << value;
    std::string fnum=number.str();
    std::string::size_type flen=fnum.size();
    std::string::size_type ndot=fnum.find('.');
    if( ndot != std::string::npos )
    {
        ndot=ndot; // Want at least 1 digit after DP
        std::string::size_type endpos = fnum.find_last_not_of('0')+1;
        if( endpos == ndot+1 ) endpos=ndot;
        if( endpos < flen )
        {
            fnum=fnum.substr(0,endpos);

        }
    }
    return fnum;
}

bool grid::writefile( const char *filename, const char *delim, const char *crlf, std::vector<int> *colids , bool markedonly, int valuendp )
{
    if( ! delim ) delim = "\t";
    if( ! crlf ) crlf = "\n";

    ofstream f(filename, std::ios_base::out | std::ios_base::binary );
    m_filename=std::string(filename);
    int crdprec = 12;
    int dataprec = valuendp > 0 ? valuendp : m_dataprec > 0 ? m_dataprec : crdprec;
    double mult=0.5;
    for( int i=0; i<dataprec; i++ ) mult *= 10.0;


    std::vector<int> columns;
    int nvalues=0;
    if( colids )
    {
        nvalues=colids->size();
    }
    else
    {
        colids = &columns;
        nvalues=m_nvalue;
        for( int i = 0; i < m_nvalue; i++ )
        {
            columns.push_back(i);
        }
    }

    f.unsetf(std::ios_base::floatfield);   // f << setiosflags(ios::fixed);
    if( ! f )
    {
        throw runtime_error(string("Invalid grid filename specified ") + filename );
    }

    if( m_header != "" ) 
    {
        f << m_lon << delim << m_lat;
        for( int iv = 0; iv < nvalues; iv++ )
        {
            f << delim << m_fields[colids->at(iv)];
        }
        f << crlf;
    }

    for( int ir = 0; ir < m_nrow; ir++ )
    {
        int ipt0 = 0;
        for( int ic = 0; ic < m_ncol; ic++, ipt0+= m_nvalue )
        {
            if( markedonly && ! marked(ir,ic) ) continue;
            double x, y;
            nodexy(ir, ic, x, y );
            f << std::setprecision(crdprec);
            f << x << delim << y;
            for( int iv = 0; iv < nvalues; iv++ )
            {
                double value = m_values[ir][ipt0+colids->at(iv)];
                if( value == 0.0 )
                {
                    f << delim << "0";
                }
                else
                {
                    f << delim << format_number(value,dataprec);
                }
            }
            f << crlf;
        }
    }
    f.close();
    return true;
}

int grid::columnid( const std::string &name )
{
    int iv;
    for( iv = m_nvalue; iv--; ) { if( name == m_fields[iv] ) break; }
    return iv;
}

void grid::colstats( int icol, double *mean, double *min, double *max )
{
    *mean = 0;
    *min = 0;
    *max = 0;
    int npt = 0;
    for( int ir = 0; ir < m_nrow; ir++ )
    {
        int ipt = icol;
        for( int ic = 0; ic < m_ncol; ic++ )
        {
            double v = m_values[ir][ipt];
            *mean += v;
            if( npt == 0 ) *min = v; else if( v < *min ) *min = v;
            if( npt == 0 ) *max = v; else if( v > *max ) *max = v;
            npt++;
            ipt += m_nvalue;
        }
    }
    if( npt > 0 ) *mean /= npt;
}

void grid::gridcoords( const point &wxy, point &gxy )
{
    double x = wxy.x - m_x0;
    double y = wxy.y - m_y0;
    gxy.x = (x*m_coldx + y*m_coldy)/m_collen;
    gxy.y = (x*m_rowdx + y*m_rowdy)/m_rowlen; 
}

bool grid::nearest( double x, double y, int &row, int &col, double tolerance )
{
    node n;
    bool result = nearest(point(x,y), n, tolerance );
    row = n.row;
    col = n.col;
    return result;
}

bool grid::nearest( const point &p, node &n, double tolerance )
{
    point gxy;
    gridcoords(p,gxy);

    double gc=floor( gxy.x + 0.5 );
    double gr=floor( gxy.y + 0.5 );
    int ic = (int) gc;
    int ir = (int) gr;

    bool result = true;
    if( tolerance > 0.0 && 
            (fabs(gxy.x-ic) > tolerance || fabs(gxy.y-ir) > tolerance) )
    {
        result=false;
    }


    if( ic < 0 )
    {
        ic = 0;
        result = false;
    }
    else if( ic >= m_ncol )
    {
        ic = m_ncol;
        result = false;
    }

    if( ir < 0 )
    {
        ir = 0;
        result = false;
    }
    else if( ir >= m_nrow )
    {
        ir = m_nrow;
        result = false;
    }

    n.row = ir;
    n.col = ic;

    return result;
}

bool grid::nodexy( const point &n, point &p )
{
    p.x = m_x0 + m_coldx*n.x + m_rowdx*n.y;
    p.y = m_y0 + m_coldy*n.x + m_rowdy*n.y;
    return true;
}

bool grid::nodexy( int ir, int ic, double &x, double &y )
{
    x = m_x0 + m_coldx*ic + m_rowdx*ir;
    y = m_y0 + m_coldy*ic + m_rowdy*ir;
    return true;
}

bool grid::nodexy( const node &n, point &p )
{
    return nodexy( n.row, n.col, p.x, p.y );
}


std::vector<double>::pointer grid::values( const node &n )
{
    return &m_values[n.row][n.col*m_nvalue];
}

void grid::setValues( const node &n, const std::vector<double> &values )
{
    int iv0=n.col*m_nvalue;
    for( int iv=0; iv < m_nvalue; iv++ )
    {
        m_values[n.row][iv0+iv]=values[iv];
    }
}

bool grid::valueAt( point &xy, vector<double> &v )
{
    v.assign(m_nvalue,0.0);
    point gxy;
    gridcoords(xy,gxy);
    if( gxy.x < 0 || gxy.y < 0 || gxy.x > m_ncol-1 || gxy.y > m_nrow-1 ) return false;
    int r0 = (int)  floor(gxy.y);
    if( r0 == m_nrow-1) r0--;
    int c = (int)  floor(gxy.x);
    if( c == m_ncol-1) c--;
    int c1 = c+1;
    int r1 = r0+1, r;
    double rv0 = r1-gxy.y;
    double cv= c1-gxy.x;
    
    for( ; c <= c1; c++, cv=1-cv ) 
    {
        double rv = rv0;
        for (r = r0; r <= r1; r++, rv=1-rv )
        {
            vector<double>::pointer nodev = values(r,c);
            double f = cv*rv;
            for( int nv = 0; nv < m_nvalue; nv++ )
            {
                v[nv] += nodev[nv]*f;
            }
        }
    }
    return true;
}

void grid::clearMarked( markaction action)
{
    for( int row = 0; row < m_nrow; row++ )
    {
        switch( action )
        {
        case on:
            m_marked[row].assign(m_ncol,true);
            break;
        case off:
            m_marked[row].assign(m_ncol,false);
            break;
        case toggle:
            vector<bool> &rowv = m_marked[row];
            for( int col = 0; col < m_ncol; col++ )
            {
                rowv[col] = ! rowv[col];
            }
            break;
        }
    }
}

void grid::mark( const node &n, markaction action )
{
    switch(action)
    {
    case on: m_marked[n.row][n.col] = true; break;
    case off: m_marked[n.row][n.col] = false; break;
    case toggle :m_marked[n.row][n.col] = ! m_marked[n.row][n.col]; break; 
    }
}

int grid::markCount()
{
    int count = 0;
    for( int row = 0; row < m_nrow; row++ ) for( int col = 0; col < m_ncol; col++ )
    {
        if( marked(row,col) ) count++;
    }
    return count;
}

bool grid::markNearest( const point &p, markaction action, double tolerance )
{
    node n;
    bool result = nearest(p,n,tolerance);
    if( result ) mark(n,action);
    return result;
}

void grid::markEdge( int width, bool inside, markaction action )
{
    if( width < 0 ) return;
    if( width*2 >= m_nrow/2 || width*2 >= m_ncol  && inside ) return;
    if( width*2 >= m_nrow/2 || width*2 >= m_ncol )
    {
        for( int row = 0; row < m_nrow; row++ ) for( int col = 0; col < m_ncol; col++ )
        {
            mark(node(row,col),action);
            return;
        }
    }
    if( inside )
    {
        for( int row=width; row<m_nrow-width; row++)
            for( int col=width; col<m_ncol-width; col++)
                mark(node(row,col),action);
        return;
    }

    for( int row=0; row < width; row++ )
    {
        for( int col=0; col < m_ncol; col++ )
        {
            mark(node(row,col),action);
            mark(node(m_nrow-row-1,col),action);
        }
    } 
    for( int row=width; row < m_nrow-width; row++ )
    {
        for( int col=0; col < width; col++ )
        {
            mark(node(row,col),action);
            mark(node(row,m_ncol-col-1),action);
        }
    }
}

void grid::toggleWithin( vector<point> &polygon, vector2<bool> &markBuffer )
{
    int npt = polygon.size();
    if( npt < 3 ) return;

    // Start with the end point just in case the polygon is not closed
    point &p = polygon[npt-1];
    point gxy;
    gridcoords(p,gxy);
    double pcol0 = gxy.x;
    double prow0 = gxy.y;
    int row0 = (int) ceil(prow0);
    if( row0 < 0 ) row0 = 0;
    else if( row0 > m_nrow ) row0 = m_nrow;

    for( int ipt = 0; ipt < npt; ipt++ )
    {
        gridcoords(polygon[ipt],gxy);
        double pcol1 = gxy.x;
        double prow1 = gxy.y;
        int row1 = (int) ceil(prow1);
        if( row1 < 0 ) row1 = 0;
        else if( row1 > m_nrow ) row1 = m_nrow;

        if( row1 != row0 )
        {
            int r0, r1;
            if( row1 < row0 ){ r0=row1; r1=row0; }
            else { r0=row0;  r1=row1; }
            for( int r = r0; r < r1; r++ )
            {
                int c = (int)ceil((pcol0*(prow1-r)+pcol1*(r-prow0))/(prow1-prow0));
                if( c < 0 ) c = 0;
                while( c < m_ncol )
                {
                    markBuffer[r][c] = ! markBuffer[r][c];
                    c++;
                }
            }
        }
        prow0 = prow1;
        pcol0 = pcol1;
        row0 = row1;
    }
}


void grid::markWhere( std::vector<std::string> fields, std::string op, double value, markaction action )
{
    std::vector<int> colids;
    bool magnitude=false;
    if( fields.size() > 1 )
    {
        magnitude=true;
        value *= value;
    }
    for( int i=0; i < fields.size(); i++ )
    {
        std::string field=fields[i];
        if( field == "" ) continue;
        int iv=columnid(fields[i]);
        if( iv < 0 )
        {
            throw runtime_error(string("Invalid field ") + field + " referenced in grid::markWhere" );
        }
        colids.push_back(iv);
    }

    enum ops { LT, LE, EQ, GE, GT, NE } optype;
    if( op == "<" || op == "lt" ) optype = LT;
    else if ( op == "<=" || op == "le" ) optype = LE;
    else if ( op == "==" || op == "eq" ) optype = EQ;
    else if ( op == ">=" || op == "ge" ) optype = GE;
    else if ( op == ">" || op == "gt" ) optype = GT;
    else if ( op == "!=" || op == "ne" ) optype = NE;
    else
    {
        throw runtime_error(string("Invalid operation ") + op + " referenced in grid::markWhere" );
    }
    int nc=colids.size();
    for( int row = 0; row < m_nrow; row++ ) for( int col = 0; col < m_ncol; col++ )
    { 
        std::vector<double>::pointer vp=values(row,col);
        double v;
        if( magnitude )
        {
            v=0.0;
            for( int ic = 0; ic < nc; ic++ )
            {
                int colid=colids[ic];
                v += vp[colid]*vp[colid];
            }

        }
        else
        {
            v = vp[colids[0]];
        }
        bool match = false;
        switch (optype) 
        {
        case LT: match = v < value; break;
        case LE: match = v <= value; break;
        case EQ: match = v == value; break;
        case GE: match = v >= value; break;
        case GT: match = v > value; break;
        case NE: match = v != value; break;
        }
        if( match )
        {
            mark(node(row,col),action);
        }
    }

}

void grid::fillMarked()
{
    grid::vector2<bool> filled;
    filled.resize(m_nrow);
    for( int row = 0; row < m_nrow; row++ )
    {
        filled[row].assign(m_ncol,true);
        filled[row][0]=m_marked[row][0];
        filled[row][m_ncol-1]=m_marked[row][m_ncol-1];
    }
    filled[0]=m_marked[0];
    filled[m_nrow-1]=m_marked[m_nrow-1];

    bool finished=false;
    while( ! finished )
    {
        finished=true;
        for( int row = 0; row < m_nrow; row++ ) 
        {
            int col0 = 1;
            int col1 = m_ncol-2;
            for( ; col0 < m_ncol-1; col0++,col1-- )
            {
                if( filled[row][col0] && ! filled[row][col0-1] && ! m_marked[row][col0] )
                {
                    finished=false;
                    filled[row][col0]=false;
                }
                if( filled[row][col1] && ! filled[row][col1+1] && ! m_marked[row][col1] )
                {
                    finished=false;
                    filled[row][col1]=false;
                }
            }
        }
        for( int col = 0; col < m_ncol; col++ )
        {
            int row0 = 1;
            int row1 = m_nrow-2;
            for( ; row0 < m_nrow-1; row0++,row1-- )
            {
                if( filled[row0][col] && ! filled[row0-1][col] && ! m_marked[row0][col] )
                {
                    finished=false;
                    filled[row0][col]=false;
                }
                if( filled[row1][col] && ! filled[row1+1][col] && ! m_marked[row1][col] )
                {
                    finished=false;
                    filled[row1][col]=false;
                }
            }
        }
    }
    for( int row = 0; row < m_nrow; row++ )
    {
        m_marked[row]=filled[row];
    }
}

void grid::reverseMarked()
{
    for( int row = 0; row < m_nrow; row++ ) for( int col = 0; col < m_ncol; col++ )
    {
        m_marked[row][col] = ! m_marked[row][col];
    }
}

void grid::processMarked( void(*func)(grid &g, node &n, void *data), void *data)
{
    for( int row = 0; row < m_nrow; row++ ) for( int col = 0; col < m_ncol; col++ )
    {
        node n(row,col);
        if( m_marked[row][col] )(*func)(*this,n,data);
    }
}

void grid::multiplyBy( double factor, bool markedonly )
{
    for( int row = 0; row < m_nrow; row++ ) 
        for( int col = 0; col < m_ncol; col++ )
        {
            if( markedonly && ! marked(row,col) ) continue;
            vector<double>::pointer v = values(row,col);
            for( int iv = 0; iv <  m_nvalue; iv++ )
            {
                v[iv] *= factor;
            }
        }
}

void grid::add( double value, bool markedonly )
{
    for( int row = 0; row < m_nrow; row++ ) 
        for( int col = 0; col < m_ncol; col++ )
        {
            if( markedonly && ! marked(row,col) ) continue;
            vector<double>::pointer v = values(row,col);
            for( int iv = 0; iv <  m_nvalue; iv++ )
            {
                v[iv] += value;
            }
        }
}

void grid::add( grid &g, double factor0, double factor1, bool markedonly )
{
    int nv = nvalue();
    std::vector<int> colids;
    for( int i=0; i<m_nvalue; i++ )
    {
        int iv=g.columnid(fieldName(i));
        if( iv < 0 )
        {
            throw runtime_error("Incompatible grid values in grid::add");
        }
        colids.push_back(iv);
    }
    vector<double> gv;
    point xy;

    for( int row = 0; row < m_nrow; row++ ) 
        for( int col = 0; col < m_ncol; col++ )
        {
            if( markedonly && ! marked(row,col) ) continue;
            nodexy(row,col,xy.x,xy.y);
            if( g.valueAt(xy,gv) )
            {
                vector<double>::pointer tv = values(row,col);
                for( int iv = 0; iv < nv; iv++ )
                {
                    tv[iv] = factor0 * tv[iv] + factor1 * gv[colids[iv]];
                }
            }
        }
}

void grid::alignto( grid &g )
{
    // Determine the grid corners in terms of the alignment grid
    // Only really makes sense if grids have consistent x,y axes...
    int colmin=0;
    int colmax=m_ncol-1;
    int rowmin=0;
    int rowmax=m_nrow-1;

    point corner;
    point gcorner;
    point ncorner;

    nodexy(node(0,0),corner);
    g.gridcoords(corner,gcorner);
    g.nodexy(point(floor(gcorner.x+0.0001),floor(gcorner.y+0.0001)),corner);
    gridcoords(corner,ncorner);
    if(ncorner.x < -0.00001) colmin=int(floor(ncorner.x+0.0001));
    if(ncorner.y < -0.00001) rowmin=int(floor(ncorner.y+0.0001));

    nodexy(node(rowmax,colmax),corner);
    g.gridcoords(corner,gcorner);
    g.nodexy(point(ceil(gcorner.x-0.0001),ceil(gcorner.y-0.0001)),corner);
    gridcoords(corner,ncorner);
    if(colmax < ncorner.x-0.00001) colmax=int(ceil(ncorner.x-0.0001));
    if(rowmax < ncorner.y-0.00001) rowmax=int(ceil(ncorner.y-0.0001));

    resize(rowmin,colmin,rowmax,colmax);
}

void grid::expandto( grid &g )
{
    // Determine the grid corners in terms of the alignment grid
    // Only really makes sense if grids have consistent x,y axes...
    int colmin=0;
    int colmax=m_ncol-1;
    int rowmin=0;
    int rowmax=m_nrow-1;

    point corner;
    point gcorner;
    point ncorner;

    g.nodexy(point(0,0),corner);
    gridcoords(corner,ncorner);
    if(ncorner.x < -0.00001) colmin=int(floor(ncorner.x+0.0001));
    if(ncorner.y < -0.00001) rowmin=int(floor(ncorner.y+0.0001));

    g.nodexy(point(g.ncol()-1,g.nrow()-1),corner);
    gridcoords(corner,ncorner);
    if(colmax < ncorner.x-0.00001) colmax=int(ceil(ncorner.x-0.0001));
    if(rowmax < ncorner.y-0.00001) rowmax=int(ceil(ncorner.y-0.0001));

    resize(rowmin,colmin,rowmax,colmax);
}

void grid::trimto( grid &g, int buffer )
{
    double minx, maxx, miny, maxy;
    g.nodexy( 0, 0, minx, miny );
    g.nodexy( g.nrow()-1, g.ncol()-1, maxx, maxy );
    trimto( minx, maxx, miny, maxy, buffer );
}

void grid::trimto( double minx, double maxx, double miny, double maxy, int buffer )
{    
    int colmax=0;
    int colmin=m_ncol-1;
    int rowmax=0;
    int rowmin=m_nrow-1;

    point corner;

    for( int ix=0; ix<2; ix++ ) 
    {
        double gx = ix ? minx : maxx;
        for( int iy=0; iy<2; iy++ )
        {
            double gy=iy ? miny : maxy;
            gridcoords(point(gx,gy),corner);
            if( corner.y-0.0001 < rowmin ) rowmin=int(ceil(corner.y-0.0001));
            if( corner.y+0.0001 > rowmax ) rowmax=int(floor(corner.y+0.0001));
            if( corner.x-0.0001 < colmin ) colmin=int(ceil(corner.x-0.0001));
            if( corner.x+0.0001 > colmax ) colmax=int(floor(corner.x+0.0001));
        }
    }

    rowmin -= buffer;
    rowmax += buffer;
    colmin -= buffer;
    colmax += buffer;

    if( rowmin < 0 ) rowmin=0;
    if( rowmax >= m_nrow ) rowmax=m_nrow-1;
    if( colmin < 0 ) colmin=0;
    if( colmax >= m_ncol ) colmax=m_ncol-1;

    resize(rowmin,colmin,rowmax,colmax);

}

void grid::resize( int rowmin, int colmin, int rowmax, int colmax )
{
    if( rowmin >= rowmax || colmin >= colmax )
    {
        throw runtime_error("Invalid parameters for grid::resize");
    }
    point xy0;
    int nrow = (rowmax-rowmin)+1;
    int ncol = (colmax-colmin)+1;
    nodexy(rowmin,colmin,xy0.x,xy0.y);

    if( rowmax < m_nrow ) m_values.erase(m_values.begin()+rowmax+1,m_values.end());
    if( rowmin > 0 ) m_values.erase(m_values.begin(),m_values.begin()+rowmin);

    for( vector2<double>::iterator i = m_values.begin(); i != m_values.end(); i++ )
    {
        vector<double> &row = *i;
        if( colmax < m_ncol ) row.erase(row.begin()+(colmax+1)*m_nvalue, row.end());
        if( colmin > 0 ) row.erase(row.begin(), row.begin()+colmin*m_nvalue);
        if( colmin < 0 ) row.insert(row.begin(),(-colmin)*m_nvalue,0.0);
        if( colmax > m_ncol-1) row.insert(row.end(),(colmax-(m_ncol-1))*m_nvalue, 0.0 );
    }

    if( rowmin < 0 )
    {
        vector<double> newrow(ncol*m_nvalue,0.0);
        m_values.insert(m_values.begin(),(-rowmin),newrow);
    }

    if( rowmax > m_nrow-1 ) 
    {
        vector<double> newrow(ncol*m_nvalue,0.0);
        m_values.insert(m_values.end(),(rowmax-(m_nrow-1)),newrow );
    }

    m_marked.assign(nrow,ncol,false);
    m_nrow = nrow;
    m_ncol = ncol;
    m_x0 = xy0.x;
    m_y0 = xy0.y;
}

bool grid::valueIsZero( int row, int col, std::vector<int> &colids )
{
    vector<double>::pointer pv = values(row,col);
    for( vector<int>::iterator i = colids.begin(); i != colids.end(); i++ )
    {
        if( fabs(pv[*i]) > m_zerotol) return false;
    }
    return true;
}

bool grid::rowIsZero( int row, std::vector<int> &colids )
{
    for( int col = 0; col < m_ncol; col++ ) if( ! valueIsZero(row,col,colids) ) return false;
    return true;
}

bool grid::colIsZero( int col, std::vector<int> &colids )
{
    for( int row = 0; row < m_nrow; row++ ) if( ! valueIsZero(row,col,colids) ) return false;
    return true;
}

void grid::trim( int borderSize, double zerotol, std::vector<int> *colids, bool expand )
{
    int rowmin, rowmax, colmin, colmax;
    std::vector<int> tolcols;

    if( ! colids )
    {
        for( int i=0; i<m_nvalue; i++ ) tolcols.push_back(i);
        colids=&tolcols;
    }
    m_zerotol=zerotol;
    if( colids )
    for( rowmin = 0; rowmin < m_nrow; rowmin++ ) if( ! rowIsZero(rowmin,*colids)) break;
    if( rowmin == m_nrow )
    {
        throw runtime_error("grid::trim failed: all grid values are zero");
    }
    for( rowmax = m_nrow-1; rowmax >= 0; rowmax-- ) if( ! rowIsZero(rowmax,*colids)) break;

    for( colmin = 0; colmin < m_ncol; colmin++ ) if( ! colIsZero(colmin,*colids)) break;
    for( colmax = m_ncol-1; colmax >= 0; colmax--) if( ! colIsZero(colmax,*colids)) break;

    rowmin -= borderSize;
    rowmax += borderSize;
    colmin -= borderSize;
    colmax += borderSize;

    if( ! expand )
    {
        if( rowmin < 0 ) rowmin=0;
        if( rowmax > m_nrow-1 ) rowmax=m_nrow-1;
        if( colmin < 0 ) colmin=0;
        if( colmax > m_ncol-1 ) colmax=m_ncol-1;
    }

    if( rowmin >= rowmax || colmin >= colmax )
    {
        throw runtime_error("grid::trim failed: calculated bounds invalid");
    }

    resize(rowmin,colmin,rowmax,colmax);
}
