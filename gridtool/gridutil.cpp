#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstring>

using namespace std;

#include "math.h"
#include "grid.hpp"
#include "gridutil.hpp"

void mark_xy_file( grid &g, const char *xyfile, grid::markaction action )
{
    ifstream f(xyfile);
    if( ! f )
    {
        throw runtime_error(string("Cannot open XY point file ")+xyfile);
    }

    string buffer;
    while( getline(f,buffer) )
    {
        stringstream s(buffer);
        double x, y;
        if( s >> x >> y )
        {
            if( ! g.markNearest(grid::point(x,y), action) )
            {
                cout << "Invalid point in " << xyfile << ": " << x << " " << y << endl;
            }
        }
    }
    f.close();
}

void expand_marked( grid &g, double distance, bool expand )
{
    // Save the current state 
    grid::vector2<bool> copy(g.nrow(),g.ncol());

    for( int row=0; row < g.nrow(); row++ )
    {
        for( int col=0; col < g.ncol(); col++ )
        {
            copy(row,col) = g.marked(row,col) ^ expand;
        }
    }

    // Work out how far to expand
    if( distance < 0 ){ distance=-distance; expand=!expand; }
    distance += 0.0001;
    int maxdist = (int) distance;
    if( maxdist < 1 ) return;

    vector<int> offsets;
    offsets.assign(maxdist+1,0);
    for( int i=0; i <=maxdist; i++ )
    {
        double offset=sqrt(distance*distance-i*i);
        offsets[i] = (int) offset;
    }

    // Now apply to affected cells
    for( int row = 0; row < g.nrow(); row++ )
    {
        for( int col = 0; col < g.ncol(); col++ )
        {
            if( ! copy(row,col) ) continue;
            bool set = false;
            for( int i=0; i <= maxdist; i++ )
            {
                int row0 = max(row-i,0);
                int row1 = min(row+i,g.nrow()-1);
                for( int j=0; j <= offsets[i]; j++ )
                {
                    int col0 = max(col-j,0);
                    int col1 = min(col+j,g.ncol()-1);
                    if( copy(row0,col0) && copy(row0,col1) && copy(row1,col0) && copy(row1,col1)) continue;
                    set = true;
                    break;
                }
                if( set ) break;
            }
            if( set ) g.setMarked(row,col,expand);
        }
    }
}

void mark_wkt_file( grid &g, bool outside, const char *wkt_file, grid::markaction action  )
{
    grid::vector2<bool> inside;
    vector<grid::point> points;

    inside.resize(g.nrow());
    for( int row=0; row < g.nrow(); row++ )
    {
        inside[row].assign(g.ncol(),outside);
    }

    ifstream f(wkt_file);
    if( ! f )
    {
        throw runtime_error(string("Cannot open WKT file ")+wkt_file);
    }

    // Very crude loop to extract rings from WKT and toggle contents.
    // Assumes anything after a '(' is potentially a WKT ring.

    
    while( f )
    {
        char c;
        while(f.get(c)) { if( c == '(' ) break; }
        points.clear();
        if( ! f ) break;
        while( f )
        {
            grid::point pt;
            f >> pt.x >> pt.y;
            if (! f ) { if( ! f.eof() ) f.clear();  break; }
            while( f.get(c)) { if( c == ')' || c == ',' ) break; }
            points.push_back(pt);
        }
        if( points.size() > 2 ) g.toggleWithin(points,inside);
    }
    for( int row=0; row < g.nrow(); row++ )
    {
        for( int col=0; col < g.ncol(); col++ )
        {
            if( inside[row][col] ) g.mark(grid::node(row,col),action);
        }
    }
}

// Messy function to generate WKT polygon(s) outlining cells affected by marked nodes
// Not sorting out internal polygons - each ring is output as a separate wkt polygon.

void outline_affected_cells( grid &g, const char *wktfile )
{
    int nrow = g.nrow();
    int ncol = g.ncol();
    if( nrow < 2 || ncol < 2 )
    {
        throw runtime_error("gridutil::affected_cells: Cannot handle degenerate grid");
    }

    ofstream os(wktfile);
    if( ! os )
    {
        throw runtime_error(string("Cannot open output wkt file ").append(wktfile));
    }
    os << "wkt_geom" << endl;
    os << setprecision(12);

    // Affected grid cells, offset by 1, to allow a boundary around the edge of unaffected cells.
    grid::vector2<bool> affected(nrow+1,ncol+1,false);

    // Start off by identifing all affected cells into col_bdy

    for( int row = 0; row < nrow; row++ )
    {
        for( int col = 0; col < ncol; col++ )
        {
            if( g.marked(row,col))
            {
                affected(row,col) = true;
                affected(row,col+1) = true;
                affected(row+1,col) = true;
                affected(row+1,col+1) = true;
            }
        }
    }

    // Cells outside the grid cannot be affected.

    for( int row = 0; row <= nrow; row++ ) 
    {
        affected( row, 0 ) = false;
        affected( row, ncol ) = false;
    }
    for( int col = 0; col <= ncol; col++ )
    {
        affected( 0,col ) = false;
        affected(nrow,col) = false;
    }

    // Rows and columns that can be part of boundary .. 
    // row_bdy(i,j) is true if line from node(i,j) to (i,j+1) is
    // part of boundary. 
    // col_bdy(i,j) is true if line from node(i,j) to (i+1,j) is
    // part of boundary.

    grid::vector2<bool> row_bdy(nrow,ncol,false);
    grid::vector2<bool> col_bdy(nrow,ncol,false);

    // Compare adjacent cells to identify boundary elements

    for( int row = 0; row < nrow; row++ )
    {
        for( int col = 0; col < ncol; col++ )
        {
            row_bdy(row,col) = affected(row,col+1) ^ affected(row+1,col+1);
            col_bdy(row,col) = affected(row+1,col) ^ affected(row+1,col+1);
        }
    }

    // Search through the boundaries sequentially to find a starting point
    // Will delete edges as they are processed, so search sequentially for 
    // a starting point.

    for( int start_row = 0; start_row < nrow; start_row++ )
    {
        for( int start_col = 0; start_col < ncol; start_col++ )
        {
            if( ! row_bdy(start_row,start_col)) continue;
            int row = start_row;
            int col = start_col;
            double x, y;
            row_bdy(row,col) = false;
            os << "POLYGON((";
            g.nodexy(row,col,x,y);
            os << x << " " << y;
            col++;
            g.nodexy(row,col,x,y);
            os << "," << x << " " << y;

            while( true )
            {
                if( col > 0 && row_bdy(row,col-1)){ row_bdy(row,col-1) = false; col--;}
                else if( col < ncol-1 && row_bdy(row,col)){ row_bdy(row,col)=false; col++;}
                else if( row > 0 && col_bdy(row-1,col)){ col_bdy(row-1,col) = false; row--;}
                else if( row < nrow-1 && col_bdy(row,col)){ col_bdy(row,col) = false; row++;}
                else break;
                g.nodexy(row,col,x,y);
                os << "," << x << " " << y;
            }

            os << "))" << endl;
        }
    }
    os.close();
}

void write_grid_extents_wkt( grid &g, const char *wkt_file )
{
    ofstream os(wkt_file);
    if( ! os )
    {
        throw runtime_error(string("Cannot open output wkt file ").append(wkt_file));
    }
    os << "wkt_geom" << endl;
    os << setprecision(12);

    int nr = g.nrow();
    int nc = g.ncol();

    os << "POLYGON((";

    for( int i = 0; i < 5; i++ )
    {
        int row, col;
        switch(i)
        {
        case 0:
        case 4:
            row = 0; col = 0;
            break;
        case 1:
            row = nr; col = 0;
            break;
        case 2:
            row = nr; col = nc;
            break;
        case 3:
            row = 0; col = nc;
            break;
        }
        double x,y;
        g.nodexy(row,col,x,y);
        if(i) os << ", ";
        os << x << " " << y;
    }

    os << "))" << endl;
    os.close();

}

void grid_columns( grid &g, string columns, vector<int> &colids )
{
    colids.clear();

    if( columns == "*")
    {
        for( int i=0; i<g.nvalue(); i++ ) colids.push_back(i);
    }
    else
    {
        string::size_type pos = 0;
        string::size_type len = columns.size();
        while( pos < len )
        {
            string::size_type end=columns.find('+',pos);
            if( end == string::npos ) end=len;
            string col = columns.substr(pos,end-pos);
            pos=end+1;
            int id=g.columnid(col);
            if( id >= 0 )
            {
                colids.push_back(id);
            }
            else
            {
                throw runtime_error( string("Requested column invalid: ").append(col));
            }
        }
    }
}

void write_linz_grid_file( grid &g, string crdsys, string header1, string header2, string header3, string vres, string columns, const char *gridfile )
{
    int nrow = g.nrow();
    int ncol = g.ncol();

    int nval = g.nvalue();

    vector<int> useCols;
    grid_columns( g, columns, useCols );
    int noutval = useCols.size();

    if( ! noutval )
    {
        throw runtime_error("Cannot create LINZ grid as no output columns selected");
    }

    grid::point p0,p1x,p1y;
    g.nodexy(grid::node(0,0),p0);
    g.nodexy(grid::node(nrow-1,0),p1y);
    g.nodexy(grid::node(0,ncol-1),p1x);

    if( p0.x != p1y.x || p0.y != p1x.y )
    {
        throw runtime_error("Cannot create LINZ grid as grid not aligned with axes");
    }

    ofstream os( gridfile );
    if( ! os )
    {
        throw runtime_error(string("Cannot open output wkt file ").append(gridfile));
    }

    os << setprecision(12);
    os << "FORMAT: GRID2L" << endl;
    os << "HEADER0: " << header1 << endl;
    os << "HEADER1: " << header2 << endl;
    os << "HEADER2: " << header3 << endl;
    os << "CRDSYS: " << crdsys << endl;
    os << "NGRDX: " << ncol << endl;
    os << "NGRDY: " << nrow << endl;
    os << "XMIN: " << p0.x << endl;
    os << "XMAX: " << p1x.x << endl;
    os << "YMIN: " << p0.y << endl;
    os << "YMAX: " << p1y.y << endl;
    os << "VRES: " << vres << endl;
    os << "NDIM: " << noutval << endl;
    os << "LATLON: 1" << endl;
    os << "VALUES: REAL" << endl;

    for( int row = 0; row < nrow; row++ )
    {
        for( int col = 0; col < ncol; col++ )
        {
            os << "V" << col+1 << "," << row+1 << ":";
            vector<double>::pointer v = g.values(row,col);
            for( int i = 0; i < noutval; i++ )
            {
                os << " " << v[useCols[i]];
            }
            os << endl;
        }
    }
    os.close();
}


static double percentile( vector<double> sorted, double percentile )
{
    int size = sorted.size()-1;
    if( size < 0 ) return 0.0;
    if( size == 0 ) return sorted[0];
    double pv = size*percentile;
    if( pv < 0.0 ) pv = 0;
    if( pv > size ) pv = size;
    int v0 = (int) floor(pv);
    if( v0 == size ) v0--;
    int v1 = v0+1;
    return (pv-v0)*sorted[v1] + (v1-pv)*sorted[v0];
}

void print_grid_stats( grid &g )
{
    cout << "Grid statistics" << endl;
    cout << "Rows: " << g.nrow() << endl;
    cout << "Cols: " << g.ncol() << endl;
    cout << "Values per node: " << g.nvalue() << endl;
    if( g.nrow() == 0 || g.ncol() == 0 || g.nvalue() == 0 ) return;
    
    for( int iv = 0; iv < g.nvalue(); iv++ )
    {
            double sum = 0.0;
            int nv = 0;
            vector<double> v( g.nrow()*g.ncol());
            for( int ir = 0; ir < g.nrow(); ir++ )
            {
                for( int ic = 0; ic < g.ncol(); ic++ )
                {
                    double value =  *(g.values(ir,ic)+iv);
                    v[nv] = value;
                    sum += value;
                    nv++;
                }
            }
        sort(v.begin(),v.end());
        cout << "Column " << g.fieldName(iv) << ":" << endl;
        cout << setprecision(6);
        cout << "  Minimum: " << v[0] << endl;
        cout << "  Maximum: " << v[nv-1] << endl;
        cout << "  Mean: " << (sum/nv) << endl;
        cout << "  Median: " << percentile(v,0.5) << endl;
    }

}
