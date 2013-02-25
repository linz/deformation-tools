#include <vector>
#include <stdexcept>
using namespace std;

#include "smoothgrid.hpp"
#include "grid.hpp"
#include "gridutil.hpp"
#include "bltmatrx.hpp"

struct obseqcol
{
    int prm;
    int row;
    int col;
    double a;
    vector<double>::pointer values;
};

class obseq : public vector<obseqcol>
{
public:
    obseq(int ncol) : vector<obseqcol>(ncol) {};
};


static void set_xy_index( void *data, int r, int c )
{
    grid::vector2<int> &index = * (grid::vector2<int> *) data;
    index[r][c]=1;
}

static void row_obseq( vector<double> filter, int row, int col, obseq &obs )
{
    for( int i = 0; i < filter.size(); i++ )
    {
        obs[i].a = filter[i];
        obs[i].row = row + i;
        obs[i].col = col;
    }
}

static void col_obseq( vector<double> filter, int row, int col, obseq &obs )
{
    for( int i = 0; i < filter.size(); i++ )
    {
        obs[i].a = filter[i];
        obs[i].row = row;
        obs[i].col = col + i;
    }
}

static bool grid_obseq( grid &g, grid::vector2<int> &ids, obseq &obs)
{
    int nobs = obs.size();
    obseqcol &o = obs[nobs-1];
    if( o.row >= g.nrow() || o.col >= g.ncol() ) return false;
    bool ok = false;
    for( obseq::iterator oi = obs.begin(); oi != obs.end(); oi++ )
    {
        obseqcol &o = *oi;
        o.prm = ids[o.row][o.col];
        if( o.prm > 0 ) ok = true;
        o.values = g.values(o.row,o.col);
    }
    return ok;
}


static void smoothgrid( grid &g, vector<double> filter )
{
    int nrow = g.nrow();
    int ncol = g.ncol();
    int nval = g.nvalue();
    grid::vector2<int> ids(nrow,ncol);

    // Determine which elements are to be calculated

    int ncalc = 0;
    for( int i = 0; i < nrow; i++ ) for( int j = 0; j < ncol; j++ )
        {
            if( g.marked(i,j))
            {
                ncalc++;
                ids[i][j] = ncalc;
            }
        }

    if( ! ncalc )
    {
        // throw runtime_error("No points marked for smoothing" );
        return;
    }

    // Determine non-zero elements of normal equations

    BLT_Def size(ncalc);
    int nobscol = filter.size();
    obseq obs(nobscol);

    for( int r = 0; r < nrow; r++ ) for( int c = 0; c < ncol; c++ ) for( int mode = 0; mode < 2; mode++ )
            {
                switch(mode)
                {
                case 0:
                    row_obseq( filter, r, c, obs );
                    break;
                case 1:
                    col_obseq( filter, r, c, obs );
                    break;
                }
                if( ! grid_obseq( g, ids, obs ) ) continue;
                for( int i=1; i < nobscol; i++  ) for( int j = 0; j < i; j++ )
                    {
                        if( obs[i].prm && obs[j].prm ) size.SetNonZero(obs[i].prm,obs[j].prm);
                    }
            }

    // Form normal equations

    BLT_Matrix N(size);
	N.Zero();
    grid::vector2<double> b(nval,ncalc,0.0);
    vector<double> obsval;

    for( int r = 0; r < nrow; r++ ) for( int c = 0; c < ncol; c++ ) for( int mode = 0; mode < 2; mode++ )
            {
                switch(mode)
                {
                case 0:
                    row_obseq( filter, r, c, obs );
                    break;
                case 1:
                    col_obseq( filter, r, c, obs );
                    break;
                }
                if( ! grid_obseq( g, ids, obs ) ) continue;
                obsval.assign(nval,0.0);
                bool haveval = false;
                for( int i=0; i < nobscol; i++  )
                {
                    obseqcol &oi = obs[i];
                    if( oi.prm )
                    {
                        for( int j = 0; j <= i; j++ )
                        {
                            obseqcol &oj = obs[j];
                            if( oj.prm )
                            {
                                N(oi.prm, oj.prm ) += oi.a*oj.a;
                            }

                        }
                    }
                    else
                    {
                        haveval = true;
                        for( int j = 0; j < nval; j++ )
                        {
                            obsval[j] -= oi.values[j]*oi.a;
                        }
                    }
                }
                if( haveval )
                {
                    for( int i=0; i < nobscol; i++  )
                    {
                        obseqcol &oi = obs[i];
                        if( oi.prm )
                        {
                            for( int j = 0; j < nval; j++ )
                            {
                                b[j][oi.prm-1] += oi.a*obsval[j];
                            }
                        }
                    }
                }
            }

    // Solve the normal equations and put values back into grid...

    if( ! N.FormCholeskiDecomposition())
    {
        throw runtime_error("Failed to invert equations to smooth grid");
    }

    for( int j = 0; j < nval; j++ )
    {
        double *bj = &(b[j][0]);
        N.Solve(bj,bj);
    }

    // Install back in to the grid..

    for( int r = 0; r < nrow; r++ ) for( int c = 0; c < ncol; c++ )
        {
            int prm = ids[r][c];
            if( prm == 0 ) continue;
            prm--;
            vector<double>::pointer v = g.values(r,c);
            for( int j = 0; j < nval; j++, v++ )
            {
                *v = b[j][prm];
            }
        }
}


void smoothgrid( grid &g, int order  )
{
    vector<double> filter;
    switch (order)
    {
    case 1:
        filter.push_back(-1.0);
        filter.push_back(1.0);
        break;
    case 2:
        filter.push_back(-1.0);
        filter.push_back(2.0);
        filter.push_back(-1.0);
        break;
    case 3:
        filter.push_back(-1.0);
        filter.push_back(3.0);
        filter.push_back(-3.0);
        filter.push_back(1.0);
        break;
    default:
        throw runtime_error("Invalid order specified to grid::smoothgrid");
    }
    smoothgrid( g,filter );
}

