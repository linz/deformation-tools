
/*  This module processes banded format lower trianglar matrices for
    storage of symmetric matrices in least squares problems.  The routines
    provide for initiallization, allocation, Cholesky decomposition, and
    inversion of the matrix.  The inversion is not complete - only the
    banded components are inverted.

    The format of the banded matrix is...


       x
       x x
       0 x x
       0 0 0 x
       0 0 0 x x
       0 0 x x x x
       0 0 0 0 0 x x
       x x x x x x x x
       x x x x x x x x x

    The zero blocks are not stored.  This scheme is managed by a structure
    bltmatrix.  It consists of an array of doubles holding the non zero
    components row by row.  For each row two integers are held, a short
    integer identifies the column number of the first non-zero row for
    the array, and a long integer is a pointer into the array of doubles
    to identify the beginning of the row.

    The matrix is initiallized by creating the structure without the
    array of doubles, but with a specified number of row.  Also the number
    of the first non-sparse row is provided.  Then non-zero rows and columns
    are registered to determine their extent of the banding.  When all
    requirements have been registered, the double array and the pointers
    are initiallized.

    The initiallized matrix may then be used for summing and solving
    equations, inverting, etc.

    Note: all references to rows and columns in the matrix are 1
    based.  It is up to the calling program to ensure that calls are made
    in the correct order and that invalid rows and columns are not passed.

*/



#include <assert.h>
#include <math.h>

#include "bltmatrx.hpp"
#include "progress.hpp"

////////////////////////////////////////////////////////////////////////////

BLT_Def::BLT_Def( short nRows ) : nRows(nRows )
{
    rows = new short[nRows];
    for( short i=0; i<nRows; i++ ) rows[i] = i+1;
}

BLT_Def::~BLT_Def()
{
    if( rows ) delete [] rows;
}

void BLT_Def::SetNonZero( short row, short col )
{
    assert( row > 0 && row <= nRows );
    assert( col > 0 && col <= nRows );
    if( row < col )
    {
        if( row < rows[col-1] ) rows[col-1] = row;
    }
    else
    {
        if( col < rows[row-1] ) rows[row-1] = col;
    }
}

////////////////////////////////////////////////////////////////////////////

BLT_Matrix::BLT_Row::BLT_Row() : first(0), last(0), values(0) {};

BLT_Matrix::BLT_Row::~BLT_Row()
{
    if( values ) delete [] values;
}


// This implementation only allows the size to be set once...

void BLT_Matrix::BLT_Row::SetSize( short firstCol, short lastCol )
{
    assert( values == 0 );
    assert( last >= first );
    values = new double[lastCol - firstCol + 1];
    first = firstCol;
    last = lastCol;
}

void BLT_Matrix::BLT_Row::Zero()
{
    if( values )
    {
        double *v;
        short i;
        for( i = last - first + 1, v = values; i--; v++ )
        {
            *v = 0.0;
        }
    }
}

// double BLT_Matrix::BLT_Row::operator *( BLT_Matrix::BLT_Row &r2 ) {
//   short firstr = first > r2.first ? first : r2.first;
//   short lastr = last < r2.last ? list : r2.last;
//   if( firstr > lastr ) return 0.0;
//   double *v1 = values + firstr - first;
//   double *v2 = r2.values + firstr - r2.first;
//   double sum = 0.0;
//   while( firstr <= lastr ) { sum += *v1++ * *v2++; firstr++; }
//   return sum;
//   }

// double BLT_Matrix::BLT_Row::operator *( double *vector ) {
//   vector += first-1;
//   double sum = 0.0;
//   double *v = values;
//   for( short i = last-first+1; i--; ) {
//      sum += *vector++ * *v++;
//      }
//   return sum;
//   }


////////////////////////////////////////////////////////////////////////////

double BLT_Matrix::small = 1.0e-10;
double BLT_Matrix::absSmall = 1.0e-30;
ProgressMeter *BLT_Matrix::meter = 0;

BLT_Matrix::BLT_Matrix( short size )
{
    assert( size > 0 );
    nRows = size;
    rows = new BLT_Row[size];
    for( short i = 1; i <= nRows; i++ )
    {
        Row(i).SetSize( 1, i );
    }
    status = Entering;
    badrow = 0;
}

BLT_Matrix::BLT_Matrix( short size, short bandwidth )
{
    assert( size > 0 );
    assert( bandwidth > 0 );
    nRows = size;
    rows = new BLT_Row[size];
    for( short i = 1; i <= nRows; i++ )
    {
        short col = i - bandwidth + 1;
        if( col < 1 ) col = 1;
        Row(i).SetSize( col, i );
    }
    status = Entering;
    badrow = 0;
}

BLT_Matrix::BLT_Matrix( BLT_Def &def  )
{
    nRows = def.nRows;
    rows = new BLT_Row[nRows];
    for( short i = 1; i <= nRows; i++ )
    {
        Row(i).SetSize( def.rows[i-1], i );
    }
    status = Entering;
    badrow = 0;
}

BLT_Matrix::~BLT_Matrix()
{
    delete [] rows;
}


void BLT_Matrix::Zero()
{
    for( short i = 1; i <= nRows; i++ ) Row(i).Zero();
    status = Entering;
    badrow = 0;
}

long BLT_Matrix::NonZeroCount()
{
    long nelement = 0;
    for( short i = 1; i <= nRows; i++ )
    {
        BLT_Row &rowi = Row(i);
        nelement += rowi.last + 1 - rowi.first;
    }
    return nelement;
}

// Could be coded more cleanly using BLT_Row and defining operator *, but
// it is messy to manage leaving the last element of the sum and to retain
// the efficiencies of incrementing pointers.

short BLT_Matrix::FormCholeskiDecomposition()
{
    if( status == Choleski ) return 1;
    if( status != Entering ) return 0;
    if( meter )
    {
        meter->Start("Forming Choleski decomposition",NonZeroCount());
    }
    long ndone = 0;
    for( short i = 1; i <= nRows; i++ )
    {
        BLT_Row &rowi = Row(i);
        for( short j = rowi.first; j <= i; j++ )
        {
            BLT_Row &rowj = Row(j);
            double *vi = rowi.values;
            double *vj = rowj.values;
            short kmin = rowi.first;
            if( rowj.first < kmin )
            {
                vj += kmin - rowj.first;
            }
            else
            {
                kmin = rowj.first;
                vi += kmin - rowi.first;
            }
            double sum = 0.0;
            for( ; kmin < j; kmin++, vi++, vj++ )
            {
                sum -= *vi * *vj;
            }
            if( i == j )
            {
                sum += *vi;
                if( sum < *vi * small || sum < absSmall )
                {
                    status = Garbage;
                    badrow = i;
                    if( meter ) meter->Finish();
                    return 0;
                }
                *vi = sqrt( sum );
            }
            else
            {
                *vi = (*vi + sum) / *vj;
            }
            ndone++;
        }
        if( meter ) meter->Update( ndone );
    }
    if( meter ) meter->Finish();
    status = Choleski;
    return 1;
}

// The usual stuff, coded with lots of pointer increment, decrement, avoiding
// array indices at the expense of clarity.

short BLT_Matrix::Solve( double *b, double *s )
{
    if( status == Entering ) FormCholeskiDecomposition();
    if( status != Choleski ) return 0;
    double *bk = b;
    for ( short i=1; i <= nRows; i++, bk++ )
    {
        BLT_Row &rowi = Row(i);
        short kmin = rowi.first;
        double *sk = s + kmin-1;
        double *rk = rowi.values;
        double sum = 0.0;
        for( ; kmin < i; kmin++, rk++, sk++ ) sum -= *rk * *sk;
        *sk = (*bk + sum) / *rk;
    }

    double *si = s + nRows - 1;
    for (  int i = nRows; i; i--, si--)
    {
        double sum = 0.0;
        double *rk = s + nRows - 1;
        for (short k = nRows; k > i; k--, rk-- )
        {
            if( i >= Row(k).first ) sum -= Row(k)(i) * *rk;
        }
        *si = (*rk + sum)/Row(i)(i);
    }
    return 1;
}

/* Form the inverse within the banded structure using the
   cholesky decomposition

   NOTE: This is not a complete inverse... elements outside the
   original bandwidth should be non-zero.
*/

short BLT_Matrix::Invert()
{
    double sum;
    short i,j,k;

    if( status == Inverse ) return 1;
    if( status != Choleski ) FormCholeskiDecomposition();
    if( status != Choleski ) return 0;

    double *tmp = new double[nRows];

    if( meter )
    {
        meter->Start("Forming inverse",NonZeroCount());
    }

    long ndone = 0;
    for (i=nRows; i; i--)
    {

        /* Save a column of the decomposition */

        for (j=nRows; j >= i; j-- )
        {
            if( i >= Row(j).first )
            {
                tmp[j-1] = Row(j)(i);
            }
            else
            {
                tmp[j-1] = 0.0;
            }
        }

        for (j=nRows; j >= i; j-- )
        {

            if( i < Row(j).first ) continue;

            sum = (i==j) ? 1.0/tmp[i-1] : 0.0;

            for (k=i; k++ < nRows; )
            {
                if( k < j )
                {
                    if( k >= Row(j).first )
                        sum -= tmp[k-1]*Row(j)(k);
                }
                else
                {
                    if( j >= Row(k).first )
                        sum -= tmp[k-1]*Row(k)(j);
                }
            }

            Row(j)(i) = sum/tmp[i-1];
            ndone++;
        }
        if( meter ) meter->Update( ndone );
    }
    delete [] tmp;
    if( meter ) meter->Finish();
    status = Inverse;
    return 1;
}

