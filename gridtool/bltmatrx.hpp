#ifndef BLTMATRX_HPP
#define BLTMATRX_HPP

#include <assert.h>

// Header file for banded lower triangular matrix structure.
// Notes: All matrix operations are 1 based, not 0 based.

// These routines are very strongly focussed upon solving least squares
// equations.  In fact that's all they are good for!

class ProgressMeter;

class BLT_Def
{
    friend class BLT_Matrix;
    short nRows;
    short *rows;
public:
    BLT_Def( short nrows );
    ~BLT_Def();
    void SetNonZero( short row, short col );
    short Size()
    {
        return nRows;
    }
};

class BLT_Matrix
{
public:
    // Constructors for 1) full matrix
    //                  2) bandwidth limited matrix
    //                  3) matrix with specified first column non
    //                     zero for each row
    BLT_Matrix( short rows );
    BLT_Matrix( short rows, short bandwidth );
    BLT_Matrix( BLT_Def &def );
    ~BLT_Matrix();

    double & operator()( short i, short j );

    void Zero();
    short Size()
    {
        return nRows;
    }

    // Return 1 = OK, 0 = failure.  BadRow returns number of row in
    // which singularity is detected.
    short FormCholeskiDecomposition();
    // Note: Inverse doesn't include non-zero elements outside
    // original bandwidth - use with caution. It cannot be used
    // for solving equations.  This must be done from the Choleski
    // decomposition.
    short Invert();
    // Note: Solve can only be called when the Choleski decomposition
    // has been formed.
    short Solve( double *b, double *r );
    short BadRow()
    {
        return badrow;
    }
    long NonZeroCount();  // Total number of potentially n.z. elements in matrix
    static void SetProgressMeter( ProgressMeter *newMeter )
    {
        meter = newMeter;
    }

private:

    class BLT_Row
    {
        friend class BLT_Matrix;
        short first;
        short last;
        double *values;

        BLT_Row();
        ~BLT_Row();
        void SetSize( short first, short last );
        void Zero();
        double & operator()( short i );
        // double operator *( BLT_Row &r2 );
        // double operator *( double *vector );
    };

    short nRows;
    short badrow;

    BLT_Row *rows;
    BLT_Row &Row( short i );
    static double small;
    static double absSmall;
    static ProgressMeter *meter;

    enum { Entering, Choleski, Inverse, Garbage } status;
};

inline
double & BLT_Matrix::BLT_Row::operator ()( short i )
{
    assert( i >= first && i <= last );
    return values[i-first];
}

inline
BLT_Matrix::BLT_Row &BLT_Matrix::Row( short i )
{
    assert( i > 0 && i <= nRows );
    return rows[i-1];
}

inline
double & BLT_Matrix::operator ()( short i, short j )
{
    assert( i > 0 && i <= nRows );
    assert( j > 0 && j <= nRows );
    return ( i < j ) ? Row(j)(i) : Row(i)(j);
}

#endif
