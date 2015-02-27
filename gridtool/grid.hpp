#ifndef GRID_HPP
#define GRID_HPP

#include <string>
#include <vector>

class grid
{
public:
    enum markaction { on, off, toggle };
    class point 
    {
    public:
        point() : x(0.0), y(0.0){}
        point( double x, double y ) : x(x), y(y) {}
        double x;
        double y;
    };

    class node
    {
    public:
        node(int row = 0, int col = 0 ) : row(row), col(col) {};
        node( const node & l ): row(l.row), col(l.col) {}
        node & operator =(const node & l ){ row=l.row; col=l.col; return *this; }
        int row;
        int col;
    };

    template <class T> class vector2 : public std::vector<std::vector<T> >
    {
    public:
        vector2<T>(int rows=0, int cols=0 ) :
            std::vector<std::vector<T> >(rows)
        {
            for( int i = 0; i < rows; i++ )
            {
                this->at(i).resize(cols);
            }
        }
        vector2<T>(int rows, int cols, T initValue ) :
            std::vector<std::vector<T> >(rows)
        {
            assign( rows, cols, initValue );
        }

        void assign( int rows, int cols, T initValue )
        {
            this->resize(rows);
            for( int i = 0; i < rows; i++ )
            {
                this->at(i).assign(cols, initValue);
            }
        }

        typename std::vector<T>::reference operator()(int i, int j)
        {
            return (*this)[i][j];
        }
    };

    grid( const char * filename = 0, char delim=' ' );
    ~grid();
    bool readfile( const char * filename, char delim=' ', int maxcols=99 );
    bool writefile( const char * filename, const char *delim = 0, std::vector<int> *colids=0, bool markedonly=false );
    void setprecision( int dataprec ){ m_dataprec = dataprec; }
    int nrow() { return m_nrow; }
    int ncol() { return m_ncol; }
    int nvalue() { return m_nvalue; }
    std::vector<double>::pointer values( int row, int col )
        { return values(node(row,col));}
    std::vector<double>::pointer values( const node &n );
    int columnid( const std::string &colname );
    void colstats( int icol, double *mean, double *min, double *max );
    bool nearest( double x, double y, int &row, int &col );
    bool nearest( const point &p, node &n );
    bool nodexy( int row, int col, double &x, double &y );
    bool nodexy( const node &n, point &p );
    bool nodexy( const point &n, point &p );

    void clearMarked( markaction action = off );
    void mark( const node &n, markaction action=on );
    bool markNearest( const point &p, markaction action=on );
    void markEdge( int width, bool inside, markaction action=on);
    void toggleWithin( std::vector<point> &polygon, vector2<bool> &markBuffer );
    void markWhere( std::string field, std::string op, double value, markaction action=on );
    bool marked( node n ){ return marked(n.row, n.col ); }
    bool marked( int row, int col ){ return m_marked[row][col]; }
    void setMarked( int row, int col, bool value ){ m_marked[row][col] = value; }
    int markCount();
    void processMarked( void(*func)(grid &g, node &n, void *data), void *data = 0);
    bool valueAt( point &xy, std::vector<double> &values );
    void add( grid &g, double factor0=1.0, double factor1=1.0, bool markedonly=false );
    void alignto( grid &g );
    void trimto( grid &g );
    void multiplyBy( double factor );
    void resize( int rowmin, int colmin, int rowmax, int colmax );
    void create( double minx, double maxx, double incx,
            double miny, double maxy, double incy,
            std::vector<std::string> columns );
    void extents( double &minx, double &maxx, double &miny, double &maxy );
    void trim( int borderSize = 0 );
    std::string fieldName( int iv );
private:
    void initiallize();
    void create( int nrow, int ncol, int nvalue );
    void gridcoords( const point &wxy, point &gxy );
    bool valueIsZero( int row, int col );
    bool rowIsZero( int row );
    bool colIsZero( int col );
    void setupFields();

    std::string m_header;
    int m_nrow;
    int m_ncol;
    int m_nvalue;
    int m_dataprec;
    double m_x0;
    double m_y0;
    double m_coldx;
    double m_coldy;
    double m_collen;
    double m_rowdx;
    double m_rowdy;
    double m_rowlen;
    vector2<bool> m_marked;
    vector2<double> m_values;
    std::string m_lon;
    std::string m_lat;
    std::vector<std::string> m_fields;
};

#endif

