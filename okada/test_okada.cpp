#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "okada.h"

using namespace std;

bool ReadFaultModel( istream &str, SegmentedFault &f )
{
    enum { header, rows, cols, dislocations, done, error } status;
    string buffer;
    int nsegs, nsegd;
    int isegs=0, isegd=0;
    int nval, ival;
    double *splits;
    double slipvec[3];

    status = header;
    while( status != done && status != error)
    {
        getline(str,buffer);
        if( ! str.good()) break;


        // Trim comment and leading spaces
        int p = buffer.find('#');
        if( p != string::npos){ buffer.erase(p); }
        if( buffer.find_first_not_of(' ') == string::npos ) continue;

        stringstream s(buffer);

        switch( status )
        {
            case header:
                {
                    double x, y, d, strike, dip;
                    s >> x >> y >> d >> strike >> dip >> nsegs >> nsegd;
                    if( ! s )
                    {
                        cerr << "Failed to read fault header line: " << buffer << endl;
                        status = error;
                    }
                    else if( nsegs < 1 || nsegd < 1 || d < 0.0 )
                    {
                        cerr << "Invalid value in header line: ",buffer,"\n";
                        status = error;
                    }
                    else
                    {
                        f.SetFaultRefSys( FaultRefSys(x, y, d, strike, dip ));
                        status = rows;

                        nval = max(nsegs,nsegd);
                        if( nval > 0  ) splits = new double[nval+1];
                        nval = nsegs+1;
                        ival = 0;
                    }
                    
                    break;
                }
            case rows:
            case cols:
                {
                    while( ival < nval )
                    {
                        double d;
                        s >> d;
                        if( ! s )
                        {
                            if( s.eof()) break;
                            cerr << "Invalid value for row/column break: " << buffer << endl;
                            status = error;
                            break;
                        }
                        splits[ival] = d;
                        ival++;
                        if( ival == nval ) break;
                    }
                    if( ival == nval )
                    {
                        if( status == rows )
                        {
                            f.SetStrikeSegBreaks( nsegs, splits );
                            nval = nsegd+1;
                            ival = 0;
                            status = cols;
                        }
                        else
                        {
                            f.SetDipSegBreaks( nsegs, splits );
                            status = dislocations;
                            isegs = isegd = 0;
                        }
                    }
                    break;
                }
            case dislocations:
                {
                    double us, ud, ut;
                    s >> us >> ud >> ut;
                    if( ! s )
                    {
                        cerr << "Invalid value for dislocation: " << buffer << endl;
                        status = error;
                        break;
                    }
                    f.SetFaultSlip( isegs, isegd, us, ud, ut );
                    isegs++;
                    if( isegs == nsegs ) { isegs = 0; isegd++; }
                    if( isegd == nsegd ) { status = done; }
                }

        }
    }
    if( splits ) delete [] splits;
    return status == done;
}


ostream &writeFaultWkt( ostream &os, SegmentedFault &f, int style = 3 )
{
    bool havez = (style & 1) ? true : false;
    bool linestring = (style & 2) ? true : false;

    string prefix, zstring, suffix;
    zstring = havez ? " Z" : "";
    if( linestring ){ prefix="LINESTRING"+zstring+"("; suffix=")"; }
    else { prefix="POLYGON"+zstring+"(("; suffix="))"; }

        os << prefix
            << setiosflags( ios::fixed )
            << setprecision(7);

        bool firstpt = true;
        for( int ipt = 0; ipt < 5; ipt++ )
        {
            int is = 0, id = 0;
            if( ipt == 2 || ipt == 3 ) id = f.NDipSegments();
            if( ipt == 1 || ipt == 2 ) is = f.NStrikeSegments();
            double x, y, z, lon, lat;
            f.FaultLocation(is, id, x, y, z );
            z = -z;
            if( firstpt ) firstpt = false;
            else os << ", ";
            os << x << " " << y;
            if( havez ) os << " " << z;

        }
        os << suffix << endl;
    return os;
}


int main( int argc, char *argv[] )
{
    if( argc != 3 && argc != 4 && argc != 5)
    {
        cout << "Require parameters: fault_model_file [test_point_file] output_file [fault_wkt]\n";
        return 0;
    }
    SegmentedFault fault;

    // Read the fault model
    ifstream in(argv[1]);

    if( ! in.good() )
    {
        cout << "Cannot read fault model file " << argv[1] << endl;
        return 0;
    }

    if( ! ReadFaultModel( in, fault )){ return 0; }

    // Open the input and output files..
    
    char *ofile;
    if( argc == 3 )
    {
        ofile = argv[2];
    }
    else
    {
        ofile = argv[3];
        in.close();
        in.open(argv[2]);
        if( ! in.good() )
        {
            cout << "Cannot open input test point file " << argv[2] << endl;
            return 0;
        }
    }
    if( argc == 5 )
    {
       ofstream wkt( argv[4]);
       if( ! wkt.good() )
       {
           cout << "Cannot open wkt file " << argv[4] << endl;
           return 0;
       }
       writeFaultWkt( wkt, fault );
       wkt.close();
    }
    ofstream out( ofile );
    if( ! out.good() )
    {
        cout << "Cannot open output file " << ofile << endl;
        return 0;
    }
    out << "x\ty\tdx\tdy\tdz\tuxx\tuxy\tuyx\tuyy" << endl;
    string buffer;
    while( getline(in, buffer) )
    {
        int p = buffer.find('#');
        if( p != string::npos){ buffer.erase(p); }
        if( buffer.find_first_not_of(' ') == string::npos ) continue;

        double x, y, ux, uy, uz, uxx, uxy, uyx, uyy;
        stringstream s(buffer);
        s >> x >> y;
        if( ! s )
        {
            cout << "Invalid test point: " << buffer << endl;
            continue;
        }
        fault.OkadaDislocation( x, y, ux, uy, uz );
        fault.OkadaStrain( x, y, uxx, uxy, uyx, uyy );
        // uxx *= 1.0e6;
        // uxy *= 1.0e6;
        // uyx *= 1.0e6;
        // uyy *= 1.0e6;
        out << setiosflags( ios::fixed )
            << setprecision(2)
            << x <<  "\t" << y 
            << setprecision(5) << "\t"
            << ux << "\t" << uy << "\t" << uz << "\t"
            << uxx << "\t" << uxy << "\t" << uyx << "\t" << uyy
            << endl;
    }
    in.close();
    out.close();
    return 0;
}
