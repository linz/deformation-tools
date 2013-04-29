// Make M_PI available to MS compiler
#define _USE_MATH_DEFINES
#include <math.h>
#define DTOR (M_PI/180.0)

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <list>
#include <stdlib.h>
#ifdef USE_BOOST_REGEX
#include <boost/regex.hpp>
#else
#include <regex>
#endif
#include "okada.h"
#include "get_image_path.h"

using namespace std;
#ifdef USE_BOOST_REGEX
using namespace boost;
#endif

// "Seismologist's projection" as Ian Reilly disparagingly called it.
//
// Choose (lat0, lon0) near the centre of the area of interest.
//
// y = (lat - lat0)*111.2 km
// x = (lon - lon0)/cos(lat0)*111.2 km
//
// I have never carefully investigated how quickly the errors grow with distance from (lat0, lon0) compared to a standard projection.  Usually, the geodetic data are sparse enough and the fault poorly enough known that such differences are unlikely to be important.
//
// For very big earthquakes (such as Sumatra) the half-space approximation breaks down and it is necessary to use  models of dislocations in a sphere.  Not the case for Darfield, fortunately.


class GNSProjection
{
public:
    GNSProjection()
        : lat0(0), lon0(0), isvalid(false), noproj(false) {}
    GNSProjection( double lon0, double lat0 )
        : isvalid(false), noproj(false)  { SetOrigin(lon0,lat0); }
    GNSProjection( const GNSProjection &proj )
        : lon0(proj.lon0), lat0(proj.lat0), clat0(proj.clat0), isvalid(proj.isvalid), noproj(proj.noproj) {}
    void SetOrigin( double lon0, double lat0 )
    {
        this->lon0 = lon0; this->lat0 = lat0; this->clat0=cos(DTOR*lat0); noproj = false; isvalid = true;
    }
    void SetNoProj( bool np = true ) { noproj = np; isvalid=true; }
    void XY( double lon, double lat, double &x, double &y );
    void LatLon( double x, double y, double &lon, double &lat );
    bool IsValid() { return isvalid; }
private:
    void setup(double lon0, double lat0);
    static const double scale;
    double lat0;
    double lon0;
    double clat0;
    bool isvalid;
    bool noproj;
};

const double GNSProjection::scale = 111200.0;

void GNSProjection::XY( double lon, double lat, double &x, double &y )
{
    if( noproj) { x = lon; y = lat; return; }

    y = (lat - lat0)*scale;
    x = (lon - lon0)*scale*clat0;
};

void GNSProjection::LatLon( double x, double y, double &lon, double &lat )
{
    if( noproj ) { lon = x; lat = y; return; }
    lat = y/scale + lat0;
    lon = x/(scale * cos(DTOR*lat)) + lon0;
};

class FaultSet
{
public:
    FaultSet() {};
    FaultSet(GNSProjection proj) : proj(proj) {};
    ~FaultSet();
    bool ReadGNSDefinition( istream &str, int nskip = 0 );
    bool CalcOkada( double lon, double lat, double *dislocation, double *strain, bool reset=true );
    ostream & write( ostream &os, int style=0, bool header=true );
private:
    GNSProjection proj;
    list<SegmentedFault *>faults;
    list<string> names;
};

FaultSet::~FaultSet()
{
    for( list<SegmentedFault *>::const_iterator i = faults.begin(); i != faults.end(); i++ )
    {
        delete( *i );
    }
}


bool FaultSet::ReadGNSDefinition( istream &str, int nskip )
{
    string buffer;
    bool ok = true;
    bool have_num = false;
    bool have_type = false;
    bool have_open = false;
    bool started = false;
    regex re = regex("^(\\w+)\\:\\s+(.*?)\\s*$");

    while( ok )
    {
        getline(str,buffer);
        if( ! str.good()) break;
        if( nskip > 0 )
        {
            nskip--;
            continue;
        }

        smatch match;
        if( regex_match(buffer,match,re))
        {
            if(match[1].str() == "Origin" )
            {
                double lon, lat;
                stringstream s(match[2].str());
                s >> lon >> lat;
                if( s )
                {
                    proj.SetOrigin(lon,lat);
                }
            }
            continue;
        }

        // Trim comment and leading spaces
        int p = buffer.find('#');
        if( p != string::npos) { buffer.erase(p); }
        if( buffer.find_first_not_of(' ') == string::npos ) continue;

        // Check for column header, and discard
        if( buffer.find("fault_num") != string::npos ) have_num = true;
        if( buffer.find("fault_type") != string::npos ) have_type = true;
        if( buffer.find("opening") != string::npos ) have_open = true;

        stringstream s(buffer);

        int type, fnum;
        double strike,dip,rake,length,width,slip,Uts,depth,lat,lon;
        type = 3;
        fnum = 0;
        // s >> strike >> dip >> rake >> length >> width >> slip >> Uts >> depth >> lat >> lon;
        Uts = 0.0;
        if( have_num ) s >> fnum;
        if( have_type ) s >> type;
        s >> strike >> dip >> rake >> length >> width >> slip;
        if( have_open ) s >> Uts;
        s >> depth >> lat >> lon;
        if( !s )
        {
            if( started )
            {
                cerr << "Error reading fault definition\n"
                     << buffer << endl;
                ok = false;
            }
            continue;
        }
        started = true;

        double xf, yf;
        string name("");
        s >> xf >> yf;
        if( s ) getline(s,name);
        string::size_type pos = name.find_first_not_of(" \t");
        if( pos == string::npos )
        {
            name = "";
        }
        else
        {
            name = name.substr(pos);
        }
        // Convert to standard units
        length *= 1000.0;
        width *= 1000.0;
        depth *= 1000.0;
        strike = 90 - strike;
        dip = -dip;

        if( ! proj.IsValid() )
        {
            cout << setiosflags( ios::fixed ) << setprecision(5);
            cout << "Setting projection origin to " << lon << " " << lat << endl;
            proj.SetOrigin(lon,lat);
        }

        double x, y;
        proj.XY(lon,lat,x,y);

        double Uss = slip * cos(rake*DTOR);
        double Uds = slip * sin(rake*DTOR);

        double fs0=0, fs1=length;
        double fd0=0, fd1=width;

        switch( type )
        {
        case 0:
            fs0 = -fs0; fs1 = 0.0;
            fd0 = -fd1; fd1 = 0.0;
            break;
        case 1:
            fs0 = -fs0; fs1 = 0.0;
            fd0 = depth/sin(dip*DTOR);
            fd1 += fd0;
            depth = 0.0;
            break;
        case 2:
            fs0 = -fs0; fs1 = 0.0;
            break;
        case 3:
            fs1 /= 2; fs0 = -fs1;
            fd1 /= 2; fd0 = -fd1;
            break;
        case 4:
            fs1 /= 2; fs0 = -fs1;
            fd0 = depth/sin(dip*DTOR);
            fd1 += fd0;
            depth = 0.0;
            break;
        case 5:
            fs1 /= 2; fs0 = -fs1;
            break;
        }

        double fs[2] = { fs0, fs1 };
        double fd[2] = { fd0, fd1 };

        SegmentedFault *f = new SegmentedFault( x, y, depth, strike, dip );
        f->SetStrikeSegBreaks( 1, fs );
        f->SetDipSegBreaks( 1, fd );
        f->SetFaultSlip( 0, 0, Uss, Uds, Uts );

        faults.push_back( f );
        names.push_back(name);
    }
    if( faults.size() == 0 )
    {
        cerr << "Fault definition file contains no fault definitions" << endl;
        ok = false;
    }
    return ok;
}

bool FaultSet::CalcOkada( double lon, double lat, double *dislocation, double *strain, bool reset )
{
    if( reset )
    {
        dislocation[0] = dislocation[1] = dislocation[2] = 0.0;
        if( strain ) { strain[0] = strain[1] = strain[2] = strain[3] = 0.0; }
    }
    double x, y;
    bool ok = true;
    proj.XY(lon,lat,x,y);
    for( list<SegmentedFault *>::const_iterator i = faults.begin(); i != faults.end(); i++ )
    {
        if( ! (*i)->AddOkada(x,y,dislocation,strain)) ok = false;
    }
    return ok;
}

ostream &FaultSet::write( ostream &os, int style, bool header )
{
    bool havez = (style & 1) ? true : false;
    bool linestring = (style & 2) ? true : false;

    string prefix, zstring, suffix;
    zstring = havez ? " Z" : "";
    if( linestring ) { prefix="MULTILINESTRING"+zstring+"(("; suffix="))"; }
    else { prefix="POLYGON"+zstring+"(("; suffix="))"; }

    if( header )
        os << "id\tstrike\tdip\tUss\tUds\tUts\tmin_depth\tmax_depth\tname\tshape" << endl;

    int nflt = 0;

    // Note: this code is written for the case that every fault has just
    // one segment - GNS fault model supplied as multiple faults rather than
    // single segmented fault.

    list<string>::const_iterator in = names.begin();
    list<SegmentedFault *>::const_iterator i = faults.begin();
    for( in = names.begin(), i = faults.begin(); i != faults.end(); in++, i++ )
    {
        SegmentedFault *f = *i;
        string name = *in;
        replace(name.begin(),name.end(),'\t',' ');
        nflt++;

        os << setiosflags( ios::fixed )
           << setprecision(1)
           << nflt << "\t" << f->strike() << "\t" << f->dip() << "\t";

        double *sv = f->SlipVector(0,0);
        os << setprecision(4);
        os << sv[0] << "\t" << sv[1] << "\t" << sv[2] << "\t";
        os << setprecision(7);
        double x, y, z, minz, maxz;
        f->FaultLocation(0, 0, x, y, minz );
        f->FaultLocation(0, f->NDipSegments(), x, y, maxz );
        if( minz > maxz ) { z = minz; minz = maxz; maxz = z; }
        os << minz << "\t" << maxz << "\t";
        os << name << "\t";
        os << prefix;

        bool firstpt = true;
        for( int ipt = 0; ipt < 5; ipt++ )
        {
            int is = 0, id = 0;
            if( ipt == 2 || ipt == 3 ) id = f->NDipSegments();
            if( ipt == 1 || ipt == 2 ) is = f->NStrikeSegments();
            f->FaultLocation(is, id, x, y, z );
            z = -z;
            double lat, lon;
            proj.LatLon(x,  y, lon, lat );
            if( firstpt ) firstpt = false;
            else os << ", ";
            os << lon << " " << lat;
            if( havez ) os << " " << z;

        }
        // For line string representation show the surface projection
        // of the fault.
        if( linestring )
        {
            double x, y, z, lon, lat;
            os << "),(";
            f->FaultLocation(0, -1, x, y, z );
            z = -z;
            proj.LatLon(x,  y, lon, lat );
            os << lon << " " << lat;
            if( havez ) os << " " << z;
            f->FaultLocation(f->NStrikeSegments(), -1, x, y, z );
            z = -z;
            proj.LatLon(x,  y, lon, lat );
            os << ", " << lon << " " << lat;
            if( havez ) os << " " << z;
        }
        os << suffix << endl;
    }
    return os;
}

static void calcStrainComponents( double *strain, double &dil, double &rot, double &shear, double &err )
{
    double uxx = strain[0]*1.0e6;
    double uxy = strain[1]*1.0e6;
    double uyx = strain[2]*1.0e6;
    double uyy = strain[3]*1.0e6;

    dil=(uxx+uyy)/2.0;
    rot=(uyx-uxy)/2.0;
    shear=sqrt((uxx-uyy)*(uxx-uyy) + (uxy+uyx)*(uxy+uyx))/2;

    double ea = (uxx*uxx+uyx*uyx)/2;
    double eb = (uxy*uxy+uyy*uyy)/2;
    double ec = uxx*uxy+uyx*uyy;

    err = sqrt(fabs(ea+eb) + sqrt((ea-eb)*(ea-eb)+ec*ec));
}

void split_string( const string &source, const string &delim, list<string> &parts )
{
    string::size_type pos = 0;
    while( 1 )
    {
        string::size_type end= source.find(delim,pos);
        if( end == string::npos ) break;
        if( end > pos ) parts.push_back( source.substr(pos,end-pos) );
        pos = end+1;
    }
    parts.push_back( source.substr(pos) );
}

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

int main( int argc, char *argv[] )
{
    char *wktfile = 0;
    int wkttype = 3;
    bool isvalid = true;
    bool compare = false;
    bool havenames = false;
    bool showlength = false;
    bool calcstrain = false;
    int nskip = 0;
    int llprecision = 6;
    int dxyprecision = 4;
    int strnprecision = 4;
    GNSProjection proj;

    while( argc > 1 && argv[1][0] == '-' )
    {
        double lon, lat;
        switch( argv[1][1] )
        {
        case 'w':
        case 'W':
            if( argc > 2 )
            {
                wktfile = argv[2];
                if( argv[1][2] == 'p' || argv[1][2]=='P') wkttype=1;
                argv++;
                argc--;
            }
            else
            {
                cout << "-w requires wkt_file argument" << endl;
                isvalid = false;
            }
            break;
        case 'p':
        case 'P':
            if( argc > 2 && _stricmp(argv[2],"none") == 0 )
            {
                proj.SetNoProj();
                argv++;
                argc--;
            }
            else if( argc > 3
                     && sscanf(argv[2],"%lf",&lon) == 1
                     && sscanf(argv[3],"%lf",&lat) == 1  )
            {
                proj.SetOrigin(lon,lat);
                argv+=2;
                argc-=2;
            }
            else
            {
                cout << "-p requires lon and lat arguments or \"none\"" << endl;
                isvalid = false;
                break;
            }
            break;
        case 'n':
        case 'N':
            havenames = true;
            break;
        case 'l':
        case 'L':
            showlength = true;
            break;
        case 'c':
        case 'C':
            compare = true;
            break;
        case 's':
        case 'S':
            calcstrain = true;
            break;
        case 'h':
        case 'H':
            if( argc < 2 || sscanf(argv[2],"%d",&nskip) != 1 )
            {
                cout << "-h requires number of lines to skip" << endl;
                isvalid = false;
            }
            else
            {
                argv++;
                argc--;
            }
            break;
        case 'x':
        case'X':
            llprecision += 4;
            dxyprecision += 4;
            strnprecision += 4;
            break;
        default:
            cout << "Invalid argument " << argv[1] << endl;
            isvalid = false;
        }
        argv++;
        argc--;
    }

    if( ! isvalid )
    {
        cout << "Failed to run - invalid arguments" << endl;
        return 0;
    }

    if( argc < 2 )
    {
        help();
    }

    if( argc != 4)
    {
        cout << "Require parameters: [-w[l|p] wktfile] fault_model_file test_point_file output_file\n";
        return 0;
    }

    list<FaultSet *> faultlist;
    list<string> filenames;
    split_string(argv[1],"+",filenames);
    
    FaultSet faults(proj);

    // Read the fault model

    for( list<string>::iterator i = filenames.begin(); i != filenames.end(); i++ )
    {
        string filename(*i);
        if( filename.length() == 0 ) continue;
        ifstream f(filename.c_str());

        if( ! f.good() )
        {
            filename += ".model";
            f.open(filename.c_str());
        }

        if( ! f.good() )
        {
            const char *buf = getenv("FAULT_MODEL_DIR");
            if( buf ) 
            {
                string filename(buf);
                filename += "/";
                filename += argv[1];
                f.open(filename.c_str());
                if( ! f.good())
                {
                    filename += ".model";
                    f.open(filename.c_str());
                }
            }
        }

        if( ! f.good() )
        {
            cout << "Cannot read fault model file " << argv[1] << endl;
            return 0;
        }

        FaultSet *faults = new FaultSet( proj );
        if( ! faults->ReadGNSDefinition(f, nskip) ) { return 0; }
        faultlist.push_back(faults);
    }

    // Open the input and output files..
    
    istream *in;
    string fname(argv[2]);
    if( fname.substr(0,5) == "grid:")
    {
      replace(fname.begin(),fname.end(),':',' ');
      in = new istringstream(fname);
    }
    else
    {
      in = new ifstream(argv[2]);
      if( ! in->good() )
      {
	  cout << "Cannot open input test point file " << argv[2] << endl;
	  return 0;
      }
    }

    if( wktfile )
    {
        ofstream wkt(wktfile);
        if( ! wkt.good())
        {
            cout << "Cannot open wkt output file " << argv[4] << endl;
            return 0;
        }
        bool header = true;
        for( list<FaultSet *>::iterator f = faultlist.begin(); f != faultlist.end(); f++ )
        {
            (*f)->write(wkt,wkttype,header);
            header=false;
        }
        wkt.close();
    }

    double demax = 0.0;
    double dnmax = 0.0;
    double dumax = 0.0;

    ofstream out( argv[3]);
    if( ! out.good() )
    {
        cout << "Cannot open output file " << argv[3] << endl;
        return 0;
    }
    if( havenames ) out << "name\t";
    out << "lon\tlat\tde\tdn\tdu";
    if( showlength ) out << "\tds";
    if( calcstrain ) out << "\tdil\trot\tshear\terr";
    if( compare )
    {
        out << "\tobs_de\tobs_dn\tobs_du";
        if( showlength) out << "\tobs_ds";
        out << "\tdif_de\tdif_dn\tdif_du";
        if( showlength) out << "\tdif_ds";
    }
    out << endl;
    out << setiosflags( ios::fixed );
    string buffer;
    while( getline(*in, buffer) )
    {
        // Skip comments
        int p = buffer.find('#');
        if( p != string::npos) { buffer.erase(p); }
        if( buffer.find_first_not_of(' ') == string::npos ) continue;

        double lon,lat, uxyz[3], duxy[4];
        double *strain = calcstrain ? duxy : 0;
        double &ux = uxyz[0];
        double &uy = uxyz[1];
        double &uz = uxyz[2];

        double lon0, lat0, lon1, lat1, dlon, dlat;
        int nln, nlt;
        stringstream s(buffer);
        string gridstr;

        if( ! compare && ! havenames)
        {
            s >> gridstr >> lon0 >> lat0 >> lon1 >> lat1 >> nln >> nlt;
            if( s ) std::transform(gridstr.begin(),gridstr.end(),gridstr.begin(),::tolower);
        }
        if( ! s || gridstr != "grid" )
        {
            s.clear();
            s.seekg(0);
            string name;
            if( havenames ) s >> name;
            s >> lon0 >> lat0;
            double oux, ouy, ouz;
            if( compare )
            {
                s >> oux >> ouy >> ouz;
            }
            if( s.fail() )
            {
                cout << "Invalid test point: " << buffer << endl;
                continue;
            }
            bool reset = true;
            for( list<FaultSet *>::iterator f = faultlist.begin(); f != faultlist.end(); f++ )
            {
                (*f)->CalcOkada( lon0, lat0, uxyz, strain, reset );
                reset = false;
            }
            if( havenames ) out << name << "\t";
            out << setprecision(llprecision)
                << lon0 <<  "\t" << lat0
                << setprecision(dxyprecision)
                << "\t" << ux << "\t" << uy << "\t" << uz;
            if( showlength )
            {
                double us = sqrt(ux*ux+uy*uy);
                out << "\t" << us;
            }
            if( calcstrain )
            {
                double dil, rot, shear, err;
                calcStrainComponents(strain,dil,rot,shear,err);
                out << setprecision(strnprecision)
                    << "\t" << dil << "\t" << rot
                    << "\t" << shear << "\t" << err;
            }
            if( compare )
            {
                out << setprecision(dxyprecision)
                    << "\t" << oux << "\t" << ouy << "\t" << ouz;
                if( showlength )
                {
                    double us = sqrt(oux*oux+ouy*ouy);
                    out << "\t" << us;
                }
                oux -= ux; ouy -= uy; ouz -= uz;
                if( fabs(oux) > demax ) demax = fabs(oux); 
                if( fabs(ouy) > dnmax ) dnmax = fabs(ouy); 
                if( fabs(ouz) > dumax ) dumax = fabs(ouz); 
                out << "\t" << oux << "\t" << ouy << "\t" << ouz;
                if( showlength )
                {
                    double us = sqrt(oux*oux+ouy*ouy);
                    out << "\t" << us;
                }
            }
            out << endl;
        }
        else
        {
            dlon = (lon1 - lon0)/nln;
            dlat = (lat1 - lat0)/nlt;
            nln++; nlt++;
            for( lat = lat0; nlt > 0; nlt--, lat += dlat )
            {
                int nlnrow = nln;
                for( lon = lon0; nlnrow > 0; nlnrow--, lon += dlon )
                {
                    bool reset = true;
                    for( list<FaultSet *>::iterator f = faultlist.begin(); f != faultlist.end(); f++ )
                    {
                        (*f)->CalcOkada( lon, lat, uxyz, strain, reset );
                        reset=false;
                    }
                    out << setprecision(llprecision)
                        << lon <<  "\t" << lat << "\t"
                        << setprecision(dxyprecision)
                        << ux << "\t" << uy << "\t" << uz;
                    if( showlength )
                    {
                        double us = sqrt(ux*ux+uy*uy);
                        out << "\t" << us;
                    }
                    if( calcstrain )
                    {
                        double dil, rot, shear, err;
                        calcStrainComponents(strain,dil,rot,shear,err);
                        out << setprecision(strnprecision)
                            << "\t" << dil << "\t" << rot
                            << "\t" << shear << "\t" << err;
                    }
                    out << endl;
                }
            }
        }
    }
    delete(in);
    out.close();
    if( compare )
    {
        cout << "Max displacement differences (" << demax << ", " << dnmax << ", " << dumax << ")" << endl;
    }
    return 0;
}
