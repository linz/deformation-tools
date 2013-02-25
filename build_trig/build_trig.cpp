#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string>

#include "get_image_path.h"
#include "triangle.h"

using namespace std;

///////////////////////////////////////////////////////////////////

class triangle
{
public:
    triangle( int pt1, int pt2, int pt3 );
    int pt1;
    int pt2;
    int pt3;
    char used;
    triangle *next;
    static triangle *find(int pt1, int pt2, int pt3 );
    static void set_unused();
    static triangle *first;
private:
    int compare( triangle *cmp );
};

triangle *triangle::first = 0;

triangle::triangle( int pt1, int pt2, int pt3 )
{
    int tmp;
    if( pt1 > pt2 ) { tmp = pt1; pt1 = pt2; pt2 = tmp; }
    if( pt2 > pt3 ) { tmp = pt3; pt3 = pt2; pt2 = tmp; }
    if( pt1 > pt2 ) { tmp = pt1; pt1 = pt2; pt2 = tmp; }
    this->pt1 = pt1;
    this->pt2 = pt2;
    this->pt3 = pt3;
    this->used = 0;
    this->next = 0;
}

int triangle::compare( triangle *cmp )
{
    if( pt1 < cmp->pt1 ) return 1;
    if( pt1 > cmp->pt1 ) return -1;
    if( pt2 < cmp->pt2 ) return 1;
    if( pt2 > cmp->pt2 ) return -1;
    if( pt3 < cmp->pt3 ) return 1;
    if( pt3 > cmp->pt3 ) return -1;
    return 0;
}

triangle *triangle::find( int pt1, int pt2, int pt3 )
{
    triangle tmp(pt1,pt2,pt3);

    triangle **ptr = &first;
    while( *ptr )
    {
        int cmp = tmp.compare(*ptr);
        if( cmp == 0 ) return *ptr;
        if( cmp < 0 ) break;
        ptr = &((*ptr)->next);
    }

    triangle *trg = new triangle(pt1,pt2,pt3);
    trg->next = *ptr;
    (*ptr) = trg;
    return trg;
}

void triangle::set_unused()
{
    for( triangle *trg = first; trg; trg=trg->next ) trg->used = 0;
}

///////////////////////////////////////////////////////////////////

class datapoint
{
public:
    int id;
    double lat;
    double lon;
    double denu[3];
    double res;
    triangle *trg;
    static int nextid;
};

int datapoint::nextid = 1;

datapoint *pts;
int npts;
int npts2;

trg_handle trg = 0;

void dump_line( void *data, double x1, double y1, double x2, double y2 )
{
    FILE *trg_dump = (FILE *) data;
    fprintf(trg_dump,"1 %.5lf %.5lf\n0 %.5lf %.5lf\n",x1,y1,x2,y2);
}

void dump_triangles( FILE *out )
{
    for( triangle *t = triangle::first; t; t=t->next )
    {
        int pen = t->used ? 4 : 5;
        fprintf(out,"%d %.5lf %.5lf\n",pen,pts[t->pt1].lon,pts[t->pt1].lat);
        fprintf(out,"0 %.5lf %.5lf\n",pts[t->pt2].lon,pts[t->pt2].lat);
        fprintf(out,"0 %.5lf %.5lf\n",pts[t->pt3].lon,pts[t->pt3].lat);
    }
}

void dump_triangulation( char *filename, int trglist )
{
    FILE *trg_dump = fopen(filename,"w");
    if( trg_dump )
    {
        trg_plot(trg,TRG_PLOT_TRIANGLES,trg_dump, dump_line);
        for( int i = 0; i < npts; i++ )
        {
            fprintf(trg_dump,"2 %.5lf %.5lf\n",pts[i].lon,pts[i].lat);
        }
        for( int i = npts; i < npts2; i++ )
        {
            fprintf(trg_dump,"3 %.5lf %.5lf\n",pts[i].lon,pts[i].lat);
        }
        if( trglist ) dump_triangles(trg_dump);
        fclose(trg_dump);
    }
}

void die( char *message )
{
    printf("%s",message);
    // if( trg ) dump_triangulation("trg_dump.xy",1);
    exit(2);
}

double maxval3( double v1, double v2, double v3 )
{
    v1 = fabs(v1);
    v2 = fabs(v2);
    v3 = fabs(v3);
    if( v2 > v1 ) v1 = v2;
    if( v3 > v1 ) v1 = v3;
    return v1;
}


void add_trg_point( int npt )
{
    datapoint *pt = &(pts[npt]);
    int sts = trg_add_point( trg, pt->lon, pt->lat, npt );
    if( sts != TRG_ERR_OK )
    {
        printf("Error status %d adding point %d\n",sts,npt);
        die("Triangulation error\n");
    }
}

void load_points( FILE *f, int ndim, bool gns )
{
    char buf[256];
    double lon,lat,de,dn,du;
    int npt;

    // Count the number of points in the file ...

    fseek(f,0,SEEK_SET);
    npt = 0;
    while( fgets(buf,256,f) != NULL )
    {
        buf[255] = 0;
        if(sscanf(buf,"%lf%lf%lf%lf%lf",&lon,&lat,&de,&dn,&du) == 5) npt++;
    }

    // Allocate the buffer for the points

    pts = (datapoint *) malloc( sizeof(datapoint) * (npt+8));
    if( ! pts ) die("Cannot allocate space for points");

    // Load the points into the array

    fseek(f,0,SEEK_SET);
    npt = 0;
    double maxval;
    int nmaxval;
    double lonmin,lonmax,latmin,latmax;
    double hmult = 1.0;
    if( gns ) hmult = 0.001;
    double vmult = ndim != 2 ? hmult : 0.0;
    hmult = ndim != 1 ? hmult : 0.0;
    while( fgets(buf,256,f) != NULL )
    {
        buf[255] = 0;
        if(sscanf(buf,"%lf%lf%lf%lf%lf",&lon,&lat,&de,&dn,&du) != 5) continue;
        datapoint *pt = &(pts[npt]);
        npt++;
        pt->id = 0;
        if( ! gns )
        {
            pt->lat = lat;
            pt->lon = lon;
            pt->denu[0] = de*hmult;
            pt->denu[1] = dn*hmult;
        }
        else
        {
            pt->lat = lon;
            pt->lon = lat;
            pt->denu[0] = dn*hmult;
            pt->denu[1] = -de*hmult;
        }
        pt->denu[2] = du*vmult;
        pt->res = 0.0;
        pt->trg = 0;
        if( npt == 1 )
        {
            lonmin = lonmax = lon;
            latmin = latmax = lat;
            nmaxval = 0;
            maxval = maxval3( de, dn, du );
        }
        else
        {
            if( lon < lonmin ) lonmin = lon; else if (lon > lonmax ) lonmax = lon;
            if( lat < latmin ) latmin = lat; else if (lat > latmax ) latmax = lat;
            double v = maxval3( dn, de, du );
            if( v > maxval ) { maxval = v; nmaxval = npt-1; }
        }
    }

    // Add boundary points which will have zero deformation

    double lonmid = (lonmin + lonmax)/2;
    double dlon = lonmax-lonmid;
    double latmid = (latmin + latmax)/2;
    double dlat = latmax-latmid;

    // Set up the initial triangulation

    trg_init( &trg, lonmid-3*dlon, lonmid+3*dlon, latmid-3*dlat, latmid+3*dlat );
    if( ! trg ) die("Cannot create triangulation");
    add_trg_point( nmaxval );

    npts = npt;
    for( int i = -1; i <= 1; i++ )
    {
        for( int j = -1; j <= 1; j++ )
        {
            if( i == 0 && j == 0 ) continue;
            double factor = 1.5;
            if( i == 0 || j == 0 ) factor = 2.0;
            datapoint *pt = &(pts[npt]);
            npt++;
            pt->id = 0;
            pt->lon = lonmid+dlon*factor*i;
            pt->lat = latmid+dlat*factor*j;
            pt->denu[0] = 0.0;
            pt->denu[1] = 0.0;
            pt->denu[2] = 0.0;
            pt->res = 0.0;
            pt->trg = 0;
            add_trg_point(npt-1);
        }
    }

    npts2 = npt;
}

void find_triangle( int c1, int c2, int c3, void *data )
{
    triangle *t = triangle::find(c1,c2,c3);
    t->used = 1;
}

void reload_triangles()
{
    triangle::set_unused();
    trg_process_trg( trg, 0, find_triangle );
}

void calc_residual( int ipt )
{
    int id[3];
    double f[3],enu[3];
    datapoint *pt = &(pts[ipt]);

    int sts = trg_lin_interp( trg, pt->lon, pt->lat, id,f,id+1,f+1,id+2,f+2 );
    if( sts != TRG_ERR_OK )
    {
        printf("Error %d interpolating triangulation for point %d\n",sts,ipt);
        die("Interpolation error\n");
    }
    pt->trg = triangle::find(id[0],id[1],id[2]);
    pt->res = 0.0;
    if( id[0] == ipt || id[1] == ipt || id[2] == ipt ) return;
    enu[0] = enu[1] = enu[2] = 0;
    for( int i = 0; i < 3; i++ )
    {
        datapoint *p = &(pts[id[i]]);
        enu[0] += p->denu[0]*f[i];
        enu[1] += p->denu[1]*f[i];
        enu[2] += p->denu[2]*f[i];
    }
    for( int i = 0; i < 3; i++ )
    {
        double res = fabs(enu[i] - pt->denu[i]);
        if( res > pt->res ) pt->res = res;
    }
}

int get_max_residual()
{
    int nmax = -1;
    double maxres = 0.0;
    for( int i = 0; i < npts; i++ )
    {
        datapoint *pt = &(pts[i]);
        if( ! pt->trg || ! pt->trg->used ) calc_residual( i );
        if( pt->res > maxres ) { maxres = pt->res; nmax = i; }
    }
    return nmax;
}

void identify_triangle_nodes( int c1, int c2, int c3, void *data )
{
    if( ! pts[c1].id ) pts[c1].id = datapoint::nextid++;
    if( ! pts[c2].id ) pts[c2].id = datapoint::nextid++;
    if( ! pts[c3].id ) pts[c3].id = datapoint::nextid++;
}

void print_triangle_nodes( int c1, int c2, int c3, void *data )
{
    FILE *out = (FILE *) data;
    fprintf(out,"T %d %d %d\n",pts[c1].id,pts[c2].id,pts[c3].id);
}

void print_triangle_nodes_csv( int c1, int c2, int c3, void *data )
{
    FILE *out = (FILE *) data;
    fprintf(out,"%d,%d,%d\n",pts[c1].id,pts[c2].id,pts[c3].id);
}

void print_triangles( FILE *out, int ndim, bool csv=false )
{
    const char *sep=csv ? "," : " ";
    const char *prefix = csv ? "" : "P ";
    for( int i = 0; i < npts2; i++ ) { pts[i].trg = 0; pts[i].id = 0; }
    trg_process_trg( trg, 0, identify_triangle_nodes );
    datapoint::nextid=1;
    for( int i = 0; i < npts2; i++ ) { if( pts[i].id ) pts[i].id = datapoint::nextid++; }
    for( int i = 0; i < npts2; i++ )
    {
        datapoint *pt = &(pts[i]);
        if( pt->id )
        {
            fprintf(out,"%s%d%s%.5lf%s%.5lf%s%.4lf",prefix,pt->id,sep,pt->lon,sep,pt->lat,sep,pt->denu[0]);
            if( ndim >= 2 ) fprintf(out,"%s%.4lf",sep,pt->denu[1]);
            if( ndim == 3 ) fprintf(out,"%s%.4lf",sep,pt->denu[2]);
            fprintf(out,"\n");
        }
    }


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
    time_t timeval;
    char runtime[64];
    int ndim = 3;
    int valid = 1;
    bool gns  = false;
    bool verbose = false;
    bool csvoutput = false;

    time(&timeval);
    strftime(runtime,64,"%#d-%b-%Y %H:%M:%S",localtime(&timeval));

    double tol = 0.001;

    printf("build_trig: builds a triangulated path file from a grid of deformation values\n");

    while( valid && argc > 1 && argv[1][0] == '-' )
    {
        switch( argv[1][1] )
        {
        case 'h': case 'H': ndim = 2; break;

	case 'u': case 'U': ndim = 1; break;

        case 'g': case 'G': gns = true; break;

        case 'v': case 'V': verbose = true; break;

        case 'c': case 'C': csvoutput = true; break;

        case 't': case 'T':
            if( sscanf(argv[1]+2,":%lf",&tol) != 1 || tol <= 0.0 )
            {
                printf("Input tolerance %s invalid\n",argv[1]);
                valid = 0;
            }
            break;

        default: valid = 0; break;
        }
        argc--;
        argv++;
    }

    if( argc != 3 || ! valid )
    {
        help();
        printf("\nParameters: [-h] [-g] [-t:tolerance] [-v] input_file output_file\n");
        printf("Input file should contain one deformation value per line,\n");
        printf("with longitude,latitude,de,dn,du separated by spaces\n");
        printf("or if -g (GNS) is given latitude,longitude,ds,de,du\n");
        printf("The triangulation will be formed to fit all values to within the tolerance\n\n");
        printf("If -h is specified height changes are ignored\n");
	printf("If -u is specified horizontal changes are ignored\n");
	printf("If -c is specified CSV format output files are generated\n");
        printf("If -v is specified more output is generated\n");
        return 1;
    }

    FILE *in = fopen(argv[1],"r");
    if( ! in )
    {
        printf("Cannot open input file %s\n",argv[1]);
        return 1;
    }


    load_points(in,ndim,gns);

    printf("%d points loaded from %s\n",npts,argv[1]);

    while(1)
    {
        reload_triangles();
        int nmax = get_max_residual();
        double res = pts[nmax].res;
        if( res < tol ) break;
        add_trg_point(nmax);
        if( verbose )
        {
            printf("Adding point %d to triangulation - residual %.4lf\n",nmax,res);
        }
    }

    if( ! csvoutput )
    {
        FILE *out = fopen(argv[2],"w");
        if( ! out )
        {
            printf("Cannot open output file %s\n",argv[2]);
            return 1;
        }
        fprintf(out,"FORMAT TRIG2L\n");
        fprintf(out,"HEADER0 Data from %s\n",argv[1]);
        fprintf(out,"HEADER1 Created by build_trig program on %s\n",runtime);
        fprintf(out,"HEADER2 Residual tolerance %.5lf\n",tol);
        fprintf(out,"CRDSYS NZGD2000\n");
        fprintf(out,"NDIM %d\n",ndim);
        print_triangles(out,ndim);
        trg_process_trg( trg, out, print_triangle_nodes );
    }
    else
    {
        string s("trig_pts_");
        s += argv[2];
        s += ".csv";
        FILE *out = fopen(s.c_str(),"w");
        if( ! out )
        {
            printf("Cannot open output file %s\n",s.c_str());
            return 1;
        }
        fprintf(out,"id,lon,lat");
        if( ndim != 1 ) fprintf(out,",de,dn");
        if( ndim != 2 ) fprintf(out,",du");
        fprintf(out,"\n");
        print_triangles(out,ndim,true);
        fclose(out);
        s = "trig_trg_";
        s += argv[2];
        s += ".csv";
        out = fopen(s.c_str(),"w");
        if( ! out )
        {
            printf("Cannot open output file %s\n",s.c_str());
            return 1;
        }
        fprintf(out,"id1,id2,id3\n");
        trg_process_trg( trg, out, print_triangle_nodes_csv );
        fclose(out);
    }

    // dump_triangulation("trg_dump.xy",0);

    return 0;

}

