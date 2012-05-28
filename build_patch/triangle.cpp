/*=================================================================*/
/*                                                                 */
/*     Code to form Delauney triangles for a set of irregular      */
/*     points. Based on algorithm by Green and Sibson (1979),      */
/*     Computing Dirichelet Tesselations in the Plane, The         */
/*     Computer Journal, V 21, No 2, Pp 168-173.                   */
/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Encoded by Chris Crook,                                     */
/*     Department of Survey and Land Information,                  */
/*     January 1990                                                */
/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Note that we have assumed we will be handling relatively    */
/*     small numbers of points, so that little attempt has been    */
/*     made to encode the algorithm efficiently.  For example      */
/*     maintaining the edge list doubles its storage               */
/*     requirement, since for each edge we are storing two         */
/*     addresses. Also it could lead to excessive paging in        */
/*     virtual memory, as edge lists become disjoint in memory.    */
/*     Using arrays for these lists, as described by Green and     */
/*     Sibson, would be more efficient. Also for computational     */
/*     convenience we are retaining the endpoints of each line     */
/*     of the Dirichelet tesselation in the edge lists. This       */
/*     saves some computation time but does potentially increase   */
/*     the time of run for large numbers of points on virtual      */
/*     memory machines, since the larger data structure may lead   */
/*     to more page faults. Retaining these coordinates means      */
/*     that the handling of "constraints" (i.e. boundaries) is     */
/*     differently treated to that of Green and Sibson, and in     */
/*     this case only rectangular constraints are considered.      */
/*                                                                 */
/*     A warning on terminology - the word edge has been used in   */
/*     this code and comments for contiguity in Green and          */
/*     Sibsons paper.  A contiguity is a pair of points sharing    */
/*     a common edge                                               */
/*                                                                 */
/*=================================================================*/
/*                                                                 */
/*     This code has been written to support several simultaneous  */
/*     triangulations of an area - managed via a handle to         */
/*     a data block defining the current triangulation.  This      */
/*     capability is required for DCDB.                            */
/*                                                                 */
/*=================================================================*/

/* System functions required */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "triangle.h"

/*================================================================*/
/*                                                                */
/*   Basic structures required - A point structure, and an        */
/*   edge list structure.  The point structure includes a         */
/*   code, which may be used by the calling routine to            */
/*   identify the point.  These codes should always be            */
/*   positive - negative codes are used inside the routine to     */
/*   identify "constraints" on the problem (i.e. the bounds of    */
/*   the window within which the triangulation is performed.      */
/*                                                                */
/*   The edges of the structure are held in a closed linked       */
/*   list (a ring structure) pointed to by the list item of       */
/*   the point structure.                                         */
/*                                                                */
/*================================================================*/

struct edgerec;

struct pointrec
{
    struct pointrec *nextpt;        /* Next point in the linked list */
    double    x,y;                    /* Coordinates of the point */
    int     code;                   /* A numeric code for the point */
    struct edgerec *list;           /* Pointer to a list of edges around the point */
};

struct edgerec
{
    struct edgerec *next;           /* Pointer to the next edge */
    struct pointrec *pnt;           /* Pointer to the point beyond the edge */
    double endx, endy;                /* Coordinates of the end of the edge */
};

typedef struct pointrec *point;
typedef struct edgerec  *edgelist;

/*-----------------------------------------------------------------*/
/*                                                                 */
/*   The following structure provides a handle managing several    */
/*   simultaneous triangulations...                                */
/*                                                                 */
/*-----------------------------------------------------------------*/

typedef struct triangulation_rec
{
    point lstpt, schpt, bdypt;
    double  cxmn, cxmx, cymn, cymx, cxrf, cyrf;
    int   nopts, errsts;
} trg_rec;

static trg_rec *cur_trg;

#define lastpoint  cur_trg->lstpt
#define searchpt   cur_trg->schpt
#define boundarypt cur_trg->bdypt
#define cxmin      cur_trg->cxmn
#define cxmax      cur_trg->cxmx
#define cymin      cur_trg->cymn
#define cymax      cur_trg->cymx
#define cxref      cur_trg->cxrf
#define cyref      cur_trg->cyrf
#define no_points  cur_trg->nopts
#define errstat    cur_trg->errsts

/* Definitions for internal error conditions */

#define ERR_INIT  1234  /* Arbitrary error codes - INIT = initiallized, OK */
#define ERR_ABORT 4321  /* ABORT = error incurred, aborted */
#define ABORTED   (cur_trg == NULL || cur_trg->errsts != ERR_INIT)

/* Definitions for external error conditions */

#ifndef TRG_ERR_OK

#define TRG_FATAL_ERROR      20
#define TRG_VALID_TRIANGLE   10

#define TRG_ERR_OK           0
#define TRG_ERR_BOUNDS       1
#define TRG_ERR_DUPLICATE    2
#define TRG_ERR_NOPOINTS     3
#define TRG_ERR_OUTSIDE      4
#define TRG_ERR_NOTRIANGLE   (TRG_VALID_TRIANGLE + 1)

#define TRG_ERR_ABORTED      (TRG_FATAL_ERROR + 1)
#define TRG_ERR_ALLOC        (TRG_FATAL_ERROR + 2)
#define TRG_ERR_BUG          (TRG_FATAL_ERROR + 3)

#endif

/* Internally defined procedures required before their definition */


static void define_constraints( double xmin, double xmax, double ymin, double ymax );
static int insert_one_point( double x, double y, int code );
static int make_boundary( double width, int code );
static void free_up_memory();
static int find_triangle( double xt, double yt, point *p1, point *p2, point *p3 );

/*==========================================================================*/
/*                                                                          */
/*     EXTERNAL CALLS TO THE TRIANGULATION ROUTINES.                        */
/*                                                                          */
/*     The calls available are                                              */
/*                                                                          */
/*     trg_init  -  Initiallizes the triangulation routines and             */
/*                  defines the rectangular region within which             */
/*                  triangulation is to be performed.                       */
/*                  Parameters                                              */
/*                  handle  trg_handle*  Handle used to reference triangulation  */
/*                  wxmin   double    Extent of window within which           */
/*                  wxmax           triangulation is performed              */
/*                  wymin                                                   */
/*                  wymax                                                   */
/*                                                                          */
/*     trg_add_point - Adds a point to the existing triangulation.          */
/*                  Parameters                                              */
/*                  handle trg_handle  Handle used to reference triangulation    */
/*                  x,y    double   coordinates of the point                  */
/*                  code   int    a code for use by the calling routine     */
/*                                to identify the point                     */
/*                                                                          */
/*     trg_add_border - adds a border of points around the triangulation    */
/*                  (mainly for DCDB)                                       */
/*                  Parameters                                              */
/*                  handle  trg_handle Handle used to reference triangulation    */
/*                  width   double  the width of the border to add (>0)       */
/*                  code    int   the code used to identify border points   */
/*                                                                          */
/*     trg_triangle - Determines the vertices of the triangle enclosing     */
/*                  or nearest to a point.                                  */
/*                  handle  trg_handle   Handle used to reference triangulation  */
/*                  x,y     double    the coordinates fo the point            */
/*                  x1,y1   double *  coordinates of the first vertex         */
/*                  c1      int *   code of the first vertex                */
/*                  x2,y2   double *  coordinates of the second vertex        */
/*                  c2      int *   code of the second vertex               */
/*                  x3,y3   double *  coordinates of the third vertex         */
/*                  c3      int *   code of the third vertex                */
/*                                                                          */
/*     trg_lin_interp - Computes the parameters of a linear interpolation   */
/*                  across the triangle enclosing the point.                */
/*                  handle  trg_handle   Handle used to reference triangulation  */
/*                  x,y     double    the coordinates of the point            */
/*                  c1      int *   the code of the first vertex            */
/*                  m1      double *  the weighting of the first vertex       */
/*                  c2      int *   the code of the second vertex           */
/*                  m2      double *  the weighting of the second vertex      */
/*                  c3      int *   the code of the third vertex            */
/*                  m3      double *  the weighting of the third vertex       */
/*                                                                          */
/*     trg_sizes    Computes statistics of the lengths of triangle sides    */
/*                  handle  trg_handle   Handle used to reference triangulation  */
/*                  mins    double*   Minimum length of a side                */
/*                  maxs    double*   Maximum length of a side                */
/*                  means   double*   Mean length of a side                   */
/*                                                                          */
/*     trg_term     Terminates triangulation routines - frees up            */
/*                  allocated memory                                        */
/*                  handle  trg_handle*  Handle used to reference the trngl.     */
/*                                                                          */
/*     Additional calls for debugging purposes are:                         */
/*                                                                          */
/*     trg_dump -   Dumps a definition of the triangulation to a file       */
/*                  Parameter                                               */
/*                  handle   trg_handle   Handle used to reference triangulation */
/*                  out      FILE *  file to dump output to                 */
/*                                                                          */
/*     trg_plot -   Plots the Delauney triangles or the Dirichelet          */
/*                  tesselation                                             */
/*                  Parameter                                               */
/*                  handle   trg_handle   Handle used to reference triangulation */
/*                  mode      int  0 = plot tesselations, o.w. plot trngls  */
/*                  drawline  int()  function to draw a line from x1,y1     */
/*                                   to x2,y2 i.e. the function is          */
/*                                                                          */
/*                                   int drawline(x1,y1,x2,y2)              */
/*                                   double x1,y1,x2,y2;                    */
/*                                                                          */
/*     These routines all return integer error codes as follows             */
/*                                                                          */
/*     TRG_ERR_OK            No error                                       */
/*                                                                          */
/*     Warning errors                                                       */
/*                                                                          */
/*     TRG_ERR_BOUNDS        New node rejected - outside predefined         */
/*                           bounds.                                        */
/*     TRG_ERR_DUPLICATE     New node rejected - identical location to      */
/*                           an existing node                               */
/*     TRG_ERR_NOPOINTS      Call to create a boundary before any nodes     */
/*                           have been defined                              */
/*     TRG_ERR_OUTSIDE       The point referenced in trg_triangle is        */
/*                           outside the network of triangles               */
/*     TRG_ERR_NOTRIANGLE    Couldn't identify a triangle including or      */
/*                           nearest the point                              */
/*                                                                          */
/*     Fatal errors                                                         */
/*                                                                          */
/*     TRG_ERR_ABORTED       A previous call to the triangulation           */
/*                           procedures has incurred a fatal error          */
/*     TRG_ERR_ALLOC         Unable to allocate required memory             */
/*     TRG_ERR_BUG           An invalid state has been reached - due        */
/*                           to a program bug                               */
/*                                                                          */
/*     TRG_FATAL_ERROR       All error codes greater than this value        */
/*                           are fatal error conditions                     */
/*                                                                          */
/*==========================================================================*/


int trg_init( trg_handle *handle,
              double wxmin, double wxmax, double wymin, double wymax )
{
    cur_trg = (trg_rec *) malloc( sizeof(trg_rec) );
    if (cur_trg==NULL) return TRG_ERR_ALLOC;
    *handle = (trg_handle) cur_trg;
    errstat = ERR_INIT;
    lastpoint = NULL;
    searchpt = NULL;
    no_points = 0;
    cxref = (wxmin+wxmax)/2;
    cyref = (wymin+wymax)/2;
    define_constraints( wxmin-cxref, wxmax-cxref, wymin-cyref, wymax-cyref );
    return TRG_ERR_OK;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     trg_add_point may return an error code 1 if the point       */
/*     cannot be entered into the triangulation because it is      */
/*     outside the window defined in trg_init, or because it is    */
/*     coincident with an existing point in the triangulation      */
/*                                                                 */
/*-----------------------------------------------------------------*/


int trg_add_point( trg_handle handle, double x, double y, int code )
{
    int i;

    cur_trg = (trg_rec *) handle;
    if ( ABORTED ) return TRG_ERR_ABORTED;
    searchpt = NULL;
    i = insert_one_point( x-cxref, y-cyref, code );
    if (i>=TRG_FATAL_ERROR) errstat = ERR_ABORT;
    return i;
}


int trg_add_border( trg_handle handle, double width, int code )
{
    int i;

    cur_trg = (trg_rec *) handle;
    if ( ABORTED ) return TRG_ERR_ABORTED;
    i = make_boundary( width, code );
    if (i>=TRG_FATAL_ERROR) errstat = ERR_ABORT;
    return i;
}

int trg_term( trg_handle *handle )
{
    cur_trg = (trg_rec *) *handle;

    if (cur_trg==NULL) return 0;
    free_up_memory();
    free(cur_trg);
    *handle = NULL;
    return TRG_ERR_OK;
}



int trg_triangle( trg_handle handle,
                  double x,double y,
                  double *x1,
                  double *y1,
                  int *c1,
                  double *x2,
                  double *y2,
                  int *c2,
                  double *x3,
                  double *y3,
                  int *c3 )
{
    point p1,p2,p3;
    int sts;
    cur_trg = (trg_rec *) handle;
    if (ABORTED) return TRG_ERR_ABORTED;
    sts = find_triangle( x-cxref, y-cyref, &p1, &p2, &p3 );
    if (sts<=TRG_VALID_TRIANGLE)
    {
        *x1 = p1->x+cxref; *y1 = p1->y+cyref;  *c1 = p1->code;
        *x2 = p2->x+cxref; *y2 = p2->y+cyref;  *c2 = p2->code;
        *x3 = p3->x+cxref; *y3 = p3->y+cyref;  *c3 = p3->code;
    }
    else
    {
        *x1 = *y1 = *x2 = *y2 = *x3 = *y3 = 0.0;
        *c1 = *c2 = *c3 = 0;
    }
    return sts;
}

int trg_lin_interp( trg_handle handle, double x, double y,
                    int *c1, double *m1, int *c2, double *m2, int *c3, double *m3 )
{
    double x1,y1,x2,y2,x3,y3,a;
    int sts;
    sts = trg_triangle( handle, x, y, &x1,&y1,c1, &x2,&y2,c2, &x3,&y3,c3 );
    if (sts==TRG_ERR_OK)
    {
        a = (x1-x2)*(y3-y2)+(y1-y2)*(x2-x3);
        *m1 = ((x-x2)*(y3-y2) + (y-y2)*(x2-x3))/a;
        *m2 = ((x-x3)*(y1-y3) + (y-y3)*(x3-x1))/a;
        *m3 = ((x-x1)*(y2-y1) + (y-y1)*(x1-x2))/a;
    }
    else *m1 = *m2 = *m3 = 0.0;
    return sts;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Find the sizes of typical triangles in the triangulation    */
/*     (min, max, and mean lengths of triangles)                   */
/*                                                                 */
/*-----------------------------------------------------------------*/

int trg_sizes( trg_handle handle, double *mins, double *maxs, double *means )
{
    point pt;
    edgelist lst;
    double dx,dy,ds;
    int nside;

    cur_trg = (trg_rec *) handle;
    if(ABORTED) return TRG_ERR_ABORTED;
    if(no_points<2) return TRG_ERR_NOTRIANGLE;

    *mins = *maxs = *means = 0.0;
    nside=0;
    for( pt=lastpoint; pt != NULL; pt = pt->nextpt ) if (pt->code >= 0)
        {
            lst = pt->list;
            if (lst != NULL) do
                {
                    if (lst->pnt->code >= 0)
                    {
                        dx = pt->x - lst->pnt->x;
                        dy = pt->y - lst->pnt->y;
                        if (dx>0 || dx==0.0 && dy>0)
                        {
                            ds = sqrt( dx*dx + dy*dy );
                            if (ds>*maxs) *maxs = ds;
                            if (ds<*mins || *mins==0.0) *mins = ds;
                            *means += ds;
                            nside++;
                        }
                    }
                }
                while (lst != pt->list && lst != NULL);
        }
    if (nside>1) *means /= nside;
    return TRG_ERR_OK;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*  TEST OUTPUT ROUTINES:                                          */
/*                                                                 */
/*  trg_dump    - Dumps triangulation edges to a file              */
/*  trg_dumptrg - Dumps triangulation triangles to a file          */
/*  trg_plot    - Plots triangulation                              */
/*                                                                 */
/*-----------------------------------------------------------------*/

int trg_dump( trg_handle handle, FILE *out )
{
    point pt;
    edgelist lst;

    cur_trg = (trg_rec *) handle;
    fprintf( out, "\n\nDirichelet tesselation dump output\n\n");
    if (cur_trg==NULL)
    {
        fprintf(out,"There is no current triangulation\n");
        return TRG_ERR_OK;
    }
    fprintf( out, "X bounds: %10.4lf %10.4lf\n", (double) cxref, (double) cxmax );
    fprintf( out, "Y bounds: %10.4lf %10.4lf\n\n", (double) cyref, (double) cymax );

    for( pt=lastpoint; pt != NULL; pt = pt->nextpt )
    {
        fprintf( out, "Code %4d:  Position %10.4lf  %10.4lf\n", pt->code,
                 (double) (pt->x+cxref), (double) (pt->y+cyref) );
        lst = pt->list;
        if (lst != NULL) do
            {
                fprintf( out, "        Adjacent to point code %4d",
                         lst->pnt ? lst->pnt->code : 999 );
                fprintf( out, ": End at %10.4lf %10.4lf\n", (double) (lst->endx+cxref),
                         (double) (lst->endy+cyref) );
                lst = lst->next;
            }
            while (lst != pt->list && lst != NULL);
        fprintf( out, "\n");
    }
    return TRG_ERR_OK;
}


int trg_dump_trg( trg_handle handle, FILE *out )
{
    point pt;
    edgelist e1, e2;
    int ntrg;

    cur_trg = (trg_rec *) handle;
    fprintf( out, "\n\nDirichelet triangulation dump output\n\n");
    if (cur_trg==NULL)
    {
        fprintf(out,"There is no current triangulation\n");
        return TRG_ERR_OK;
    }
    fprintf( out, "X bounds: %10.4lf %10.4lf\n", (double) cxref, (double) cxmax );
    fprintf( out, "Y bounds: %10.4lf %10.4lf\n\n", (double) cyref, (double) cymax );

    ntrg = 0;
    for( pt=lastpoint; pt != NULL; pt = pt->nextpt )
    {
        if( pt->code < 0 ) continue;
        e1 = pt->list;
        while( 1 )
        {
            if( e1->pnt->code > pt->code )
            {
                e2 = e1->next;
                if( e2->pnt->code > pt->code )
                {
                    ntrg++;
                    fprintf(out,"%4d %4d %4d %4d\n",ntrg,pt->code,
                            e1->pnt->code,e2->pnt->code);
                }
            }
            e1 = e1->next;
            if( e1 == pt->list ) break;
        }
    }
    return TRG_ERR_OK;
}

int trg_process_trg( trg_handle handle, void *data, void (*process)( int c1, int c2, int c3, void *data))
{
    point pt;
    edgelist e1, e2;
    int ntrg;

    cur_trg = (trg_rec *) handle;
    if (cur_trg==NULL)
    {
        return TRG_ERR_OK;
    }

    for( pt=lastpoint; pt != NULL; pt = pt->nextpt )
    {
        if( pt->code < 0 ) continue;
        e1 = pt->list;
        while( 1 )
        {
            if( e1->pnt->code > pt->code )
            {
                e2 = e1->next;
                if( e2->pnt->code > pt->code )
                {
                    (*process)( pt->code, e1->pnt->code, e2->pnt->code, data );
                }
            }
            e1 = e1->next;
            if( e1 == pt->list ) break;
        }
    }
    return TRG_ERR_OK;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Very inefficient plot routine                               */
/*     Parameters are mode, 0 = plot tesselations, 1 =             */
/*     plot triangles, and drawline a function with parameters     */
/*     x1, y1, x2, y2 which draws a line from x1,y1 to x2,y2       */
/*                                                                 */
/*-----------------------------------------------------------------*/

int trg_plot( trg_handle handle, int mode, void *data, void (*drawline)( void *data, double x1, double y1, double x2, double y2) )
{
    point pt;
    edgelist lst;
    double x1,y1,x2,y2;

    cur_trg = (trg_rec *) handle;
    if(ABORTED) return TRG_ERR_ABORTED;
    for( pt=lastpoint; pt != NULL; pt = pt->nextpt ) if (pt->code >= 0)
        {
            lst = pt->list;
            if (lst != NULL) do
                {
                    if (mode)
                    {
                        if (lst->pnt->code >= 0)
                        {
                            x1 = pt->x+cxref; y1 = pt->y+cyref;
                            x2 = lst->pnt->x+cxref; y2 = lst->pnt->y+cyref;
                        }
                        else { x1 = 1.0; x2 = 0.0;}
                    }
                    else
                    {
                        x1 = lst->endx+cxref; y1 = lst->endy+cyref;
                        x2 = lst->next->endx+cxref; y2 = lst->next->endy+cyref;
                        if (lst->next->pnt->code < 0 ) { x1 = 1.0; x2 = 0.0; }
                    }
                    if( x1<x2 || (x1==x2 && y1<y2 )) (*drawline)(data,x1,y1,x2,y2);
                    lst = lst->next;
                }
                while (lst != pt->list && lst != NULL);
        }
    return TRG_ERR_OK;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Routines to manage an edge list - inserting and deleting      */
/*   elements from the list.  The insertion routine returns a      */
/*   pointer to the new edge, or NULL if the memory allocation     */
/*   failed.                                                       */
/*                                                                 */
/*-----------------------------------------------------------------*/

static void delete_next_edge( edgelist list )
{
    edgelist nextc;
    nextc = list->next;
    list->next = nextc->next;
    free( nextc );
}

static edgelist insert_edge( edgelist list, point pt, double ex, double ey )
{
    edgelist newedge;
    newedge = (edgelist) malloc ( sizeof( struct edgerec ) );
    if (newedge != NULL)
    {
        newedge->endx = ex;
        newedge->endy = ey;
        newedge->pnt  = pt;
        if (list != NULL)
        {
            newedge->next = list->next;
            list->next = newedge;
        }
        else newedge->next = newedge;
    }
    return newedge;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*  Create a new point - updating no_points and the lastpoint      */
/*  pointer.  Returns a pointer to the new point, or NULL if       */
/*  the memory allocation failed.                                  */
/*                                                                 */
/*-----------------------------------------------------------------*/

static point create_point( double x, double y, int code )
{
    point newpt;
    newpt = (point) malloc( sizeof( struct pointrec ));
    if (newpt!=NULL)
    {
        newpt->x = x;
        newpt->y = y;
        newpt->code = code;
        newpt->nextpt = lastpoint;
        newpt->list = NULL;
        lastpoint = newpt;
        if(code>=0) no_points++;
    }
    return newpt;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Procedure to find the nearest existing point to a new       */
/*     point, performing a walk through the tesselation. This      */
/*     implementation checks all neighbours of a point before      */
/*     deciding which one to go to in its search for the closest   */
/*     point.  An alternative that may be more efficient (?)       */
/*     would be to go to the first point that is closer.  Note     */
/*     that the routine avoids negative codes in order to avoid    */
/*     stepping to constraints.                                    */
/*                                                                 */
/*-----------------------------------------------------------------*/


static point find_nearest( double x, double y, point start )
{
    point seek;
    point next;
    edgelist check;
    double distance, checkdst, dx, dy;

    next = start;
    do
    {
        start = next;
        dx = x-start->x; dy = y-start->y;
        distance=dx*dx+dy*dy;
        check = start->list;
        do
        {
            seek = check->pnt;
            if (seek && seek->code>=0)
            {
                dx = x-seek->x; dy=y-seek->y;
                checkdst = dx*dx+dy*dy;
                if (checkdst<distance)
                {
                    distance=checkdst;
                    next=seek;
                }
            }
            check = check->next;
        }
        while (check != start->list);
    }
    while (next != start);
    return start;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Procedure to find the edge of a tesselation tile at which   */
/*     the bisector leaves the tile.  The routine returns a        */
/*     pointer to the edge, and returns the position at which      */
/*     the bisector exits. This routine tests for failure to       */
/*     find the edge.  If this occurs it assumes that the edge     */
/*     has been entered at the corner, and so returns the same     */
/*     point on the next edge.                                     */
/*                                                                 */
/*-----------------------------------------------------------------*/

static edgelist find_end_edge( double x, double y,
                               point start_pt, edgelist start_edge, double *ex, double *ey )
{
    edgelist edge, last_edge;
    double bx, by, bc, cc, lc;

    edge=start_edge;
    last_edge = NULL;

    if (start_pt->code >= 0)
    {
        bx = start_pt->x - x;
        by = start_pt->y - y;
        bc = (bx*bx+by*by)/2;
        cc = (edge->endx - x)*bx + (edge->endy - y)*by - bc;
    }

    for(;;)
    {
        last_edge=edge;
        edge=edge->next;
        if (start_pt->code>=0) lc=cc;
        else
        {

            /* Reached corner of a boundary */

            if (edge->pnt->code<0) { cc=1; lc=0; break; }
            bx = edge->pnt->x - x; by = edge->pnt->y - y; bc = (bx*bx+by*by)/2;
            lc = (last_edge->endx-x)*bx + (last_edge->endy-y)*by - bc;
        }
        cc = (edge->endx - x)*bx + (edge->endy - y)*by - bc;

        /* Edge has been intersected */

        if (lc<0 && cc>=0) break;

        /* Failed to find an edge - assume started at the nearest corner
           of a polygon */

        if ( edge==start_edge )
        {
            last_edge=edge;
            edge=edge->next;
            lc = 0;
            cc = 1;
            break;
        }
    }
    *ex = (lc*edge->endx - cc*last_edge->endx)/(lc-cc);
    *ey = (lc*edge->endy - cc*last_edge->endy)/(lc-cc);
    return edge;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Find the component of a list of edges which is adjacent to    */
/*   a point.                                                      */
/*                                                                 */
/*-----------------------------------------------------------------*/

static edgelist find_adjacent_edge( edgelist list, point pt )
{
    edgelist search;
    search=list;
    while (search->pnt != pt)
    {
        search=search->next;
        if (search==list) {search=NULL; break;}
    }
    return search;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Remove edges of a tile resulting from insertion of the      */
/*     new point, and create a new edge corresponding to the       */
/*     point. The list pointer for the central point is reset in   */
/*     case the corresponding element is deleted from the list.    */
/*                                                                 */
/*-----------------------------------------------------------------*/

static int remove_unwanted_edges( point centpt, edgelist start_edge,
                                  double sx, double sy, edgelist end_edge, double ex, double ey, point pt )
{

    while (start_edge->next != end_edge )
        delete_next_edge( start_edge );
    centpt->list = start_edge;
    start_edge->endx = sx;
    start_edge->endy = sy;
    return insert_edge( start_edge, pt, ex, ey) == NULL ?
           TRG_ERR_ALLOC : TRG_ERR_OK;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Routine to insert the first point, forming the "point"      */
/*     records corresponding to the boundaries within which        */
/*     triangulation is being performed, as well as that for       */
/*     the first point itself.                                     */
/*                                                                 */
/*     The constraints determine the window within which the       */
/*     triangulation will be performed.  It is assumed that this   */
/*     is defined by a rectangle with maximum and minimum x and    */
/*     y values.  We define the initial tesselation on this        */
/*     region when the first point is inserted.                    */
/*                                                                 */
/*     The tesselation being created looks like ...                */
/*                                                                 */
/*       \     c1      /                                           */
/*         \ --------/                                             */
/*      c4  |  pt    |  c2                                         */
/*         / --------\                                             */
/*       /      c3     \                                           */
/*                                                                 */
/*                                                                 */
/*-----------------------------------------------------------------*/

static void define_constraints( double xmin, double xmax, double ymin, double ymax )
{

    cxmin = xmin;
    cxmax = xmax;
    cymin = ymin;
    cymax = ymax;
}


static int insert_first_pt( double x, double y, int code )
{
    point c1, c2, c3, c4, pt;

    c1 = create_point( 0.0, 0.0, -1);
    c2 = create_point( 0.0, 0.0, -2);
    c3 = create_point( 0.0, 0.0, -3);
    c4 = create_point( 0.0, 0.0, -4);
    pt = create_point( x, y, code );

    if ( !(c1 && c2 && c3 && c4 && pt) ) return TRG_ERR_ALLOC;

    c1->list = insert_edge( NULL, c4, cxmin, cymax );
    c1->list = insert_edge( c1->list, pt, cxmax, cymax );
    c1->list = insert_edge( c1->list, c2, cxmax, cymax );

    c2->list = insert_edge( NULL, c1, cxmax, cymax );
    c2->list = insert_edge( c2->list, pt, cxmax, cymin );
    c2->list = insert_edge( c2->list, c3, cxmax, cymin );

    c3->list = insert_edge( NULL, c2, cxmax, cymin );
    c3->list = insert_edge( c3->list, pt, cxmin, cymin);
    c3->list = insert_edge( c3->list, c4, cxmin, cymin );

    c4->list = insert_edge( NULL, c3, cxmin, cymin );
    c4->list = insert_edge( c4->list, pt, cxmin, cymax );
    c4->list = insert_edge( c4->list, c1, cxmin, cymax );

    pt->list = insert_edge( NULL, c1, cxmin, cymax );
    pt->list = insert_edge( pt->list, c4, cxmin, cymin );
    pt->list = insert_edge( pt->list, c3, cxmax, cymin );
    pt->list = insert_edge( pt->list, c2, cxmax, cymax );

    boundarypt = c1;
    return pt->list==NULL ? TRG_ERR_ALLOC : TRG_ERR_OK;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*   Basic routine to insert a new point into the triangulation    */
/*                                                                 */
/*-----------------------------------------------------------------*/


int insert_one_point( double x, double y, int code )
{
    double ex, ey, sx, sy;
    point start_pt, pt, newpt;
    edgelist edge, next_edge, newlist;
    int err;

    /* First check that the point is in range, and that the code is valid */

    if( x<=cxmin || x>=cxmax || y<=cymin || y>=cymax ) return TRG_ERR_BOUNDS;

    if (code<0) code=0;

    /* If it is the first point then use the insert_first_pt routine
       to enter the point and set up the boundary points */

    if( lastpoint==NULL ) return insert_first_pt( x, y, code );


    /* Find nearest point and check that it is not coincident with the
       new point */

    start_pt = find_nearest( x, y, lastpoint );
    sx = start_pt->x-x; sy=start_pt->y-y;
    if (sx*sx+sy*sy<=0) return TRG_ERR_DUPLICATE;

    /* Allocate memory for the new point, and the start of its
       edge list */

    newpt = create_point( x, y, code );
    if (newpt==NULL) return TRG_ERR_ALLOC;

    /* Find the first edge cut by the bisector of the line between the new
       point and the nearest point */

    next_edge = find_end_edge( x, y, start_pt, start_pt->list, &ex, &ey );
    pt = start_pt;

    /* Now chase round the adjacent edges until we get back to where we
       started */

    do
    {
        edge = find_adjacent_edge( next_edge->pnt->list, pt );
        if (edge==NULL) return TRG_ERR_BUG;
        sx = ex; sy = ey;
        pt = next_edge->pnt;
        next_edge = find_end_edge( x, y, pt, edge, &ex, &ey );
        err = remove_unwanted_edges( pt, edge, sx, sy, next_edge, ex, ey, newpt);
        if (err != TRG_ERR_OK) return err;
        newlist = insert_edge( newpt->list, pt, sx, sy );
        if (newlist == NULL) return TRG_ERR_ALLOC;
        if (newpt->list == NULL) newpt->list = newlist;
    }
    while( pt != start_pt );
    return TRG_ERR_OK;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*  Release the memory held by a list of points and their          */
/*  associated lists of edges.                                     */
/*                                                                 */
/*-----------------------------------------------------------------*/

static void free_up_memory( )
{
    edgelist first, curr, nextc;
    point currpt;

    while (lastpoint != NULL )
    {
        first = lastpoint->list;
        nextc = first;
        do
        {
            curr = nextc;
            nextc = curr->next;
            free (curr );
        }
        while (nextc != first && nextc != NULL);
        currpt = lastpoint;
        lastpoint = lastpoint->nextpt;
        free( currpt );
    }
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*    Routine to find the triangle enclosing a point.  Uses        */
/*    routine set_search_point to define the start point for       */
/*    searches to be the point nearest the centroid of all points. */
/*                                                                 */
/*-----------------------------------------------------------------*/

static int set_search_point( )
{
    double cx,cy;
    int npt;
    point pt;

    cx=cy=0;
    npt=0;
    for (pt=lastpoint; pt!=NULL; pt=pt->nextpt) if (pt->code>=0)
        {
            cx += pt->x;
            cy += pt->y;
            npt++;
        }
    if(npt>0) { cx /= npt; cy /= npt;}
    searchpt = find_nearest( cx, cy, lastpoint );
    return npt;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     To find the triangle enclosing a point, first find the      */
/*     nearest point to the test point.  Then search round the     */
/*     adjoining points until the angle between two points         */
/*     includes the angle to the test point.  This is done by      */
/*     looking at the dot product of the vectors to adjacent       */
/*     points with that at 90deg to the test point.  Since we      */
/*     are going anti-clockwise around the point, this should      */
/*     change to negative to positive when we pass the required    */
/*     angle.  The only conditions in which this may not occur     */
/*     are if the test point is coincident with its nearest        */
/*     point, or if it is included in an angle greater than or     */
/*     equal to 180 deg. Having found the angle, we then check     */
/*     whether the point lies inside the triangle formed by the    */
/*     vertices.  If not we repeat the search based on the         */
/*     opposite vertex of the adjacent triangle.                   */
/*                                                                 */
/*-----------------------------------------------------------------*/

static int find_triangle( double xt, double yt, point *p1, point *p2, point *p3 )
{
    double x, y, xp, yp, x1, y1, x2, y2, dp1, dp2;
    edgelist e1, e2, ep1, ep2, end;
    int outside;

    if (no_points<3) return TRG_ERR_NOTRIANGLE;
    if (searchpt == NULL ) set_search_point();
    *p1 = find_nearest( xt, yt, searchpt );

    if (xt==(*p1)->x && yt==(*p1)->y )
    {
        for( e1=(*p1)->list; e1->pnt->code<0; e1=e1->next );
        for( e2=e1->next; e2->pnt->code<0; e2=e2->next );
        outside=TRG_ERR_OK;
    }

    else do
        {

            for( e2=(*p1)->list; e2->pnt->code<0; e2=e2->next );
            xp = (*p1)->x; yp = (*p1)->y;
            x = xt-xp;   y = yt-yp;
            e1 = ep1 = ep2 = NULL;

            end = e2;
            x2 = e2->pnt->x - xp; y2 = e2->pnt->y - yp; dp2 = x*y2 - y*x2;
            for (;;)
            {

                /* If end of loop reached without exiting, use provisional
                   endpoints */

                if( e2 == end && e1 != NULL )
                {
                    e1=ep1; e2=ep2;
                    outside=TRG_ERR_OUTSIDE;
                    break;
                }

                e1 = e2; x1 = x2; y1 = y2; dp1 = dp2;
                for (e2 = e1->next; e2->pnt->code<0; e2=e2->next );
                x2 = e2->pnt->x - xp; y2 = e2->pnt->y - yp; dp2 = x*y2 - y*x2;

                /* If the angle between the two points includes the test point,
                   check whether it is inside the triangle, and exit the loop.
                   If the search needs to proceed to an adjacent triangle,
                   then reset *p1 and flag the condition by NULLing e1 */

                if( dp1 <=0 && dp2 > 0 )
                {
                    if (e2 == e1->next )
                    {
                        x2 -= x1; y2 -= y1; x -= x1; y -= y1;
                        if (x2*y-y2*x >= 0) outside = TRG_ERR_OK;
                        else { *p1 = e1->pnt; e1 = e2 = NULL; }
                    }
                    else outside = TRG_ERR_OUTSIDE;
                    break;
                }

                /* Check for angle>180, and set provisional endpoints if found */

                if( x1*y2 - y1*x2 <= 0.0 ) { ep1 = e1; ep2 = e2; }
            }
        }
        while (e1 == NULL);

    if ( e1 && e2 ) { *p2 = e1->pnt; *p3 = e2->pnt; return outside; }
    else return TRG_ERR_NOTRIANGLE;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*     The following routines are used to define a boundary of     */
/*     points around an existing triangulation.  This code is      */
/*     specifically for the DCDB application.                      */
/*                                                                 */
/*     The boundary is identified by travelling around the         */
/*     edgelist of the boundaries.  It is defined as a list of     */
/*     points, which are subsequently deleted as they are added    */
/*     to the triangulation.                                       */
/*                                                                 */
/*-----------------------------------------------------------------*/


static edgelist bdyedge;

static edgelist first_boundary_edge()
{
    for ( bdyedge=boundarypt->list; bdyedge->pnt->code<0; bdyedge=bdyedge->next);
    return bdyedge;
}

static edgelist next_boundary_edge ( )
{
    point bp;
    bdyedge = bdyedge->next;
    while (bdyedge->pnt->code<0)
    {
        bp = bdyedge->pnt;
        bdyedge=find_adjacent_edge( bp->list, boundarypt );
        boundarypt = bp;
        bdyedge=bdyedge->next->next;
    }
    return bdyedge;
}


/*-----------------------------------------------------------------*/
/*                                                                 */
/*     Two routines to insert boundary points, mainly to ensure    */
/*     that they are not separated by more than the boundary       */
/*     width.  This ensures that they will be adjacent in the      */
/*     Dirichelet tesselation.  The first deals with the           */
/*     straight line segment between two points, and the second    */
/*     with curves.  The latter recursively bisects and angle      */
/*     until the included angle is less than 60 degrees (tested    */
/*     using the dot product - the cosine of an angle less than    */
/*     60 degrees must be greater than 0.5). In the latter the     */
/*     two initial directions are supplied as unit vectors         */
/*     (dx1,dy1) and (dx2,dy2).  The routine assumes the           */
/*     included angle between these two is in the range 0-180.     */
/*                                                                 */
/*-----------------------------------------------------------------*/


static int join_points( double x1, double y1, double x2, double y2,
                        double width, int code )
{
    double dst;
    int nint;
    dst = hypot( x1-x2, y1-y2 );
    nint = ceil(dst/width);
    x2 = (x2-x1)/nint; y2 = (y2-y1)/nint;
    for (; --nint;)
    {
        x1 += x2; y1 += y2;
        if (create_point( x1, y1, code )==NULL) return TRG_ERR_ALLOC;
    }
    return TRG_ERR_OK;
}

static int bisect_angle( double x0, double y0, double dx1, double dy1, double dx2, double dy2, double width, int code )
{
    double dx3, dy3, dst;
    int sts;
    if (dx1*dx2+dy1*dy2 < 0.5)
    {
        dx3 = dx1+dx2; dy3 = dy1+dy2;
        dst = hypot(dx3, dy3);
        if(dst>0.0001) { dx3 /= dst; dy3 /= dst;}
        else { dx3 = dy1; dy3 = -dx1; }
        sts = bisect_angle( x0, y0, dx1, dy1, dx3, dy3, width, code );
        if (sts != TRG_ERR_OK) return sts;
        if (create_point( x0+dx3*width, y0+dy3*width, code)==NULL) return TRG_ERR_ALLOC;
        sts = bisect_angle( x0, y0, dx3, dy3, dx2, dy2, width, code );
        if (sts != TRG_ERR_OK) return sts;
    }
    return TRG_ERR_OK;
}

/*-----------------------------------------------------------------*/
/*                                                                 */
/*  The main routine for creating the boundary.                    */
/*                                                                 */
/*-----------------------------------------------------------------*/

#define COS45 0.707107

static int make_boundary( double width, int code )
{
    edgelist lastedge, edge, nextedge, firstedge;
    point thispt, nextpt;
    double dxl,dyl,dxn,dyn,dsn,xv,yv,dprd;
    point savelast, firstpt;
    int   savenpt;
    int   status;
    int   pnt;
    static double cosa[] = { 1.0, COS45, 0.0, -COS45, -1.0, -COS45, 0.0, COS45 };
    static double sina[] = { 0.0, COS45, 1.0, COS45, 0.0, -COS45, -1.0, -COS45 };

    /* If there are no points there is no border to add */

    if( no_points<1 ) return TRG_ERR_NOPOINTS;

    /* With just one point - add a border around that point */

    if( no_points==1)
    {
        xv = lastpoint->x;
        yv = lastpoint->y;
        for ( pnt=0; pnt<sizeof(cosa)/sizeof(double); pnt++)
        {
            status = insert_one_point(xv+width*cosa[pnt],yv+width*sina[pnt],code);
            if (status>=TRG_FATAL_ERROR) return status;
        }
        return TRG_ERR_OK;
    }

    /* Form a list of boundary points - save the lastpoint and no_points
       variables, so that they can be restored after the list is completed */

    if(code<0) code=0;

    savelast = lastpoint;
    savenpt = no_points;
    lastpoint = NULL;
    firstpt = NULL;

    /* Initiallize the loop that will follow around the boundary -
       note that the edge returned by first_boundary_edge is not
       used as this may missed by next_boundary_edge when we have
       travelled around the boundary.  Instead we start at the
       first edge returned by next_boundary_edge. */

    lastedge = NULL;
    first_boundary_edge();
    firstedge = edge = next_boundary_edge();
    nextedge = next_boundary_edge();
    thispt = edge->pnt;
    nextpt = nextedge->pnt;
    dxn = nextpt->x - thispt->x; dyn = nextpt->y - thispt->y;
    dsn = hypot ( dxn, dyn ); dxn /= dsn; dyn /= dsn;

    /* Loop to chase around boundary */

    for(;;)
    {
        if ( edge==firstedge && lastedge != NULL ) break;
        lastedge = edge;
        edge = nextedge;
        thispt = nextpt;
        nextedge = next_boundary_edge();
        nextpt = nextedge->pnt;
        dxl = dxn; dyl = dyn;
        dxn = nextpt->x - thispt->x; dyn = nextpt->y - thispt->y;
        dsn = hypot( dxn, dyn ); dxn /= dsn; dyn /= dsn;

        /* Check that this is not a re-entrant corner */

        if (dxl*dyn-dyl*dxn < 0.0 || nextpt == lastedge->pnt )
        {
            xv = thispt->x - width*dyl; yv =thispt->y + width*dxl;

            /* Connect from the previous point (if any) to the first point
               at the vertex */

            if ( lastpoint != NULL )
            {
                status = join_points(lastpoint->x,lastpoint->y,xv,yv,width,code);
                if (status != 0) break;
            }

            /* If the change of angle at the vertex is small, insert a
               single point bisecting the angle.  Otherwise insert a range
               of points */

            dprd = dxl*dxn+dyl*dyn;

            if (dprd>0.75)
            {
                dxl+=dxn; dyl+=dyn;
                dprd = width/sqrt(dxl*dxl+dyl*dyl);
                if (create_point( thispt->x-dyl*dprd,
                                  thispt->y+dxl*dprd,
                                  code) == NULL) { status=2; break;}
                if (firstpt==NULL)  firstpt=lastpoint;
            }

            else
            {

                /* Insert the first point of the vertex */

                if (create_point( xv, yv, code )==NULL) { status=2; break;}
                if(firstpt==NULL) firstpt=lastpoint;

                /* Create the boundary points around the angle of the vertex */

                status = bisect_angle(thispt->x,thispt->y,-dyl,dxl,-dyn,dxn,width,code);
                if (status != 0) break;

                /* Create the last point of the vertex */

                xv = thispt->x - width*dyn; yv = thispt->y + width*dxn;
                if (create_point( xv, yv, code )==NULL) { status=2; break;}
            }
        }
    }

    /* Now join the last point of the boundary back to the first */

    if( status==0 && firstpt != NULL )
        status = join_points( lastpoint->x, lastpoint->y,
                              firstpt->x, firstpt->y, width, code );

    /* If an error has been generated - dispose of the list of boundary points */

    if(status != 0) free_up_memory();

    /* Restore the original lastpoint */

    firstpt = lastpoint;
    lastpoint = savelast;
    no_points = savenpt;

    /* Dispose of the boundary list, adding the points into the triangulation
       as they are deleted */

    if (status==0) while ( firstpt != NULL )
        {
            xv = firstpt->x; yv = firstpt->y;
            thispt = firstpt->nextpt;
            free(firstpt);
            firstpt=thispt;
            if(status<TRG_FATAL_ERROR) status=insert_one_point( xv, yv, code );
        }
    return status<TRG_FATAL_ERROR ? TRG_ERR_OK : status;
}

