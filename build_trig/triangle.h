/* TRIANGLE.H - header file for triangle.c */

/* Definition of handle for triangulation, and arithmetic precision */

typedef void *trg_handle;

/* Procedures defined in the triangulation routines */

int trg_init( trg_handle *handle, double wxmin, double wxmax, double wymin, double wymax);
int trg_add_point( trg_handle handle, double x, double y, int code );
int trg_add_border( trg_handle handle, double width, int code );
int trg_triangle( trg_handle handle, double x, double y,
                  double *x1, double *y1, int *c1,
                  double *x2, double *y2, int *c2,
                  double *x3, double *y3, int *c3);
int trg_lin_interp( trg_handle handle, double x, double y,
                    int *c1, double *m1,
                    int *c2, double *m2,
                    int *c3, double *m3 );
int trg_sizes( trg_handle handle, double *mins, double *maxs, double *means);
int trg_term( trg_handle *handle );
int trg_dump( trg_handle handle, FILE *out );
int trg_dump_trg( trg_handle handle, FILE *out );
int trg_process_trg( trg_handle handle, void *data, void (*process)( int c1, int c2, int c3, void *data));
int trg_plot( trg_handle handle, int mode, void *data, void (*drawline)( void *data, double x1, double y1, double x2, double y2) );

/* Status codes returned by the routines */

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


#define TRG_PLOT_TESSELATION 0
#define TRG_PLOT_TRIANGLES   1
