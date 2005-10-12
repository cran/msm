#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Applic.h>

 /* index to treat a vector as a matrix. ith row, jth column. Fills columns first, as in R */
#define MI(i, j, nrows) ( (int) ((j)*(nrows) + (i)) )
/* index to treat a vector as a 3-dimensional array. Left-most index varies fastest, as in R */
#define MI3(i, j, k, n1, n2) ( (int) ((k)*(n1*n2) + (j)*(n1) + (i)) )

/* Macros to switch quickly between C and S memory handling. Currently not used */

#define USE_CALLOC
/* #define USE_SALLOC */

#ifdef USE_CALLOC
#define MSM_ALLOC(length, type) Calloc((length), type)
#define MSM_FREE(var) Free((var))
#else
#define MSM_ALLOC(length, type) (type *) S_alloc((length), sizeof(type))
#define MSM_FREE(var) 
#endif 

typedef double * Array3;
typedef double * Matrix;
typedef int * iMatrix;
typedef double * vector;
typedef int * ivector;

struct msmdata {
    /* for non-hidden model */
    int *fromstate;
    int *tostate;
    double *timelag;
    double *cov;
    int *whichcov;
    int *nocc;
    int *whicha;
    int *obstype;

    /* for hidden model */
    int *subject;
    double *time;
    double *obs; /* observed state or any other HMM observation */
    int *firstobs;
    int *whichcovh;

    int nobs;
    int npts;
};

struct qmodel {
    int nst;
    int npars;
    int *ivector;
    double *intens;
};

struct qcmodel {
    int *ncovs;
    double *coveffect;
};

struct cmodel {
    int ncens;
    int *censor;
    int *censstind;
    int *censstates;
};

struct hmodel {
    int hidden;
    int *models;
    int *npars;
    int *firstpar;
    int *ncovs;
    double *pars;
    int totpars; 
    double *coveffect;
    int *links;
    double *initp;
};

typedef struct msmdata msmdata;
typedef struct qmodel qmodel;
typedef struct qcmodel qcmodel;
typedef struct cmodel cmodel;
typedef struct hmodel hmodel;

double qij(int i, int j, vector intens, ivector qvector, int nstates);
double pijdeath(int r, int s, Matrix pmat, vector intens, ivector qvector, int n);
void Pmat(Matrix pmat, double t, vector intens, int *qvector, int nstates, int exacttimes, int debug);
void DPmat(Array3 dpmat, double t, vector intens, ivector qvector, int n, int np, int exacttimes, int debug);
void dpijdeath(int r, int s, Array3 dpmat, Matrix pmat, vector intens, ivector qvector, int n, int np, vector dcontrib);
int repeated_entries(vector vec, int n);

double logit(double x);
double expit(double x);
double identity(double x);
int all_equal(double x, double y);

void MatrixExpPadeR(double *ExpAt, double *A, int *n, double *t);
