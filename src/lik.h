#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Applic.h>
/* #include <stdio.h> */

#define MI(i, j, ncols) ( (int) ((i)*(ncols) + (j)) ) /* index to treat a vector as a matrix. Fills rows first  */
#define logit(x) (log ( (x) / (1 - (x))))
#define expit(x) (exp (x) / (1 + exp(x)))

typedef double * Matrix;
typedef int * iMatrix;
typedef double * vector;
typedef int * ivector;

struct data {
    int *subject;
    double *time;
    int *state;
    int *tostate;
    int fromto;
    double *cov;
    double *misccov;
    int nobs;
    int npts;
    int ncovs;
    int nmisccovs;
};

typedef struct data data;

struct model {
    int *qvector;
    int *evector;
    int *constraint;
    int *miscconstraint;
    int *baseconstraint;
    int *basemiscconstraint;
    int nst;
    int nms;
    int nintens;
    int nintenseffs;
    int nmisc;
    int nmisceffs;
    int ncoveffs;
    int nmisccoveffs;
    int covmatch;
    int ndeath;
    int *death;
    int ncens; int *censor; int *censstates; int *censstind;
    int exacttimes;
    double *intens;
    double *coveffect;
    double *miscprobs;
    double *misccoveffect;
    double *initprobs;
};

typedef struct model model;

void msmCEntry(int *do_what,      /* 1 = eval likelihood, 2 = Viterbi */
	       double *params,    /* full parameter vector */
	       double *allinits,  /* all initial values */
	       int *misc,
	       int *p,            /* number of parameters */
	       int *subjvec,      /* vector of subject IDs */
	       double *timevec,   /* vector of observation times (fromto==F), or time lags (fromto==T) */
	       int *statevec,     /* vector of observed states (fromto==F), or from-states (fromto==T) */
	       int *tostatevec,     /* vector of observed to-states (fromto==T) */
	       int *fromto,       /* if data is of the alternative form "from state, to state, time interval" */
	       int *qvector,      /* vectorised matrix of allowed transition indicators */
	       int *evector,      /* vectorised matrix of allowed misclassification indicators */
	       double *covvec,    /* vectorised matrix of covariate values */
	       int *constrvec,    /* list of constraints for each covariate */
	       double *misccovvec,/* vectorised matrix of misclassification covariate values */
	       int *miscconstrvec,/* list of constraints for each misclassification covariate */
	       int *baseconstraint, /* constraints on baseline transition intensities */
	       int *basemiscconstraint, /* constraints on baseline misclassification probabilities */
	       double *initprobs, /* initial state occupancy probabilities */
	       int *nstates,      /* number of Markov states */
	       int *nms,          /* number of underlying states which can be misclassified */
	       int *nintens,      /* number of intensity parameters */
	       int *nintenseffs,  /* number of distinct intensity parameters */
	       int *nmisc,        /* number of misclassification rates */
	       int *nmisceffs,    /* number of distinct misclassification rates */
	       int *nobs,         /* number of observations in data set */
	       int *npts,         /* number of individuals in data set */
	       int *ncovs,        /* number of covariates on transition rates */
	       int *ncoveffs,     /* number of distinct covariate effect parameters */
	       int *nmisccovs,    /* number of covariates on misclassification probabilities */
	       int *nmisccoveffs, /* number of distinct misclassification covariate effect parameters */
	       int *covmatch,     /* use the covariate value from the previous or next observation */
	       int *ndeath,        /* number of death states */
	       int *death,        /* vector of indices of death states */
	       int *ncens,
	       int *censor,
	       int *censstates,
	       int *censstind,
	       int *exacttimes,   /* indicator for exact transition times */
	       int *nfix,    /* number of fixed parameters */
	       int *fixedpars,    /* which parameters to fix */
	       double *returned   /* returned -2 log likelihood */
	       );

void msmLikelihood (data *d, model *m, int misc, double *returned);

/* Fill a parameter vector with either the current values from the optimisation or the fixed inital values */

void fillparvec(double *parvec, /* named vector to fill (e.g. intens = baseline intensities) */
		int p, /* number of parameters optimised over */
		double *params, /* current values of parameters being optimised over */
		double *allinits, /* full vector of initial values */
		int nfix,   /* number of fixed parameters, ie length of fixedpars */
		int *fixedpars,  /* indices of allinits which are fixed */
		int ni,    /* length of parvec */
		int *ifix,  /* current index into fixedpars */
		int *iopt,  /* current index into params */
		int *iall   /* current index into allinits */
		);

double likmisc(int pt, data *d, model *m);
void UpdateLik(int *current, int nc, double dt, int k, int last, data *d, model *m, 
	       double *cumprod, double *newprod, double *lweight);
void AddCovs(int obs, data *d, model *m, double *newintens);
void AddMiscCovs(int obs, data *d, model *m, double *newp);
void GetCensored (int obs, model *m, int *nc, int **states);
double PObsTrue(int obst,      /* observed state */
		int tst,       /* true state */
		double *miscprobs, /* misclassification probabilities */
		model *m
		);
double liksimple(data *d, model *m);
double liksimple_fromto(data *d, model *m);
void Viterbi(data *d, model *m, double *fitted);


double pijt(int i, int j, double t, vector intens, int *qvector, int nstates, int exacttimes);
double qij(int i, int j, vector intens, ivector qvector, int nstates);
void Pmat(Matrix pmat, double t, vector intens, int *qvector, int nstates, int exacttimes);
void FillQmatrix(int *qvector, vector intens, Matrix qmat, int nstates);
void MatrixExp(Matrix mat, int n, Matrix expmat, double t);
int repeated_entries(vector vec, int n);
int is_element(int a, int *b, int n);
void MatrixExpSeries(Matrix mat, int n, Matrix expmat, double t);
void MatTranspose(Matrix A, Matrix AT, int n);
void MatInv(Matrix A, Matrix Ainv, int n);
void MultMat(Matrix A, Matrix B, int arows, int acols, int bcols, Matrix AB);
void MultMatDiag(Matrix A, vector diag, int n, Matrix AB);
void FormIdentity(Matrix A, int n);
