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
  double tunit;
};

typedef struct data data;

struct model {
  int *qvector;
  int *evector;
  int *constraint;
  int *miscconstraint;
  int nst;
  int nms;
  int nintens;
  int nmisc;
  int ncoveffs;
  int nmisccoveffs;
  int covmatch;
  int death;
  int exacttimes;
  double *intens;
  double *coveffect;
  double *miscprobs;
  double *misccoveffect;
  double *initprobs;
};

typedef struct model model;

void msmCEntry(int *do_what,      /* 1 = eval likelihood, 2 = Viterbi, 3 = prediction */
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
	       double *initprobs, /* initial state occupancy probabilities */
	       int *nstates,      /* number of Markov states */
	       int *nms,          /* number of underlying states which can be misclassified */
	       int *nintens,      /* number of intensity parameters */
	       int *nmisc,        /* number of misclassification rates */
	       int *nobs,         /* number of observations in data set */
	       int *npts,         /* number of individuals in data set */
	       int *ncovs,        /* number of covariates on transition rates */
	       int *ncoveffs,     /* number of distinct covariate effect parameters */
	       int *nmisccovs,    /* number of covariates on misclassification probabilities */
	       int *nmisccoveffs, /* number of distinct misclassification covariate effect parameters */
	       int *covmatch,     /* use the covariate value from the previous or next observation */
	       int *death,        /* indicator for death state */
	       double *tunit,     /* time unit in days */
	       int *exacttimes,   /* indicator for exact transition times */
	       int *fixedpars,    /* which parameters to fix */
	       double *predtimes, /* prediction times for one-step-ahead prediction */
	       int *npreds,       /* number of prediction times */
	       double *returned   /* returned -2 log likelihood */
  );

void msmLikelihood (data *d, model *m, int misc, double *returned);

/* Fill a parameter vector with either the current values from the optimisation or the fixed inital values */

void fillparvec(double *parvec, /* named vector to fill (e.g. intens = baseline intensities) */
		double *params, /* current values of parameters being optimised over */
		double *allinits, /* full vector of initial values */
		int *fixi,  /* indices of allinits which are fixed */
		int ni,    /* length of parvec */
		int *ifix,  /* current index into fixi */
		int *iopt,  /* current index into params */
		int *iall   /* current index into allinits */
  );

double likmisc(int pt, data *d, model *m);
void UpdateLik(int state, double dt, int k, int last, int predict_death, data *d, model *m, 
	       double *cumprod, double *newprod, double lweight_old, double *lweight_new);
void AddCovs(int obs, data *d, model *m, double *newintens);
void AddMiscCovs(int obs, data *d, model *m, double *newp);
double PObsTrue(int obst,      /* observed state */
		int tst,       /* true state */
		double *miscprobs, /* misclassification probabilities */
		model *m
  );
double liksimple(data *d, model *m);
double liksimple_fromto(data *d, model *m);
void Viterbi(data *d, model *m, double *fitted);


double pijt(int i, int j, double t, vector intens, int *qvector, int nstates, int exacttimes);
void Pmat(Matrix pmat, double t, vector intens, int *qvector, int nstates, int exacttimes);
void FillQmatrix(int *qvector, vector intens, Matrix qmat, int nstates);
void MatrixExp(Matrix mat, int n, Matrix expmat, double t);
int repeated_entries(vector vec, int n);
void MatrixExpSeries(Matrix mat, int n, Matrix expmat, double t);
void MatTranspose(Matrix A, Matrix AT, int n);
void MatInv(Matrix A, Matrix Ainv, int n);
void MultMat(Matrix A, Matrix B, int arows, int acols, int bcols, Matrix AB);
void MultMatDiag(Matrix A, vector diag, int n, Matrix AB);
void FormIdentity(Matrix A, int n);
