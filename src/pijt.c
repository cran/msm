/* *****************************************************************
   PROGRAM: pijt.c 
   AUTHOR:  Chris Jackson
   DATE:    November 2001

   Linear algebra routines for obtaining Markov transition probabilities 
   (the P-matrix) from transition intensities (the Q-matrix)

   ******************************************************************  */ 

#include "lik.h"

/* Calculates i-j transition probability in time t given an intensity matrix */

double pijt(int i, int j, double t, vector intens, ivector qvector, int nstates, int exacttimes)
{
    Matrix qmat = (double *) S_alloc( (nstates)*(nstates), sizeof(double));
    Matrix exptq = (double *) S_alloc( (nstates)*(nstates), sizeof(double));
    double pr, pii;
    FillQmatrix(qvector, intens, qmat, nstates);
    MatrixExp(qmat, nstates, exptq, t);
    if (exacttimes) {
	pii = exp(t * qmat[MI(i, i, nstates)] );
	pr = ( i==j  ?  pii  : pii * qmat[MI(i, j, nstates)] );
    }
    else 
	pr = exptq[MI(i, j, nstates)];
    return ((pr < 0) ? 0 : pr); /* avoid problems with returning tiny negative numbers */
}

/* Returns i-j transition intensity time t given an intensity matrix */

double qij(int i, int j, vector intens, ivector qvector, int nstates)
{
    Matrix qmat = (double *) S_alloc( (nstates)*(nstates), sizeof(double));
    FillQmatrix(qvector, intens, qmat, nstates);
    return qmat[MI(i,j,nstates)];
}

/* Calculates the whole transition matrix in time t given an intensity matrix */

void Pmat(Matrix pmat, double t, vector intens, ivector qvector, int nstates, int exacttimes)
{
    int i,j;
    double pii;
    Matrix qmat = (double *) S_alloc( (nstates)*(nstates), sizeof(double));
    FillQmatrix(qvector, intens, qmat, nstates);
    if (exacttimes) { 
	for (i=0; i<nstates; ++i) 
	    for (j=0; j<nstates; ++j) {
		pii = exp(t * qmat[MI(i, i, nstates)] );
		pmat[MI(i, j, nstates)] = ( i==j  ?  pii  : pii * qmat[MI(i, j, nstates)] );
	    }
    }
    else {
	MatrixExp(qmat, nstates, pmat, t);
    }
}

/* Fills the required entries of the intensity-matrix with the current intensities */
/* Also works in the same way for the misclassification matrix */

void FillQmatrix(ivector qvector, vector intens, Matrix qmat, int nstates)
{
    int i, j, k=0;
    for (i=0; i<nstates; ++i) {
	qmat[MI(i, i, nstates)] = 0;
	for (j=0; j<nstates; ++j) {
	    if (j != i) {
		qmat[MI(i, j, nstates)] = 0;
		if (qvector[MI(i, j, nstates)] == 1) {
		    qmat[MI(i, j, nstates)] = intens[k];
		    qmat[MI(i, i, nstates)] -= intens[k];
		    ++k;
		}
	    }
	}
    }
}

/* Calls EISPACK/LINPACK routines to find exponential of a matrix */
/* If matrix has repeated eigenvalues, then use power series approximation instead of eigensystem decomposition */

void MatrixExp(Matrix mat, int n, Matrix expmat, double t)
{  
    Matrix work = (double *) S_alloc(n*n, sizeof(double));
    iMatrix worki = (int *) S_alloc(n*n, sizeof(int));
    vector revals = (double *) S_alloc(n, sizeof(double));
    vector ievals = (double *) S_alloc(n, sizeof(double));
    Matrix evecst = (double *) S_alloc(n*n, sizeof(double));
    Matrix evecs = (double *) S_alloc(n*n, sizeof(double));
    Matrix evecsinv = (double *) S_alloc(n*n, sizeof(double));
    Matrix matt = (double *) S_alloc(n*n, sizeof(double));
    int i, err, matz = 1;

    MatTranspose(mat, matt, n);
    /* calculate eigensystem */
    F77_CALL(rg) (&n, /* leading dim */
		  &n, /* no of rows */
		  matt, 
		  revals, 
		  ievals,
		  &matz,
		  evecst,
		  worki,  /* workspace */
		  work,  /* workspace */
		  &err);
    if (repeated_entries (revals, n)){
	MatrixExpSeries(mat, n, expmat, t);
    }
    else {
	MatTranspose(evecst, evecs, n);
	for (i=0; i<n; ++i)
	    revals[i] = exp(revals[i] * t);
	MatInv(evecs, evecsinv, n);
	MultMatDiag(revals, evecsinv, n, work);
	MultMat(evecs, work, n, n, n, expmat);
    }
}

/* Tests if a vector has any non-unique entries */

int repeated_entries(vector vec, int n)
{
    int i, j;
    for (i=1; i<n; ++i)
	for (j=0; j<i; ++j)
	    if (vec[j] == vec[i]) 
		return 1;
    return 0;
}

/* Tests if an integer a is a member of the integer set b */

int is_element(int a, int *b, int n)
{
    int i;
    for (i=0; i<n; ++i)
	if (b[i] == a)
	    return 1;
    return 0;
}

/* Calculate a matrix exponential using a power series approximation */
/* Adapted from mexp in Jim Lindsey's rmutil library */

void MatrixExpSeries(Matrix A, int n, Matrix expmat, double t)
{
    int i, j;
    int order = 20;   /* number of terms in series */
    int underflow_correct = 3;
    Matrix Apower = (double *) S_alloc(n*n, sizeof(double));
    Matrix Temp = (double *) S_alloc(n*n, sizeof(double));
    Matrix AA = (double *) S_alloc(n*n, sizeof(double));
    for (i=0; i<(n*n); ++i)
	AA[i] = A[i] * (t / pow(2, underflow_correct));
    FormIdentity(expmat, n);
    FormIdentity(Apower, n);
    for (i=1; i<=order; i++) {
	MultMat(AA, Apower, n, n, n, Temp);
	for (j=0; j<(n*n); ++j){
	    Apower[j] = Temp[j] / i;
	    expmat[j] += Apower[j];
	}
    }
    for (i=0; i<underflow_correct; ++i){
	MultMat(expmat, expmat, n, n, n, Temp);
	for (j=0; j<(n*n); ++j)
	    expmat[j] = Temp[j];
    }
}

/* Transpose a matrix */

void MatTranspose(Matrix A, Matrix AT, int n)
{
    int i, j;
    for (i=0; i<n; ++i)
	for (j=0; j<n; ++j)
	    AT[MI(i,j,n)] = A[MI(j,i,n)];
}

/* Invert a matrix by calling LINPACK QR decomposition routines */

void MatInv(Matrix A, Matrix Ainv, int n)
{
    int i, j, rank;
    double tol=1e-07;
    Matrix work = (double *) S_alloc(n*n, sizeof(double));
    Matrix qraux = (double *) S_alloc(n*n, sizeof(double));
    int info, *pivot = (int *) S_alloc(n, sizeof(int));
    Matrix ident = (double *) S_alloc(n*n, sizeof(double));
    Matrix temp = (double *) S_alloc(n*n, sizeof(double));
    for (i=0; i<(n*n); ++i)
	temp[i] = A[i];
    F77_CALL(dqrdc2) (temp, &n, &n, &n, &tol, &rank, qraux, pivot, work);
    for (i=0; i<n; ++i)
	for (j=0; j<n; ++j){
	    if (i==j)
		ident[MI(i,j,n)] = 1;
	    else
		ident[MI(i,j,n)] = 0;
	}
    F77_CALL(dqrcf) (temp, &n, &rank, qraux, ident, &n, Ainv, &info);
}

/* Multiplies two matrices together */

void MultMat(Matrix A, Matrix B, int arows, int acols, int bcols, Matrix AB)
{
    int i, j, k;  
    for (i = 0; i < arows; i++) {
	for (j = 0; j < bcols; j++) {
	    AB[MI(i, j, bcols)] = 0;
	    for (k = 0; k < acols; k++) 
		AB[MI(i, j, bcols)] += A[MI(i, k, acols)] * B[MI(k, j, bcols)];
	}
    }
}

/* Pre-multiplies a general matrix by a diagonal matrix (given by a vector) */

void MultMatDiag(vector diag, Matrix B, int n, Matrix AB)
{
    int i, j;
    for (i = 0; i < (n*n); ++i)
	AB[i] = 0;
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    AB[MI(i, j, n)] += diag[i] * B[MI(i, j, n)];
	}
    }
}

/* Set A to be an n x n identiy matrix */

void FormIdentity(Matrix A, int n)
{
    int i;
    for (i = 0; i < (n*n); ++i)
	A[i] = 0;
    for (i = 0; i < n; ++i)
	A[MI(i, i, n)] = 1;
}
