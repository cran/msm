/* *****************************************************************

   pijt.c: Linear algebra routines for msm. To obtain Markov
   transition probabilities (the P-matrix) from transition intensities
   (the Q-matrix).

   Copyright (c) Chris Jackson 2001-2005

   Pade approximants method of calculating matrix exponentials is
   based on matexp.cc from JAGS, Copyright (c) 2004, Martyn Plummer.  
   http://www-fis.iarc.fr/~martyn/software/jags
   
   *  This program is free software; you can redistribute it and/or modify
   *  it under the terms of the GNU General Public License as published by
   *  the Free Software Foundation; either version 2 of the License, or
   *  (at your option) any later version.
   
   ******************************************************************  */ 

#include "msm.h"
/* Crude benchmarks have shown that using the eigensystem routines
   from LINPACK, and the matrix inversion routines from LAPACK, is the
   faster combination */
#define _USE_LINPACK_EIGEN_
#define _USE_LAPACK_INVERSE_
#define _MEXP_METHOD_ 1 /* 1 for Pade approximation, 2 for series. Pade is more robust. */
#include "R_ext/Lapack.h"


/* Set A to be an n x n identity matrix */

void FormIdentity(Matrix A, int n)
{
    int i;
    for (i = 0; i < (n*n); ++i)
	A[i] = 0;
    for (i = 0; i < n; ++i)
	A[MI(i, i, n)] = 1;
}

/* Invert a matrix by calling LINPACK QR decomposition routines */
/* DGETRF/DGETRI method (current)  is from LAPACK. This is the fastest.   */
/* DQR method (used in <=0.4.1) is from LINPACK */ 

void MatInv(Matrix A, Matrix Ainv, int n)
{
    int i, j;
    Matrix temp = (Matrix) Calloc(n*n, double);
    Matrix work = (Matrix) Calloc(n*n, double);
    int nsq=n*n, info, *pivot = Calloc(n, int);
    for (i=0; i<nsq; ++i)
      temp[i] = A[i];
    F77_CALL(dgetrf) (&n, &n, temp, &n, pivot, &info);
    if (info < 0)
	REprintf("error code %d from Lapack routine dgetrf\n", info);
    F77_CALL(dgetri) (&n, temp, &n, pivot, work, &nsq, &info);
    if (info < 0)
	REprintf("error code %d from Lapack routine dgetri\n", info);
    for (i=0; i<n; ++i)
	for (j=0; j<n; ++j)
	    Ainv[MI(i,j,n)] = temp[MI(i,j,n)];
    Free(work); Free(pivot); Free(temp); 
}

void MatInvDQR(Matrix A, Matrix Ainv, int n)
{
    int i, rank;
    Matrix temp = (Matrix) Calloc(n*n, double);
    Matrix work = (Matrix) Calloc(n*n, double);
    Matrix qraux = (Matrix) Calloc(n*n, double);
    Matrix ident = (Matrix) Calloc(n*n, double); 
    int nsq=n*n, info, *pivot = Calloc(n, int);
    double tol=1e-07;
    for (i=0; i<nsq; ++i)
      temp[i] = A[i];
    F77_CALL(dqrdc2) (temp, &n, &n, &n, &tol, &rank, qraux, pivot, work); 
    FormIdentity(ident, n); 
    F77_CALL(dqrcf) (temp, &n, &rank, qraux, ident, &n, Ainv, &info); 
    if (info < 0)
	REprintf("error code %d from Linpack routine dqrcf\n", info);
    Free(temp); Free(work); Free(qraux); Free(ident);Free(pivot); 
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

/* Calculate a matrix exponential using a power series approximation */
/* Adapted from mexp in Jim Lindsey's rmutil library */

void MatrixExpSeries(Matrix A, int n, Matrix expmat, double t)
{
    int i, j;
    int order = 20;   /* number of terms in series */
    int underflow_correct = 3;
    Matrix Apower = Calloc(n*n, double);
    Matrix Temp = Calloc(n*n, double);
    Matrix AA = Calloc(n*n, double);
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
    Free(Apower); Free(Temp); Free(AA);
}

/* Pade approximants method of calculating matrix exponentials 
   From matexp.cc from JAGS by Martyn Plummer  Copyright (c) 2004
   http://www-fis.iarc.fr/~martyn/software/jags
*/

static void solve(double *X, double const *A, double const *B, int n)
{
    /* Solve AX = B, where all matrices are square */
    int N = n*n;
    double *Acopy = Calloc(N, double);
    double *work = Calloc(N, double);
    int *ipiv = Calloc(N, int);
    int info = 0;
    static int c_1 = 1;
    F77_CALL(dcopy)(&N, A, &c_1, Acopy, &c_1);
    F77_CALL(dcopy)(&N, B, &c_1, X, &c_1);
    F77_CALL(dgesv)(&n, &n, Acopy, &n, ipiv, X, &n, &info);
    if (info < 0)
	REprintf("argument %d of Lapack routine dgesv had illegal value\n", -info);
    if (info > 0)
	REprintf("Lapack routine dgesv: system is exactly singular\n");
    Free(Acopy); 
    Free(ipiv);
    Free(work);
}

static void
padeseries (double *Sum, double *A, int m, int order, 
            double scale, double *Temp)
{
    int i, j, r; 
    int N = m*m;
    FormIdentity(Sum, m);
    for (j = order; j >= 1; --j) {
	double s = (order-j+1) / (j*(2*order-j+1) * scale);
	MultMat(Sum, A, m, m, m, Temp);
	for (i = 0; i < N; ++i) {
	    Sum[i] = Temp[i] * s;
	}
	for (r = 0; r < m; ++r) {
	    Sum[r*m+r] += 1;
	}
    }
}

void 
MatrixExpPade(double *ExpAt, double *A, int n, double t)
{
  /* Calculate exp(A*t) by diagonal Pade approximation with scaling and
     squaring */
    int i, j; 
  int order = 8;
  int N = n*n;
  double *workspace =  Calloc( 4*N, double);
  double * Temp = workspace;
  double * At = workspace + N;
  double * Num = workspace + 2*N;
  double * Denom = workspace + 3*N;
  double l1 = F77_CALL(dlange)("1", &n, &n, At, &n, 0); /* L-1 norm */
  double linf = F77_CALL(dlange)("i", &n, &n, At, &n, Temp); /* L-Infinity norm */
  double K = (log(l1) + log(linf))/log(4);
  int npower =  (int)(K) + 4;
  double scale = 1;

  /* Multiply by t */

  for (i = 0; i < N; ++i) {
    At[i] = A[i] * t;  
  }

  /* Scale the matrix by a power of 2 */
  /* 
     The expression below is not clear because it is optimized.  The
     idea is that sqrt(l1 * linf) is an upper bound on the L2 norm of
     the matrix At (i.e the largest eigenvalue). We want to take the
     log, to base 2 of this to get the smallest K, st ||At/2^K|| <= 1.
  */
  if (npower < 0) {
      npower = 0;
  }
  for (i = 0; i < npower; ++i) {
    scale *= 2;
  }

  /* Calculate exp(A/scale) by Pade series  */

  padeseries (Num, At, n, order, scale, Temp);
  for (i = 0; i < N; ++i) {
    At[i] = -At[i];
  }
  padeseries (Denom, At, n, order, scale, Temp);
  solve(ExpAt, Denom, Num, n);

  /* Now repeatedly square the result */
  for (i = 0; i < npower; ++i) {
    for (j = 0; j < N; ++j) {
      Temp[j] = ExpAt[j];
    }
    MultMat(Temp, Temp, n, n, n, ExpAt);
  }
  Free(workspace);
}

void MatrixExpPadeR(double *ExpAt, double *A, int *n, double *t)
{
    MatrixExpPade(ExpAt, A, *n, *t);
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

/* Calls EISPACK/LINPACK routines to find exponential of a matrix by eigensystem decomposition */
/* If matrix has repeated eigenvalues, then use Pade approximation instead */

void MatrixExp(Matrix mat, int n, Matrix expmat, double t, int debug)
{  
    int i, err=0, nsq=n*n, matz = 1;
    Matrix work = (Matrix) Calloc(nsq, double);
    iMatrix worki = Calloc(nsq, int);
    vector revals = (vector) Calloc(n, double);
    vector ievals = (vector) Calloc(n, double);
    Matrix evecs = (Matrix) Calloc(nsq, double);
    Matrix evecsinv = (Matrix) Calloc(nsq, double);
    Matrix temp = (Matrix) Calloc(nsq, double);
    for (i=0; i<nsq; ++i)
	temp[i] = mat[i];
    /* calculate eigensystem */
    F77_CALL(rg) (&n, &n, temp, revals, ievals, &matz, evecs, worki, work, &err);
    if (repeated_entries (revals, n) || (err != 0)){
#if _MEXP_METHOD_==1 
	MatrixExpPade(expmat, mat, n, t);
#elif _MEXP_METHOD_==2 
	MatrixExpSeries(mat, n, expmat, t);
#endif
    }
    else {
	for (i=0; i<n; ++i)
	    revals[i] = exp(revals[i] * t);
#ifdef _USE_LAPACK_INVERSE_
	MatInv(evecs, evecsinv, n);
#endif
#ifdef _USE_LINPACK_INVERSE_
	MatInvDQR(evecs, evecsinv, n);
#endif
	MultMatDiag(revals, evecsinv, n, work);
	MultMat(evecs, work, n, n, n, expmat);
    }
    Free(work);  Free(worki);  Free(revals);  Free(ievals); Free(evecs);  Free(evecsinv); Free(temp);
}


/* Calls LAPACK routines to find exponential of a matrix by eigensystem decomposition */
/* If matrix has repeated eigenvalues, then use Pade approximation instead */

void MatrixExpDG(Matrix mat, int n, Matrix expmat, double t)
{
    Matrix work;
    iMatrix worki = Calloc(n*n, int);
    vector revals = (vector) Calloc(n, double);
    vector ievals = (vector) Calloc(n, double);
    Matrix evecs = (Matrix) Calloc(n*n, double);
    Matrix evecsinv = (Matrix) Calloc(n*n, double);
    Matrix temp = (Matrix) Calloc(n*n, double);
    int lwork = -1, i, err, nsq=n*n;
    char jobVL[1], jobVR[1];
    double *left=0, tmp;
    for (i=0; i<nsq; ++i)
	temp[i] = mat[i];
    jobVL[0] = 'N'; jobVR[0] = 'V';
    /* calculate optimal size of workspace */
    F77_CALL(dgeev)(jobVL, jobVR, &n, temp, &n, revals, ievals, left, &n, evecs, &n, &tmp, &lwork, &err);
    lwork = (int) tmp;
    work = (Matrix) Calloc(lwork, double);
    /* calculate eigensystem */
    F77_CALL(dgeev)(jobVL, jobVR, &n, temp, &n, revals, ievals, left, &n, evecs, &n, work, &lwork, &err);
    if (repeated_entries (revals, n) || (err != 0)){
#if _MEXP_METHOD_==1 
	MatrixExpPade(expmat, mat, n, t);
#elif _MEXP_METHOD_==2 
	MatrixExpSeries(mat, n, expmat, t);
#endif
    }
    else {
	for (i=0; i<n; ++i)
	    revals[i] = exp(revals[i] * t);
#ifdef _USE_LAPACK_INVERSE_
	MatInv(evecs, evecsinv, n);
#endif
#ifdef _USE_LINPACK_INVERSE_
	MatInvDQR(evecs, evecsinv, n);
#endif
	MultMatDiag(revals, evecsinv, n, work);
	MultMat(evecs, work, n, n, n, expmat);
    }
    Free(work);  Free(worki);  Free(revals);  Free(ievals);  Free(evecs);  Free(evecsinv);  Free(temp);
}

/* Fills the required entries of the intensity matrix with the current intensities */
/* qvector supplied filled by row. intens ordered by row, qmat filled by column */

void FillQmatrix(ivector qvector, vector intens, Matrix qmat, int nstates)
{
    int i, j, k=0;
    for (i=0; i<nstates; ++i) {
	qmat[MI(i, i, nstates)] = 0;
	for (j=0; j<nstates; ++j) {
	    if (j != i) {
		qmat[MI(i, j, nstates)] = 0;
		if (qvector[MI(j, i, nstates)] == 1) {
		    qmat[MI(i, j, nstates)] += intens[k];
		    qmat[MI(i, i, nstates)] -= intens[k];
		    ++k;
		}
	    }
	}
    }
}

/* Returns i-j transition intensity time t given vectors of intensities and transition indicators */

double qij(int i, int j, vector intens, ivector qvector, int nstates)
{
    double q; 
    Matrix qmat = Calloc( (nstates)*(nstates), double);
    FillQmatrix(qvector, intens, qmat, nstates);
    q = qmat[MI(i,j,nstates)];
    Free(qmat);
    return q;
}

/* Calculates the whole transition matrix in time t given an intensity matrix */

void Pmat(Matrix pmat, double t, vector intens, ivector qvector, int nstates, int exacttimes, int debug)
{
    int i,j;
    double pii;
    Matrix qmat = (Matrix) Calloc( (nstates)*(nstates), double);
    FillQmatrix(qvector, intens, qmat, nstates);
    if (exacttimes) { 
	for (i=0; i<nstates; ++i) 
	    for (j=0; j<nstates; ++j) {
		pii = exp(t * qmat[MI(i, i, nstates)] );
		pmat[MI(i, j, nstates)] = ( i==j  ?  pii  : pii * qmat[MI(i, j, nstates)] );
	    }
    }
    else {
#ifdef _USE_LINPACK_EIGEN_
	MatrixExp(qmat, nstates, pmat, t, debug);
#endif
#ifdef _USE_LAPACK_EIGEN_
	MatrixExpDG(qmat, nstates, pmat, t);
#endif
    }
    /* Floating point fuzz sometimes causes trouble */
    for (i=0; i<nstates; ++i) 
	for (j=0; j<nstates; ++j) {
	    if (pmat[MI(i, j, nstates)] < DBL_EPSILON) pmat[MI(i, j, nstates)] = 0;
	    if (pmat[MI(i, j, nstates)] > 1 - DBL_EPSILON) pmat[MI(i, j, nstates)] = 1;
	}
    Free(qmat);
}

/* Derivatives of P matrix - Todo. */

#if DERIV

void FormDQ(Matrix DQ, int u) 
{
    
}

void DPmat(Matrix pmat, double t, vector intens, ivector qvector, int nstates, int exacttimes)
{
    Matrix DQ = (Matrix) Calloc(nstates*nstates, double);
    vector revals = (vector) Calloc(n, double);
    vector ievals = (vector) Calloc(n, double);
    Matrix evecs = (Matrix) Calloc(n*n, double);
    Matrix evecsinv = (Matrix) Calloc(n*n, double); 
    int i, j; 
    FormDQ(DQ, u);
    MultMat(DQ, evecs, work);
    MultMat(evecsinv, work, G);
    for (i=0; i<n; ++i) {
	eit = exp(revals[i] * t);
	for (j=0; i<n; ++j) 
	    {
		if (i==j) 
		    V[MI(i, j, nstates)] = G[MI(i,i,nstates)] * t * eit;
		else { 
		    ejt = exp(revals[j] * t);
		    V[MI(i, j, nstates)] = G[MI(i,j,nstates)] * (eit - ejt) / (revals[i] - revals[j]);
		}
	    }    
    MultMat(V, evecsinv, work);
    MultMat(evecs, work, pmat);
    Free(DQ); Free(revals); Free(ievals); Free(evecs); Free(evecsinv);
}
#endif


