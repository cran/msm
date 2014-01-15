/* *****************************************************************
   PROGRAM: lik.c
   AUTHOR:  Chris Jackson
   DATE:    July 2004

   Routines for calculating likelihoods for multi-state Markov and
   hidden Markov models.

   ******************************************************************  */

#include "msm.h"
#include "hmm.h"
#include <Rmath.h>
#define NODEBUG
#define NODEBUG2
#define NOVITDEBUG
#define NOFISHDEBUG
#define NODERIVDEBUG

linkfn LINKFNS[3][2] = {
    {identity, identity},
    {log, exp},
    {logit, expit}
};

/* MUST KEEP THIS IN SAME ORDER AS .msm.HMODELPARS IN R/constants.R */
hmmfn HMODELS[] = {
    hmmCat,
    hmmIdent,
    hmmUnif,
    hmmNorm,
    hmmLNorm,
    hmmExp,
    hmmGamma,
    hmmWeibull,
    hmmPois,
    hmmBinom,
    hmmTNorm,
    hmmMETNorm,
    hmmMEUnif,
    hmmNBinom,
    hmmBeta,
    hmmT
};

/* MUST MATCH order of .msm.CTASKS in R/constants.R */
#define DO_LIK 0
#define DO_DERIV 1
#define DO_INFO 2
#define DO_VITERBI 3
#define DO_LIK_SUBJ 4
#define DO_DERIV_SUBJ 5
#define DO_DPMAT 6

#define OBS_SNAPSHOT 1
#define OBS_PANEL 1 /* preferred term now */
#define OBS_EXACT 2
#define OBS_DEATH 3

double logit(double x)
{
    return log(x / (1 - x));
}

double expit(double x)
{
    return exp(x) / ( 1 + exp(x) );
}

double identity(double x)
{
    return x;
}

/* Good-enough floating point equality comparison */

int all_equal(double x, double y)
{
    return fabs (x - y) <= DBL_EPSILON * fabs(x);
}

/* Return a vector of the nc possible true states that a censored state could represent */
/* These will be summed over when calculating the likelihood */
/* Compare one-indexed obs against one-indexed cm->censor. Return one-indexed current (*states) */

void GetCensored (double obs, cmodel *cm, int *nc, double **states)
{
    int j, k=0, n, cens=0;
    if (cm->ncens == 0)
	n = 1;
    else {
	while (!all_equal(obs, cm->censor[k]) && k < cm->ncens)
	    ++k;
	if (k < cm->ncens) {
	    cens = 1;
	    n =  cm->censstind[k+1] - cm->censstind[k];
	}
	else n = 1;
    }
    if (cm->ncens == 0 || !cens)
	(*states)[0] = obs;
    else { for (j = cm->censstind[k]; j < cm->censstind[k+1]; ++j)
	(*states)[j - cm->censstind[k]] = cm->censstates[j]; }
    *nc = n;
}

/*
  Calculate p (obs curr | true i) for hidden Markov models
  If observation is not necessarily of the true state (!obstrue),
  then this is just the HMM outcome probability (summed over censor set if necessary)
  If obstrue, this observation is not misclassified.
  e.g. censor set 1,2,3,    state set 1,2,3,4,
  pout =   if i in curr 1, else 0
*/

void GetOutcomeProb(double *pout, double *curr, int nc, double *hpars, hmodel *hm, qmodel *qm, int obstrue)
{
    int i, j;
    for (i=0; i<qm->nst; ++i) {
	pout[i] = 0;
	if (!obstrue) {
	    for (j=0; j<nc; ++j)
		pout[i] += (HMODELS[hm->models[i]])(curr[j], &(hpars[hm->firstpar[i]]));
	}
	else {
	    for (j=0; j<nc; ++j)
		if ((int) curr[j] == i+1)
		    pout[i] = 1;
	}
    }
}

void normalize(double *in, double *out, int n, double *lweight)
{
    int i; double ave;
    for (i=0, ave=0; i<n; ++i)
	ave += in[i];
    ave /= n;
    if (ave == 0) ave = 1;
    for (i=0; i<n; ++i)
	out[i] = in[i] / ave;
    *lweight -= log(ave);
}

/* Post-multiply the row-vector cump by matrix T to accumulate the likelihood */

void update_likhidden(double *curr, int nc, int obsno, msmdata *d, qmodel *qm,
		      hmodel *hm, double *cump, double *newp, double *lweight)
{
    int i, j, ideath;
    double *pout = Calloc(qm->nst, double);
    double *T           = Calloc((qm->nst)*(qm->nst), double);
    double *pmat        = Calloc((qm->nst)*(qm->nst), double);
    double *qmat = &(qm->intens[MI3(0, 0, obsno-1, qm->nst, qm->nst)]);
    double *hpars = &(hm->pars[MI(0, obsno, hm->totpars)]);

    GetOutcomeProb(pout, curr, nc, hpars, hm, qm, d->obstrue[obsno]);
#ifdef DEBUG2
    printf("hpars: ");
    for (i=0; i<hm->totpars; ++i)
	printf("%f,",hpars[i]);
    printf("\n");
    printf("pout: ");
    for (i = 0; i < qm->nst; ++i) {
	printf("%4.3f, ", pout[i]);
    }
    printf("\n: ");
#endif
    /* calculate the transition probability (P) matrix for the time interval dt */
    Pmat(pmat, d->time[obsno] - d->time[obsno-1], qmat, qm->nst,
	 (d->obstype[obsno] == OBS_EXACT), qm->iso, qm->perm, qm->qperm, qm->expm);
    for(j = 0; j < qm->nst; ++j)
	{
#ifdef DEBUG2
	    printf("pmat: ");
#endif
	    newp[j] = 0.0;
	    for(i = 0; i < qm->nst; ++i)
		{
 		    if (d->obstype[obsno] == OBS_DEATH) {
			/* Find the true state that obs represents.  This should be the state with outcome model hmmIdent(obs) */
 			if (d->obstrue[obsno])
 			    ideath = curr[0] - 1;
 			else
			    for (ideath=0; ideath < qm->nst; ++ideath)
				if (hm->models[ideath] == 1 && hmmIdent(curr[0], &(hpars[hm->firstpar[ideath]])))
				    break;
			T[MI(i,j,qm->nst)] = pmat[MI(i,j,qm->nst)] * qmat[MI(j,ideath,qm->nst)];
		    }
		    else {
			T[MI(i,j,qm->nst)] = pmat[MI(i, j, qm->nst)] * pout[j];
		    }
		    if (T[MI(i,j,qm->nst)] < 0) T[MI(i,j,qm->nst)] = 0;
		    newp[j] = newp[j] + cump[i]*T[MI(i,j,qm->nst)];
#ifdef DEBUG2
   		    printf("%4.3f, ", pmat[MI(i,j,qm->nst)]);
#endif
		}
#ifdef DEBUG2
   	    printf("\n");
	    printf("newp: ");
	    printf("%lf, ", newp[j]);
#endif
	}
    /* re-scale the likelihood at each step to prevent it getting too small and underflowing */
    /*  while cumulatively recording the log scale factor   */
    normalize (newp, cump, qm->nst, lweight);
    Free(pout); Free(T); Free(pmat);
}

/* Likelihood for the hidden Markov model for one individual */

double likhidden(int pt, /* ordinal subject ID */
		 msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm)
{
    double *curr = Calloc (qm->nst, double);
    /* no more than nst states allowed */
    double *cump     = Calloc(qm->nst, double);
    double *newp     = Calloc(qm->nst, double);
    double *pout = Calloc(qm->nst, double);
    double lweight, lik, *hpars;
    int i, obsno, nc=1;
    if (d->firstobs[pt] + 1 == d->firstobs[pt+1])
      return 0; /* individual has only one observation */
    /* Likelihood for individual's first observation */
    hpars = &(hm->pars[MI(0, d->firstobs[pt], hm->totpars)]);
    GetCensored((double)d->obs[d->firstobs[pt]], cm, &nc, &curr);
    GetOutcomeProb(pout, curr, nc, hpars, hm, qm, d->obstrue[d->firstobs[pt]]);
    /* Likelihood contribution for initial observation */
#ifdef DEBUG2
    printf("hpars: ");
    for (i=0; i<hm->totpars; ++i)
	printf("%f,",hpars[i]);
    printf("\n");
    for (i = 0; i < qm->nst; ++i) {
	printf("initp=%lf, ", hm->initp[MI(pt,i,d->npts)]);
	printf("pout=%lf, ", pout[i]);
    }
    printf("cump\n: ");
#endif
    for (i = 0; i < qm->nst; ++i) {
      cump[i] = pout[i];
      /* Ignore initprobs if observation is known to be the true state
	 or TODO, can we set it in R to one for obs state, zero for others? */
      if (!d->obstrue[d->firstobs[pt]]) cump[i] = cump[i]*hm->initp[MI(pt,i,d->npts)];
#ifdef DEBUG2
      printf("%lf, ", cump[i]);
#endif
    }
#ifdef DEBUG2
      printf("\n");
#endif
    lweight=0;
    /* Matrix product loop to accumulate the likelihood for subsequent observations */
    for (obsno = d->firstobs[pt]+1; obsno <= d->firstobs[pt+1] - 1; ++obsno)
	{
	    R_CheckUserInterrupt();
	    GetCensored((double)d->obs[obsno], cm, &nc, &curr);
	    update_likhidden(curr, nc, obsno, d, qm, hm, cump, newp, &lweight);
	}
#ifdef DEBUG2
      printf("cump: ");
#endif
    for (i = 0, lik = 0; i < qm->nst; ++i) {
#ifdef DEBUG2
      printf("%lf, ", cump[i]);
#endif
	lik = lik + cump[i];
    }
#ifdef DEBUG2
    printf("\nlik=%lf,lweight=%lf\n",lik,lweight);
#endif
      Free(curr); Free(cump);  Free(newp); Free(pout);
    /* Transform the likelihood back to the proper scale */
    return -2*(log(lik) - lweight);
}


void update_likcensor(int obsno, double *prev, double *curr, int np, int nc,
			msmdata *d, qmodel *qm,  hmodel *hm,
			double *cump, double *newp, double *lweight)
{
    double *pmat = Calloc((qm->nst)*(qm->nst), double);
    double *qmat = &(qm->intens[MI3(0, 0, obsno-1, qm->nst, qm->nst)]);
    double contrib;
    int i, j, k;
    Pmat(pmat, d->time[obsno] - d->time[obsno-1], qmat, qm->nst,
	 (d->obstype[obsno] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
    for(i = 0; i < nc; ++i)
	{
	    newp[i] = 0.0;
	    for(j = 0; j < np; ++j) {
		if (d->obstype[obsno] == OBS_DEATH) {
		    contrib = 0;
		    for (k = 0; k < qm->nst; ++k)
			if (k != curr[i]-1)
			    contrib += pmat[MI((int) prev[j]-1, k, qm->nst)] *
				qmat[MI(k, (int) curr[i]-1, qm->nst)];
		    newp[i] += cump[j] * contrib;

		}
		else {
		    newp[i] += cump[j] * pmat[MI((int) prev[j]-1, (int) curr[i]-1, qm->nst)];
		}
	    }
	}
    normalize(newp, cump, nc, lweight);
    Free(pmat);
}

double likcensor(int pt, /* ordinal subject ID */
		 msmdata *d, qmodel *qm,
		 cmodel *cm, hmodel *hm
		 )
{
    double *cump     = Calloc(qm->nst, double);
    double *newp     = Calloc(qm->nst, double);
    double *prev     = Calloc(qm->nst, double);
    double *curr     = Calloc(qm->nst, double);
    double lweight = 0, lik;
    int i, obs, np=0, nc=0;
    if (d->firstobs[pt] + 1 == d->firstobs[pt+1])
      return 0; /* individual has only one observation */
    for (i = 0; i < qm->nst; ++i)
	cump[i] = 1;
    GetCensored((double)d->obs[d->firstobs[pt]], cm, &np, &prev);
    for (obs = d->firstobs[pt]+1; obs <= d->firstobs[pt+1] - 1; ++obs)
	{
	    /* post-multiply by sub-matrix of P at each obs */
	    GetCensored((double)d->obs[obs], cm, &nc, &curr);
	    update_likcensor(obs, prev, curr, np, nc, d, qm, hm,
			     cump, newp, &lweight);
	    np = nc;
	    for (i=0; i<nc; ++i) prev[i] = curr[i];
	}
    for (i = 0, lik = 0; i < nc; ++i)
	lik = lik + cump[i];
    Free(cump);  Free(newp);  Free(prev); Free(curr);
    return -2*(log(lik) - lweight);
}

/* Likelihood for the non-hidden multi-state Markov model. Data of
   form "time-difference, covariates, from-state, to-state, number of
   occurrences" */

double liksimple(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm)
{
    int i;
    double lik=0, contrib=0;
    double *pmat = Calloc((qm->nst)*(qm->nst), double), *qmat;
    qmat = &(qm->intens[MI3(0, 0, 0, qm->nst, qm->nst)]);
    for (i=0; i < d->nobs; ++i)
	{
	    R_CheckUserInterrupt();
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstype[i] != d->obstype[i-1])) {
		/* we have a new timelag/covariates/obstype combination. Recalculate the
		   P matrix for this */
		/* pointer to Q matrix for ith datapoint */
		qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		Pmat(pmat, d->timelag[i], qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
	    }
	    if (d->obstype[i] == OBS_DEATH) {
		contrib = pijdeath(d->fromstate[i], d->tostate[i], pmat, qmat, qm->nst);
	    }
	    else
		contrib = pmat[MI(d->fromstate[i], d->tostate[i], qm->nst)];
	    lik += d->nocc[i] * log(contrib);
#ifdef DEBUG
/*	    printf("obs %d, from %d, to %d, time %lf, obstype %d, ", i, d->fromstate[i], d->tostate[i], d->timelag[i], d->obstype[i]);
	    printf("nocc %d, con %lf, lik %lf\n", d->nocc[i], log(contrib), lik);*/
	    printf("%d-%d in %lf, q=%lf,%lf, ll=%lf\n",d->fromstate[i], d->tostate[i], d->timelag[i],qmat[0],qmat[1], d->nocc[i] * log(contrib));
#endif
	}
    Free(pmat);
    return (-2*lik);
}

/* Likelihood for the non-hidden multi-state Markov model, by subject */

double liksimple_subj(int pt, /* ordinal subject ID */
		      msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm)
{
    int i, from, to;
    double lik=0, pm=0, dt;
    double *pmat = Calloc((qm->nst)*(qm->nst), double), *qmat;
    for (i = d->firstobs[pt]+1; i < d->firstobs[pt+1]; ++i) {
	R_CheckUserInterrupt();
	dt = d->time[i] - d->time[i-1];
	from = fprec(d->obs[i-1] - 1, 0); /* convert state outcome to integer */
	to = fprec(d->obs[i] - 1, 0);
	qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]); /* use covariate at start of transition */
	Pmat(pmat, dt, qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
	if (d->obstype[i] == OBS_DEATH)
	    pm = pijdeath(from, to, pmat, qmat, qm->nst);
	else
	    pm = pmat[MI(from, to, qm->nst)];
#ifdef DEBUG
	printf("i=%d, %d-%d in %lf, q=%lf,%lf, ll=%lf\n",i,from,to,dt,qmat[0],qmat[1],log(pm));
#endif
	lik += log(pm);
    }
    Free(pmat);
    return (-2*lik);
}

void msmLikelihood (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    int pt;
    double likone;
    *returned = 0;
    /* Likelihood for hidden Markov model */
    if (hm->hidden)
	{
	    for (pt = 0;  pt < d->npts; ++pt){
		likone = likhidden (pt, d, qm, cm, hm);
#ifdef DEBUG
		printf("pt %d, lik %lf\n", pt, likone);
#endif
		*returned += likone;
	    }
	}
    /* Likelihood for Markov model with censored outcomes */
    else if (cm->ncens > 0)
	{
	    for (pt = 0;  pt < d->npts; ++pt){
		likone = likcensor (pt, d, qm, cm, hm);
		*returned += likone;
	    }
	}
    /* Likelihood for simple non-hidden, non-censored Markov model */
    else {
	*returned = liksimple (d, qm, cm, hm);
    }
}

/* First derivatives of the log-likelihood for the non-hidden
   multi-state Markov model. */

void derivsimple(msmdata *d, qmodel *qm,  cmodel *cm, hmodel *hm, double *deriv)
{
    int i, p, np=qm->ncpars;
    double *qmat, *dqmat;
    double *pmat = Calloc(qm->nst * qm->nst, double);
    double *dpmat = Calloc(qm->nst * qm->nst * np, double);
    double *dp = Calloc(np, double);
    double pm;
    qmat = &(qm->intens[MI3(0, 0, 0, qm->nst, qm->nst)]);
    dqmat = &(qm->dintens[MI4(0, 0, 0, 0, qm->nst, qm->nst, np)]);
    for (p = 0; p < np; ++p) deriv[p] = 0;
    for (i=0; i < d->nobs; ++i)
    {
	    R_CheckUserInterrupt();
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstype[i] != d->obstype[i-1])) {
		/* we have a new timelag/covariates/obstype combination. Recalculate the
		   P matrix and its derivatives for this */
		qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		Pmat(pmat, d->timelag[i], qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
		dqmat = &(qm->dintens[MI4(0, 0, 0, i, qm->nst, qm->nst, np)]);
 		DPmat(dpmat, d->timelag[i], dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));
	    }
	    if (d->obstype[i] == OBS_DEATH) {
		pm = pijdeath(d->fromstate[i], d->tostate[i], pmat, qmat, qm->nst);
		dpijdeath(d->fromstate[i], d->tostate[i], dpmat, pmat, dqmat, qmat, qm->nst, np, dp);
	    }
	    else {
		pm = pmat[MI(d->fromstate[i], d->tostate[i], qm->nst)];
		for (p = 0; p < np; ++p)
		    dp[p] = dpmat[MI3(d->fromstate[i], d->tostate[i], p, qm->nst, qm->nst)];
	    }
	    for (p = 0; p < np; ++p) {
		if (pm > 0)
		    deriv[p] += d->nocc[i] * dp[p] / pm;
#ifdef DERIVDEBUG
		printf("%d, %d, %d, %d, %6.4f, %d, %d, %lf, %lf\n", i, p, d->fromstate[i], d->tostate[i], d->timelag[i], d->obstype[i], d->nocc[i], dp[p], -2 * d->nocc[i] * dp[p]/pm);
#endif
	    }
    }
    for (p = 0; p < np; ++p) {
	deriv[p] *= -2;   /* above is dlogL/dtu as in Kalb+Law, we want
			     deriv of -2 loglik  */
    }
    Free(pmat); Free(dpmat); Free(dp);
}

/* First derivatives of the likelihood for the non-hidden multi-state
   Markov model. Uses data of form "subject, time, state, covariates".
   Returns derivatives by individual and parameter for use in the
   score residual diagnostic.
*/

void derivsimple_subj(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *deriv)
{
    int pt, i, p, np=qm->ncpars;
    double *qmat, *dqmat;
    double *pmat = Calloc(qm->nst * qm->nst, double);
    double *dpmat = Calloc(qm->nst * qm->nst * np, double);
    double *dp = Calloc(np, double);
    double pm=0, dt;
    int from, to;
    for (pt = 0;  pt < d->npts; ++pt){
	{
	    R_CheckUserInterrupt();
	    if (d->firstobs[pt+1] > d->firstobs[pt] + 1) { /* individual has more than one observation? */
	      for (p = 0; p < np; ++p) {
		deriv[MI(pt,p,d->npts)] = 0;
	      }
	      for (i = d->firstobs[pt]+1; i < d->firstobs[pt+1]; ++i) {
		    dt = d->time[i] - d->time[i-1];
		    from = fprec(d->obs[i-1] - 1, 0); /* convert state outcome to integer */
		    to = fprec(d->obs[i] - 1, 0);
		    qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]); /* use covariate at start of transition */
		    Pmat(pmat, dt, qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
		    dqmat = &(qm->dintens[MI4(0, 0, 0, i-1, qm->nst, qm->nst, np)]);
		    DPmat(dpmat, dt, dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));
		    if (d->obstype[i] == OBS_DEATH) {
			pm = pijdeath(from, to, pmat, qmat, qm->nst);
			dpijdeath(from, to, dpmat, pmat, dqmat, qmat, qm->nst, np, dp);
		    }
		    else {
			pm = pmat[MI(from, to, qm->nst)];
			for (p = 0; p < np; ++p)
			    dp[p] = dpmat[MI3(from, to, p, qm->nst, qm->nst)];
		    }
		    for (p = 0; p < np; ++p) {
		      deriv[MI(pt,p,d->npts)] += dp[p] / pm; /* on loglik scale not -2*loglik */
#ifdef DEBUG
		      printf("%d, %d, %d, %d, %d, %6.4f, %d, %lf, %lf\n", i, p, pt, from, to, dt, d->obstype[i], -2 * dp[p] / pm, -2*deriv[MI(pt,p,d->npts)]);
#endif
		    }
		}
		for (p = 0; p < np; ++p)
		    deriv[MI(pt,p,d->npts)] *= -2;
	    }
	    else for (p = 0; p < np; ++p)
		deriv[MI(pt,p,d->npts)] = 0;
	}
    }
    Free(pmat); Free(dpmat); Free(dp);
}

void infosimple(msmdata *d, qmodel *qm,  cmodel *cm, hmodel *hm, double *info)
{
    int i, j, p, q, np=qm->ncpars;
    double *qmat, *dqmat;
    double *pmat = Calloc(qm->nst * qm->nst, double);
    double *dpmat = Calloc(qm->nst * qm->nst * np, double);
    double *dpm = Calloc(qm->nst* np, double);
    double *pm = Calloc(qm->nst, double);
    for (i=0; i < d->nobs; ++i)
	{
	    R_CheckUserInterrupt();
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstype[i] != d->obstype[i-1])) {
		/* we have a new timelag/covariates/obstype combination. Recalculate the
		   P matrix and its derivatives for this */
		qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		Pmat(pmat, d->timelag[i], qmat, qm->nst, (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);
		dqmat = &(qm->dintens[MI4(0, 0, 0, i, qm->nst, qm->nst, np)]);
 		DPmat(dpmat, d->timelag[i], dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));
	    }
	    if (d->obstype[i] != OBS_PANEL)
		error("Fisher information only available for panel data\n");
	    for (j=0; j < qm->nst; ++j)  {
		pm[j] = pmat[MI(d->fromstate[i], j, qm->nst)];
		for (p = 0; p < np; ++p)
		    dpm[MI(j,p,qm->nst)] = dpmat[MI3(d->fromstate[i], j, p, qm->nst, qm->nst)];
	    }
	    if ((i==0) || (d->whicha[i] != d->whicha[i-1]) || (d->obstype[i] != d->obstype[i-1]) ||
		(d->fromstate[i] != d->fromstate[i-1])) {
		for (p = 0; p < np; ++p) {
		    for (q = 0; q < np; ++q) {
			/* for expected information, sum over all possible destination states for this fromstate/time/cov combination */
			for(j = 0; j<qm->nst; ++j)  {
			    if (pm[j] > 0)
				info[MI(p,q,np)] +=  d->noccsum[i] * dpm[MI(j,p,qm->nst)] * dpm[MI(j,q,qm->nst)] /  pm[j];
			}
		    }
		}
	    }
	}
    for (p = 0; p < np; ++p) {
	for (q = 0; q < np; ++q) {
	    info[MI(p,q,np)] *= 2; /* above is E(-d2logL/dtudtv) as in Kalb+Law.
					   we want second deriv of of -2 loglik */
	}
    }
    Free(pm); Free(dpm); Free(dpmat); Free(pmat);
}

/* Derivatives of the P matrix, used for the Pearson test p-value */
/* Returns a ntrans * ntostates * npars matrix */
/* Panel data only (obstype 1) */

void dpmat_obs(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *deriv)
{
    int pt, i, j, k, p, from, np = qm->ncpars;
    double *dpmat = Calloc(qm->nst * qm->nst * np, double), *qmat, *dqmat;
    double dt;

    j=0;
    for (pt = 0;  pt < d->npts; ++pt)
	{
	    R_CheckUserInterrupt();
	    if (d->firstobs[pt+1] > d->firstobs[pt] + 1) { /* individual has more than one observation? */
		for (i = d->firstobs[pt]+1; i < d->firstobs[pt+1]; ++i) {
		    ++j;
		    qmat = &(qm->intens[MI3(0, 0, i, qm->nst, qm->nst)]);
		    dqmat = &(qm->dintens[MI4(0, 0, 0, i, qm->nst, qm->nst, np)]);
		    dt = d->time[i] - d->time[i-1];
		    from = fprec(d->obs[i-1] - 1, 0); /* convert state outcome to integer */
		    DPmat(dpmat, dt, dqmat, qmat, qm->nst, np, (d->obstype[i] == OBS_EXACT));
		    for (p = 0; p < np; ++p) {
			for (k=0; k < qm->nst; ++k) {
			    deriv[MI3(j-1,k,p,d->ntrans,qm->nst)] = dpmat[MI3(from, k, p, qm->nst, qm->nst)];
//			    printf("%d %lf, ",MI3(j-1,k,p,d->ntrans,qm->nst),dpmat[MI3(from, k, p, qm->nst, qm->nst)]);
			}
		    }
//		    printf("\n");
		}
	    }
	}
    Free(dpmat);
}

/* Return vector of subject-specific log likelihoods */

void msmLikelihood_subj (msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *returned)
{
    int pt;
    for (pt = 0;  pt < d->npts; ++pt){
	if (hm->hidden)
	    returned[pt] = likhidden (pt, d, qm, cm, hm);
	else if (cm->ncens > 0)
	    returned[pt] = likcensor (pt, d, qm, cm, hm);
	else
	    returned[pt] = liksimple_subj (pt, d, qm, cm, hm);
    }
}

/* Find zero-based index of maximum element of a vector x */

void pmax(double *x, int n, int *maxi)
{
    int i=0;
    *maxi = i;
    for (i=1; i<n; ++i) {
	if (x[i] > x[*maxi]) {
	    *maxi = i;
	}
    }
}

/* Calculates the most likely path through underlying states */

void Viterbi(msmdata *d, qmodel *qm, cmodel *cm, hmodel *hm, double *fitted)
{
    int i, j, tru, k, kmax, obs, nc = 1;
/*
    double *pmat = (double *) S_alloc((qm->nst)*(qm->nst), sizeof(double));
    int *ptr = (int *) S_alloc((d->nobs)*(qm->nst), sizeof(int));
    double *lvold = (double *) S_alloc(qm->nst, sizeof(double));
    double *lvnew = (double *) S_alloc(qm->nst, sizeof(double));
    double *lvp = (double *) S_alloc(qm->nst, sizeof(double));
    double *curr = (double *) S_alloc (qm->nst, sizeof(double));
    double *pout = (double *) S_alloc(qm->nst, sizeof(double));
*/

    double *pmat = Calloc((qm->nst)*(qm->nst), double);
    int *ptr = Calloc((d->nobs)*(qm->nst), int);
    double *lvold = Calloc(qm->nst, double);
    double *lvnew = Calloc(qm->nst, double);
    double *lvp = Calloc(qm->nst, double);
    double *curr = Calloc (qm->nst, double);
    double *pout = Calloc(qm->nst, double);

    double dt;
    double *qmat, *hpars;

    i = 0;
    if (d->obstrue[i]) {
      for (k = 0; k < qm->nst; ++k)
	lvold[k] = (k+1 == d->obs[i] ? 0 : R_NegInf);
    }
    else {
      GetCensored(d->obs[i], cm, &nc, &curr);
      /* initial observation is a censored state. No HMM here, so initprobs not needed */
      if (nc > 1) {
	for (k = 0, j = 0; k < qm->nst; ++k) {
	  if (k+1 == curr[j]) {
	    lvold[k] = 0;
	    ++j;
	  }
	  else lvold[k] = R_NegInf;
	}
      }
      else { /* use initprobs */
	for (k = 0; k < qm->nst; ++k)
	    lvold[k] = log(hm->initp[MI(0, k, d->npts)]);
      }
    }

    for (i = 1; i <= d->nobs; ++i)
	{
	    R_CheckUserInterrupt();
#ifdef VITDEBUG
	    printf("obs %d\n", i);
#endif
	    if ((i < d->nobs) && (d->subject[i] == d->subject[i-1]))
		{
#ifdef VITDEBUG
		    printf("subject %d\n", d->subject[i]);
#endif
		    dt = d->time[i] - d->time[i-1];
		    qmat = &(qm->intens[MI3(0, 0, i-1, qm->nst, qm->nst)]);
		    hpars = &(hm->pars[MI(0, i, hm->totpars)]); /* not i-1 as pre 1.2.3 */
		    GetCensored(d->obs[i], cm, &nc, &curr);
		    GetOutcomeProb(pout, curr, nc, hpars, hm, qm, d->obstrue[i]);
#ifdef VITDEBUG
		    for (tru=0;tru<nc;++tru) printf("curr[%d] = %1.0lf, ",tru, curr[tru]); printf("\n");
#endif
		    Pmat(pmat, dt, qmat, qm->nst,
			 (d->obstype[i] == OBS_EXACT), qm->iso, qm->perm,  qm->qperm, qm->expm);

		    for (tru = 0; tru < qm->nst; ++tru)
			{
/* lvnew =  log prob of most likely path ending in tru at current obs.
   kmax  = most likely state at the previous obs
*/
			    for (k = 0; k < qm->nst; ++k) {
				lvp[k] = lvold[k] + log(pmat[MI(k, tru, qm->nst)]);
			    }
			    if (d->obstrue[i-1])
				kmax = d->obs[i-1] - 1;
			    else pmax(lvp, qm->nst, &kmax);
			    lvnew[tru] = log ( pout[tru] )  +  lvp[kmax];
			    ptr[MI(i, tru, d->nobs)] = kmax;
#ifdef VITDEBUG
			    printf("true %d, pout[%d] = %lf, lvold = %lf, pmat = %lf, lvnew = %lf, ptr[%d,%d]=%d\n",
				   tru, tru, pout[tru], lvold[tru], pmat[MI(kmax, tru, qm->nst)], lvnew[tru], i, tru, ptr[MI(i, tru, d->nobs)]);
#endif
			}
		    for (k = 0; k < qm->nst; ++k)
			lvold[k] = lvnew[k];
		}
	    else
		{
		    /* Traceback for current individual */
		    pmax(lvold, qm->nst, &kmax);
		    obs = i-1;
		    fitted[obs] = (d->obstrue[obs] ? d->obs[obs]-1 : kmax);
#ifdef VITDEBUG
		    printf("traceback for subject %d\n", d->subject[i-1]);
		    printf("obs=%d,fitted[%d]=%1.0lf\n",obs,obs,fitted[obs]);
#endif
		    while   ( (obs > 0) && (d->subject[obs] == d->subject[obs-1]) )
			{
			    fitted[obs-1] = ptr[MI(obs, fitted[obs], d->nobs)];
#ifdef VITDEBUG
			    printf("fitted[%d] = ptr[%d,%1.0lf] = %1.0lf\n", obs-1, obs, fitted[obs], fitted[obs-1]);
#endif
			    --obs;
			}
#ifdef VITDEBUG
		    printf("\n");
#endif
		    if (i < d->nobs) {
			if (d->obstrue[i]) {
			    for (k = 0; k < qm->nst; ++k)
				lvold[k] = (k+1 == d->obs[i] ? 0 : R_NegInf);
			}
			else {
			    GetCensored(d->obs[i], cm, &nc, &curr);
			    /* initial observation is a censored state. No HMM here, so initprobs not needed */
			    if (nc > 1) {
				for (k = 0, j = 0; k < qm->nst; ++k) {
				    if (k+1 == curr[j]) {
					lvold[k] = 0;
					++j;
				    }
				    else lvold[k] = R_NegInf;
				}
			    }
			    else { /* use initprobs */
				for (k = 0; k < qm->nst; ++k)
				    lvold[k] = log(hm->initp[MI(d->subject[i]-1, k, d->npts)]);
			    }
			}
		    }
		}
	}
    Free(pmat); Free(ptr); Free(lvold); Free(lvnew); Free(lvp); Free(curr); Free(pout);
}

/* This function is called from R to provide an entry into C code for
   evaluating likelihoods and related quantities  */

void msmCEntry(
	       int *do_what,
	       /* Parameters */
	       double *Q,  /* Array nstates x nstates x nobs */
	       double *DQ,  /* Array nstates x nstates x npars x nobs */
	       double *hmmpars,

	       /* Data for non-HMM multi-state models */
	       int *fromstate,  /* Distinct combinations of from, to, time lag, covariates */
	       int *tostate,
	       double *timelag,
	       int *nocc,       /* Number of occurrences of each distinct combination */
	       int *noccsum,    /* Number of those summed over tostate and replicated (for Fisher info) */
	       int *whicha,   /* indicator for the from, to, time lag, covs combination corresponding to the current obs */
	       int *obstype,   /* observation scheme, 1 snapshot, 2 exact, 3 death */

	       /* Data for HMMs */
	       int *subjvec,      /* vector of subject IDs */
	       double *timevec,   /* vector of observation times */
	       double *obsvec,     /* vector of observations (observed states in the case of misclassification models) */
	       int *firstobs,   /* 0-based index into data of the first observation for each individual */
	       int *obstrue,   /* which observations in a HMM represent the true underlying state (none by default) */

	       /* HMM specification */
	       int *hidden,       /* well, is it or isn't it? */
	       int *hmodels,      /* which hidden Markov distribution */
	       int *htotpars,      /* total number of HMM parameters */
	       int *hfirstpar,    /* index into hmmpars of first parameter for each state */
	       double *initprobs, /* initial state occupancy probabilities, as a npts x nstates matrix (since v1.2.1) */

	       int *nst,      /* number of Markov states */
	       int *iso, /* graph isomorphism ID */
	       int *perm, /* permutation to the base isomorphism */
	       int *qperm, /* permutation of intensity parameters from the base isomorphism */
	       int *expm, /* use expm package for matrix exponential */
	       int *nintens,      /* number of intensity parameters */
	       int *nobs,         /* number of observations in data set (hmm/cens) or number of aggregated transitions (standard) */
	       int *n,         /* number of observations in data set (used for derivs by individual in standard models. no of rows of covobsvec) */
	       int *npts,         /* number of individuals in data set */
	       int *ntrans,       /* number of (disaggregated) transitions, equal to n - npts if at least two obs (one trans) per individual */
	       int *ncoveffs,     /* number of covariate effects on transition rates (excluding intercept, used for dim of derivatives) */

	       int *ncens,     /* number of distinct forms of censoring */
	       int *censor,    /* censoring indicators in data, of length ncens */
	       int *censstates, /* list of possible states represented by censoring indicators */
	       int *censstind,  /* starting index into censstates for each censoring indicator */

	       double *returned,   /* returned -2 log likelihood , derivatives or Viterbi fitted values */
	       double *returned2   /* Fisher information if required */
	       )
{
    msmdata d; qmodel qm; cmodel cm; hmodel hm;

    d.fromstate = fromstate;     d.tostate = tostate;     d.timelag = timelag;
    d.nocc = nocc;  d.noccsum=noccsum;
    d.whicha = whicha; d.obstype = obstype; d.obstrue = obstrue;
    d.subject = subjvec; d.time = timevec; d.obs = obsvec; d.firstobs = firstobs;
    d.nobs = *nobs;  d.n = *n; d.npts = *npts; d.ntrans = *ntrans;

    qm.nst = *nst; qm.npars = *nintens; qm.intens = Q;
    qm.dintens = DQ; qm.ncpars = *nintens + *ncoveffs;
    qm.iso = *iso; qm.perm = perm; qm.qperm = qperm; qm.expm = *expm;

    cm.ncens = *ncens; cm.censor = censor; cm.censstates=censstates; cm.censstind=censstind;

    hm.hidden = *hidden;  hm.models = hmodels;  hm.totpars = *htotpars;
    hm.firstpar = hfirstpar;  hm.pars = hmmpars;
    hm.initp = initprobs;
    if (*do_what == DO_LIK) {
	msmLikelihood(&d, &qm, &cm, &hm, returned);
    }

    else if (*do_what == DO_DERIV) {
	derivsimple(&d, &qm, &cm, &hm, returned);
    }

    else if (*do_what == DO_INFO) {
	derivsimple(&d, &qm, &cm, &hm, returned);
	infosimple(&d, &qm, &cm, &hm, returned2);
    }

    else if (*do_what == DO_LIK_SUBJ) {
	msmLikelihood_subj(&d, &qm, &cm, &hm, returned);
    }

    else if (*do_what == DO_DERIV_SUBJ) {
	derivsimple_subj(&d, &qm, &cm, &hm, returned); /* derivative of loglik by subject. used for score residuals */
    }

    else if (*do_what == DO_VITERBI) {
	Viterbi(&d, &qm, &cm, &hm, returned);
    }

    else if (*do_what == DO_DPMAT) {
	dpmat_obs(&d, &qm, &cm, &hm, returned);
    }

    else error("Unknown C task.\n");
}
