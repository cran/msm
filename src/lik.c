/* *****************************************************************
   PROGRAM: lik.c 
   AUTHOR:  Chris Jackson
   DATE:    November 2001

   Routines for calculating likelihoods for Markov models with and without misclassification

   ******************************************************************  */ 


#include "lik.h"

/* This function is called from R to provide an entry into C code for
   evaluating likelihoods, doing Viterbi state reconstruction and short
   term prediction */

void msmCEntry( 
    int *do_what,      /* 1 = eval likelihood, 2 = Viterbi */
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
    int *constraint,    /* list of constraints for each covariate */
    double *misccovvec,/* vectorised matrix of misclassification covariate values */
    int *miscconstraint,/* list of constraints for each misclassification covariate */
    double *initprobs, /* initial state occupancy probabilities */
    int *nst,      /* number of Markov states */
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
    double *returned   /* returned -2 log likelihood , Viterbi fitted values, or predicted values */
    )
{
    int ifix = 0, iopt = 0, iall = 0;
    double *intens = (double *) S_alloc(*nintens, sizeof(double));
    double *coveffect = (double *) S_alloc(*ncoveffs, sizeof(double));
    double *miscprobs = (double *) S_alloc(*nmisc, sizeof(double));
    double *misccoveffect = (double *) S_alloc(*nmisccoveffs, sizeof(double));
    data d; 
    model m;

    fillparvec ( intens,        params, allinits, fixedpars, *nintens,      &ifix, &iopt, &iall ) ;
    fillparvec ( coveffect,     params, allinits, fixedpars, *ncoveffs,     &ifix, &iopt, &iall ) ;
    fillparvec ( miscprobs,     params, allinits, fixedpars, *nmisc,        &ifix, &iopt, &iall ) ;
    fillparvec ( misccoveffect, params, allinits, fixedpars, *nmisccoveffs, &ifix, &iopt, &iall ) ;

    d.subject = subjvec;     d.time = timevec;     d.state = statevec;     d.tostate = tostatevec;
    d.cov = covvec;     d.misccov = misccovvec;    d.fromto = *fromto;   d.nobs = *nobs;   
    d.npts = *npts;    d.ncovs = *ncovs;    d.nmisccovs = *nmisccovs;    d.tunit = *tunit; 

    m.qvector = qvector;    m.evector = evector;    m.constraint = constraint;    m.miscconstraint = miscconstraint;
    m.nst = *nst;  m.nms = *nms;  m.nintens = *nintens;  m.nmisc = *nmisc;  m.ncoveffs = *ncoveffs;
    m.nmisccoveffs = *nmisccoveffs;  m.covmatch = *covmatch;  m.death = *death;  m.exacttimes = *exacttimes;
    m.intens = intens;  m.coveffect = coveffect;  m.miscprobs = miscprobs;  m.misccoveffect = misccoveffect;
    m.initprobs = initprobs;

    if (*do_what == 1) {
	msmLikelihood(&d, &m, *misc, returned);
    }

    else if (*do_what == 2) {
	Viterbi(&d, &m, returned);
    }

}

void msmLikelihood (data *d, model *m, int misc, double *returned) {
    int pt;
    double likone;
    /* Likelihood for misclassification model */
    if (misc) 
    {
	*returned = 0;
	for (pt = 0;  pt < d->npts; ++pt){
	    likone = likmisc (pt, d, m);
	    *returned += likone;
	}
    }
  
    /* Likelihood for simple model without misclassification */
    else 
    {
	if (d->fromto)
	    *returned = liksimple_fromto (d, m);				     
	else
	    *returned = liksimple (d, m);
    }
}


/* Fill a parameter vector with either the current values from the optimisation or the fixed inital values */

void fillparvec(double *parvec, /* named vector to fill (e.g. intens = baseline intensities) */
		double *params, /* current values of parameters being optimised over */
		double *allinits, /* full vector of initial values */
		int *fixi,  /* indices of allinits which are fixed */
		int ni,    /* length of parvec */
		int *ifix,  /* current index into fixi */
		int *iopt,  /* current index into params */
		int *iall   /* current index into allinits */
    )
{
    int i;
    for (i=0; i<ni; ++i, ++(*iall)){
	if (*iall == fixi[*ifix]){
	    parvec[i] = allinits[*iall];
	    ++(*ifix);
	}
	else {
	    parvec[i] = params[*iopt];
	    ++(*iopt);
	}
    }
}



/* Likelihood for the misclassification model for one individual */

double likmisc(int pt, /* ordinal subject ID */
	       data *d, model *m
    )
{
    double *cumprod     = (double *) S_alloc(m->nst, sizeof(double)); 
    double *newprod     = (double *) S_alloc(m->nst, sizeof(double));  
    double *newmisc     = (double *) S_alloc(m->nmisc, sizeof(double));
    double lweight, lik;
    int i, pti, k, first=0, last=0;
    for (i = 1, pti = 0; i < d->nobs ; i++)
    {  /* find index in data set of individual's first and last observations */
	if (pti==pt) {
	    first = i-1;
	    while ((d->subject[i] == d->subject[i-1]) && (i < d->nobs)) ++i;
	    last = i-1;
	    break;
	}
	else if (d->subject[i] != d->subject[i-1]) ++pti;
    }
    AddMiscCovs(first, d, m, newmisc); 
    for (i = 0; i < m->nst; ++i)
	cumprod[i] = PObsTrue(d->state[first], i, newmisc, m) * m->initprobs[i]; /* cumulative matrix product */
    lweight=0;

    /* Matrix product loop to accumulate the likelihood */
    for (k = first+1; k <= last; ++k)
    {
	UpdateLik(d->state[k], d->time[k] - d->time[k-1],
		  k, last, 0, d, m, cumprod, newprod, lweight, &lweight);
	for (i = 0; i < m->nst; ++i)
	    cumprod[i] = newprod[i];
    }
  
    lik=0;
    for (i = 0; i < m->nst; ++i)
	lik = lik + cumprod[i];
    /* Transform the likelihood back to the proper scale */
    return -2*(log(lik) - lweight); 
}


/* Post-multiply the row-vector cumprod by matrix T to accumulate the likelihood */

void UpdateLik(int state, double dt, int k, int last, int predict_death, data *d, model *m, 
	       double *cumprod, double *newprod, double lweight_old, double *lweight_new)
{
    int i, j;
    double newprod_ave;
    double *T           = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *newintens   = (double *) S_alloc(m->nintens, sizeof(double) );
    double *newmisc     = (double *) S_alloc(m->nmisc, sizeof(double));
    double *pmat        = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *pmatdeath   = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *pmattodeath = (double *) S_alloc((m->nst)*(m->nst), sizeof(double)); 
    AddCovs(k - 1 + m->covmatch, d, m, newintens);
    AddMiscCovs(k - 1 + m->covmatch, d, m, newmisc);
    /* calculate the transition probability (P) matrix for the time interval dt */
    Pmat(pmat, dt, newintens, m->qvector, m->nst, m->exacttimes);
    if ((state == m->nst - 1) && (m->death) && (!m->exacttimes))
    {
	/* or dt - 1 day and 1 day, if death time is known within a day */
	Pmat(pmattodeath, dt - 1 / d->tunit, newintens, m->qvector, m->nst, m->exacttimes);
	Pmat(pmatdeath, 1/d->tunit, newintens, m->qvector, m->nst, m->exacttimes);
    }
    for(j = 0; j < m->nst; ++j)
    {
	newprod[j] = 0.0;
	for(i = 0; i < m->nst; ++i)
	{
	    if (predict_death)
		/* Only for doing prediction - used to calculate probability of death by t conditionally on history */
		T[MI(i,j,m->nst)] = (j==0  ?  pmat[MI(i, m->nst - 1, m->nst)]  :  0);
	    else if ((k == last) && (state == m->nst - 1) && (m->death))
		/* last observation was death and death time known exactly */
		T[MI(i,j,m->nst)] = pmattodeath[MI(i,j,m->nst)] * pmatdeath[MI(j, m->nst - 1, m->nst)];
	    else
		T[MI(i,j,m->nst)] = pmat[MI(i, j, m->nst)] * PObsTrue(state, j, newmisc, m);
	    if (T[MI(i,j,m->nst)] < 0) T[MI(i,j,m->nst)] = 0;
	    newprod[j] = newprod[j] + cumprod[i]*T[MI(i,j,m->nst)];
	}
    }
    /* re-scale the likelihood at each step to prevent it getting too small and underflowing */
    /*  while cumulatively recording the log scale factor   */
    for(i = 0, newprod_ave=0.0; i < m->nst ; ++i)
	newprod_ave += newprod[i];  
    newprod_ave /=  m->nst;
    if (newprod_ave == 0)
	newprod_ave = 1;
    for(i = 0; i < m->nst ; ++i)
	newprod[i] = newprod[i] / newprod_ave;
    *lweight_new = lweight_old  -  log(newprod_ave);
}



/* Add effects of covariates to the transition rates at an observation number obs */

void AddCovs(int obs, data *d, model *m, double *newintens)
{
    double *lintens =  (double *) S_alloc ( m->nintens, sizeof(double) );
    int j, k;
    for (j = 0; j < m->nintens; ++j)
    { 
	lintens[j] = log(m->intens[j]);
	for (k = 0; k < d->ncovs; ++k)
	    lintens[j] += m->coveffect[m->constraint[MI(k, j, m->nintens)] - 1] * d->cov[MI(d->nobs, obs, k)];
	newintens[j] = exp(lintens[j]);
    }
}


/* Add effects of misclassification covariates to the misclassification rates at an observation number obs */

void AddMiscCovs(int obs, data *d, model *m, double *newp)
{
    double *logitp = (double *) S_alloc ( m->nmisc , sizeof(double) );
    int j, k;
    for (j = 0; j < m->nmisc; ++j)
    { 
	logitp[j] = logit(m->miscprobs[j]);
	for (k = 0; k < d->nmisccovs; ++k)
	    logitp[j] += m->misccoveffect[m->miscconstraint[MI(k, j, m->nmisc)] - 1] * d->misccov[MI(d->nobs, obs, k)];
	newp[j] = expit(logitp[j]);
    }
}


/* Misclassification model: 
   calculate P(obs state | true state, covariates) */

double PObsTrue(int obst,      /* observed state */
		int tst,       /* true state */
		double *miscprobs, /* misclassification probabilities */
		model *m
    )
{
    int s;
    double this_miscprob, probsum;
    double *emat = (double *) S_alloc( (m->nst)*(m->nst), sizeof(double));
    /* construct the misclassification prob matrix from the parameter vector */
    FillQmatrix(m->evector, miscprobs, emat, m->nst); /* extend this so they sum to 1 instead of 0 ? */
    if (obst == tst){
	probsum = 0;
	for (s = 0; s < m->nst; ++s)
	    if (s != tst)
		probsum += emat[MI(tst, s, m->nst)];
	this_miscprob = 1 - probsum;
    }
    else 
	this_miscprob = emat[MI(tst, obst, m->nst)];
    return this_miscprob;
}


/* Likelihood for the simple model. Data of form "subject ID, obs time, obs state" */

double liksimple(data *d, model *m)
{
    int i,j;
    double dt, lik=0, contrib, *newintens = (double *) S_alloc ( m->nintens , sizeof(double) );
    for (i = 1; i < d->nobs; ++i){
	AddCovs(i - 1 + m->covmatch, d, m, newintens);
	if (d->subject[i-1] == d->subject[i]){ 
	    dt = d->time[i] - d->time[i-1];
	    if ((m->death) && (d->state[i] == m->nst - 1) && (!m->exacttimes)) 
	    {
		/* if final state is death, then death date is known within a day */
		contrib=0;
		for (j = d->state[i-1]; j < d->state[i]; ++j)
		  contrib += 
		      pijt(d->state[i-1], j,  dt - 1 / d->tunit, newintens, m->qvector, m->nst, m->exacttimes)*
		      pijt(j, d->state[i], 1 / d->tunit,  newintens, m->qvector, m->nst, m->exacttimes);
		lik += log(contrib);
	    }
	    else 
		lik += log(pijt(d->state[i-1], d->state[i], dt, newintens, m->qvector, m->nst, m->exacttimes));
	    /*	    printf("%d, lik=%lf\n", d->subject[i], lik); */
	}
    } 
    return (-2*lik); 
}

/* Likelihood for the simple model. Data of form "from-state, to-state, time-difference" */

double liksimple_fromto(data *d, model *m)
{
    int i,j;
    double lik=0, contrib;
    double *newintens = (double *) S_alloc ( m->nintens , sizeof(double) );
    for (i=0; i < d->nobs; ++i)
    {
	AddCovs(i, d, m, newintens);
	if ((m->death) && (d->tostate[i] == m->nst - 1) && (!m->exacttimes)) 
	{
	    contrib=0;
	    for (j = d->state[i]; j < d->tostate[i]; ++j)
		contrib += 
		    pijt(d->state[i], j, d->time[i] - 1 / d->tunit, newintens, m->qvector, m->nst, m->exacttimes)*
		    pijt(j, d->tostate[i], 1 / d->tunit, newintens, m->qvector, m->nst, m->exacttimes);
	    lik += log(contrib);
	}
	else 
	    lik += log(pijt(d->state[i], d->tostate[i], d->time[i], newintens, m->qvector, m->nst, m->exacttimes));
    }

    return (-2*lik); 
}

/* Calculates the most likely path through underlying states */ 

void Viterbi(data *d, model *m, double *fitted)
{
    int i, true, k, kmax, obs;
    double *newintens   = (double *) S_alloc(m->nintens, sizeof(double) );
    double *newmisc     = (double *) S_alloc(m->nmisc, sizeof(double));
    double *pmat = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *pmatdeath   = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *pmattodeath = (double *) S_alloc((m->nst)*(m->nst), sizeof(double)); 
    int *ptr = (int *) S_alloc((d->nobs)*(m->nst), sizeof(int));
    double *lvold = (double *) S_alloc(m->nst, sizeof(double));
    double *lvnew = (double *) S_alloc(m->nst, sizeof(double));
    double maxk, try, dt;

    for (k = 0; k < m->nst; ++k) 
	lvold[k] = 0;

    for (i = 1; i <= d->nobs; ++i)
    {
	if ((d->subject[i] == d->subject[i-1]) && (i < d->nobs))
	{
	    dt = d->time[i] - d->time[i-1];
	    AddCovs(i-1 + m->covmatch, d, m, newintens);
	    AddMiscCovs(i-1 + m->covmatch, d, m, newmisc);
	    Pmat(pmat, dt, newintens, m->qvector, m->nst, m->exacttimes);
	    if ((d->state[k] == m->nst - 1) && (m->death)){
		Pmat(pmattodeath, dt - 1 / d->tunit, newintens, m->qvector, m->nst, m->exacttimes);
		Pmat(pmatdeath, 1 / d->tunit, newintens, m->qvector, m->nst, m->exacttimes);
	    }
	    /* TODO: some sort of utility function for maxima and positional maxima ? */
	    for (true = 0; true < m->nst; ++true)
	    {
		kmax = 0;
		maxk = lvold[0] + log( pmat[MI(0, true, m->nst)] );
		if (true > 0) 
		{
		    for (k = 1; k < m->nst; ++k){
			try = lvold[k] + log( pmat[MI(k, true, m->nst)] );
			if (try > maxk) {
			    maxk = try;
			    kmax = k;
			}
		    }
		}
		lvnew[true] = log( PObsTrue(d->state[i], true, newmisc, m) )  +  maxk;
		ptr[MI(i, true, m->nst)] = kmax;
	    }
	    for (k = 0; k < m->nst; ++k)
		lvold[k] = lvnew[k];
	}
	else 
	{
	    /* Traceback for current patient */
	    maxk = lvold[0];
	    kmax = 0;
	    for (k=1; k < m->nst; ++k)
	    {
		try = lvold[k];
		if (try > maxk) {
		    maxk = try;
		    kmax = k;
		}
	    }
	  
	    obs = i-1;
	    fitted[obs] = kmax;
	    while   ( (d->subject[obs] == d->subject[obs-1]) && (obs > 0) ) 
	    {
		fitted[obs-1] = ptr[MI(obs, fitted[obs], m->nst)];
		--obs;
	    }
	    fitted[obs] = 0;

	    for (k = 0; k < m->nst; ++k) 
		lvold[k] = 0;
	}
    }
}
