/* *****************************************************************
   PROGRAM: lik.c 
   AUTHOR:  Chris Jackson
   DATE:    November 2001

   Routines for calculating likelihoods for Markov models with and without misclassification

   ******************************************************************  */ 


#include "lik.h"
#define NODEBUG
#define NOLIKDEBUG

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
    int *baseconstraint, /* constraints on baseline transition intensities */
    int *basemiscconstraint, /* constraints on baseline misclassification probabilities */
    double *initprobs, /* initial state occupancy probabilities */
    int *nst,      /* number of Markov states */
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

    fillparvec ( intens,        *p,  params, allinits, *nfix, fixedpars, *nintenseffs,  &ifix, &iopt, &iall ) ;
    fillparvec ( coveffect,     *p,  params, allinits, *nfix, fixedpars, *ncoveffs,     &ifix, &iopt, &iall ) ;
    fillparvec ( miscprobs,     *p,  params, allinits, *nfix, fixedpars, *nmisceffs,    &ifix, &iopt, &iall ) ;
    fillparvec ( misccoveffect, *p,  params, allinits, *nfix, fixedpars, *nmisccoveffs, &ifix, &iopt, &iall ) ;

    d.subject = subjvec;     d.time = timevec;     d.state = statevec;     d.tostate = tostatevec;
    d.cov = covvec;     d.misccov = misccovvec;    d.fromto = *fromto;   d.nobs = *nobs;   
    d.npts = *npts;    d.ncovs = *ncovs;    d.nmisccovs = *nmisccovs;     

    m.qvector = qvector;    m.evector = evector;    m.constraint = constraint;    m.miscconstraint = miscconstraint;
    m.baseconstraint = baseconstraint;    m.basemiscconstraint = basemiscconstraint;
    m.nst = *nst;  m.nms = *nms;  m.nintens = *nintens;  m.nmisc = *nmisc;  m.ncoveffs = *ncoveffs;
    m.nintenseffs = *nintenseffs; m.nmisceffs = *nmisceffs; 
    m.nmisccoveffs = *nmisccoveffs;  m.covmatch = *covmatch; m.ndeath=*ndeath; m.death = death;  
    m.ncens = *ncens; m.censor = censor; m.censstates=censstates; m.censstind=censstind;
    m.exacttimes = *exacttimes;
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
#ifdef DEBUG
	    printf("pt %d, lik %lf\n", pt, likone);
#endif
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
		int p, /* number of parameters optimised over */
		double *params, /* current values of parameters being optimised over */
		double *allinits, /* full vector of initial values */
		int nfix,   /* number of fixed parameters, ie length of fixedpars */
		int *fixedpars,  /* indices of allinits which are fixed */
		int ni,    /* length of parvec */
		int *ifix,  /* current index into fixedpars */
		int *iopt,  /* current index into params */
		int *iall   /* current index into allinits */
    )
{
    int i;
    for (i=0; i<ni; ++i, ++(*iall)){
	if ((*ifix < nfix) && (*iall == fixedpars[*ifix])) {
	    parvec[i] = allinits[*iall];
	    ++(*ifix);
	}
	else if (*iopt < p) {
	    parvec[i] = params[*iopt];
	    ++(*iopt);
	}
    }
}

/* Form vector of probabilities of a given outcome conditionally on each underlying state */
/* If outcome has nc=1 (not censored) then outcome probs are defined by misclassification matrix */
/* If outcome has nc>1 (censored) then no misclassification  */

void GetCensoredPObsTrue(double *pout, int *current, int nc, double *newmisc, model *m)
{
    int i, j;
    for (i=0; i<m->nst; ++i) {
	pout[i] = 0;
	if (nc == 1) 
	    pout[i] =  PObsTrue(current[0], i, newmisc, m);
	else { 
	    for (j=0; j<nc; ++j) 
		if (current[j] == i) 
		    pout[i] = 1;
	}
    }
}

/* Likelihood for the misclassification model for one individual */

double likmisc(int pt, /* ordinal subject ID */
	       data *d, model *m
    )
{
    int *current = (int *) S_alloc (1, sizeof(int));
    double *pout = (double *) S_alloc(m->nst, sizeof(double)); 
    double *cumprod     = (double *) S_alloc(m->nst, sizeof(double)); 
    double *newprod     = (double *) S_alloc(m->nst, sizeof(double));  
    double *newmisc     = (double *) S_alloc(m->nmisc, sizeof(double));
    double lweight, lik;
    int i, pti, k, first=0, last=0, nc=1;
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
    GetCensored(d->state[first], m, &nc, &current);
    GetCensoredPObsTrue(pout, current, nc, newmisc, m);
    for (i = 0; i < m->nst; ++i)
	cumprod[i] = pout[i] * m->initprobs[i]; /* cumulative matrix product */
    lweight=0;
    /* Matrix product loop to accumulate the likelihood */
    for (k = first+1; k <= last; ++k)
    {
	GetCensored(d->state[k], m, &nc, &current);
#ifdef DEBUG
	if (nc > 0) {
	    for (i=0; i<nc; ++i)
		printf("curr %d = %d, ", i, current[i]);
	    printf("\n");
	}
#endif
	UpdateLik(current, nc, d->time[k] - d->time[k-1],
		  k, last, d, m, cumprod, newprod, &lweight);
	for (i = 0; i < m->nst; ++i)
	    cumprod[i] = newprod[i];
#ifdef LIKDEBUG
	    for (i = 0; i < m->nst; ++i) {
		printf("cump %d = %lf, ", i, cumprod[i]);
	    }
	    printf("lweight = %lf\n", lweight);
#endif
    }
  
    lik=0;
    for (i = 0; i < m->nst; ++i)
	lik = lik + cumprod[i];

    /* Transform the likelihood back to the proper scale */
    return -2*(log(lik) - lweight); 
}


/* Post-multiply the row-vector cumprod by matrix T to accumulate the likelihood */

void UpdateLik(int *current, int nc, double dt, int k, int last, data *d, model *m, 
	       double *cumprod, double *newprod, double *lweight)
{
    int i, j;
    double newprod_ave;
    double *T           = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *newintens   = (double *) S_alloc(m->nintens, sizeof(double));
    double *newmisc     = (double *) S_alloc(m->nmisc, sizeof(double)); 
    double *pmat        = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    double *pout = (double *) S_alloc(m->nst, sizeof(double)); 
    AddCovs(k - 1 + m->covmatch, d, m, newintens);
    AddMiscCovs(k - 1 + m->covmatch, d, m, newmisc);
    GetCensoredPObsTrue(pout, current, nc, newmisc, m);
    
    /* calculate the transition probability (P) matrix for the time interval dt */
    Pmat(pmat, dt, newintens, m->qvector, m->nst, m->exacttimes);
    for(j = 0; j < m->nst; ++j)
    {
	newprod[j] = 0.0;
	for(i = 0; i < m->nst; ++i)
	{
	    if ((k == last) && (m->ndeath > 0) && is_element(current[0], m->death, m->ndeath) && !m->exacttimes) {
		/* last observation was death and death time known exactly */
		T[MI(i,j,m->nst)] = pmat[MI(i,j,m->nst)] * qij(j, current[0], newintens, m->qvector, m->nst);
#ifdef LIKDEBUG	       
		printf("obs %d, death i=%d, j=%d, state=%d, pmat=%lf, HM=%lf, HM=%lf, res=%lf\n", k, i, j, current[0], pmat[MI(i, j, m->nst)],
		       pout[j], 
		       PObsTrue(current[0], j, newmisc, m),
		       T[MI(i,j,m->nst)]);
#endif
	    }
	    else {
		T[MI(i,j,m->nst)] = pmat[MI(i, j, m->nst)] * pout[j]; 
#ifdef LIKDEBUG	       
		printf("obs %d, i=%d, j=%d, state=%d, pmat=%lf, HM=%lf, HM=%lf, res=%lf\n", k, i, j, current[0], pmat[MI(i, j, m->nst)],
		       pout[j], 
		       PObsTrue(current[0], j, newmisc, m),
		       T[MI(i,j,m->nst)]);
#endif
	    }
	    if (T[MI(i,j,m->nst)] < 0) T[MI(i,j,m->nst)] = 0;
	    newprod[j] = newprod[j] + cumprod[i]*T[MI(i,j,m->nst)];
	}
    }
    /* re-scale the likelihood at each step to prevent it getting too small and underflowing */
    /*  while cumulatively recording the log scale factor   */
#ifdef LIKDEBUG
    for (i = 0; i < m->nst; ++i) {
	printf("newp %d = %lf, ", i, newprod[i]);
    }
    printf("\n");
#endif
    for(i = 0, newprod_ave=0.0; i < m->nst ; ++i)
	newprod_ave += newprod[i];  
    newprod_ave /=  m->nst;
    if (newprod_ave == 0)
	newprod_ave = 1;
    for(i = 0; i < m->nst ; ++i)
	newprod[i] = newprod[i] / newprod_ave;
    *lweight -= log(newprod_ave);
}



/* Add effects of covariates to the transition rates at an observation number obs */

void AddCovs(int obs, data *d, model *m, double *newintens)
{
    double *lintens =  (double *) S_alloc ( m->nintens, sizeof(double));
    int j, k;
    for (j = 0; j < m->nintens; ++j)
    { 
	lintens[j] = log(m->intens[m->baseconstraint[j] - 1]);
	for (k = 0; k < d->ncovs; ++k)
	    lintens[j] += m->coveffect[m->constraint[MI(k, j, m->nintens)] - 1] * d->cov[MI(d->nobs, obs, k)];
	newintens[j] = exp(lintens[j]);
    }
}


/* Add effects of misclassification covariates to the misclassification rates at an observation number obs */

void AddMiscCovs(int obs, data *d, model *m, double *newp)
{
    double *logitp = (double *) S_alloc ( m->nmisc , sizeof(double));
    int j, k;
    for (j = 0; j < m->nmisc; ++j)
    { 
	logitp[j] = logit(m->miscprobs[m->basemiscconstraint[j] - 1]);
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
    FillQmatrix(m->evector, miscprobs, emat, m->nst);
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

/* Return a vector of the possible states that a censored state could represent */ 
/* These will be summed over when calculating the likelihood */

void GetCensored (int obs, model *m, int *nc, int **states) 
{
    int j, k=0, n, cens=0;
    if (m->ncens == 0)
	n = 1;
    else {
	while (obs != m->censor[k] && k < m->ncens)
	    ++k;
	if (k < m->ncens) { 
	    cens = 1; 
	    n =  m->censstind[k+1] - m->censstind[k];
	}
	else n = 1;
    } 
    *states = (int *) S_realloc ((char *) *states, n, *nc, sizeof(int));
    if (m->ncens == 0 || !cens)
	(*states)[0] = obs;
    else for (j = m->censstind[k]; j < m->censstind[k+1]; ++j)
	(*states)[j - m->censstind[k]] = m->censstates[j];
    *nc = n;    
}

/* Likelihood for the simple model. Data of form "subject ID, obs time, obs state" */

double liksimple(data *d, model *m)
{
    int i,j,k, nc=1, np=1;
    int *previous = (int *) S_alloc (1, sizeof(int));
    int *current = (int *) S_alloc (1, sizeof(int));
    double dt, lik=0, contrib, *newintens = (double *) S_alloc ( m->nintens , sizeof(double));
    for (i = 1; i < d->nobs; ++i){
	AddCovs(i - 1 + m->covmatch, d, m, newintens);
	if (d->subject[i-1] == d->subject[i]){ 
	    dt = d->time[i] - d->time[i-1];
	    GetCensored(d->state[i-1], m, &np, &previous);
	    GetCensored(d->state[i], m, &nc, &current);
	    if ((m->ndeath > 0) && is_element(d->state[i], m->death, m->ndeath) && !m->exacttimes)
		/* if state is a "death" state, i.e. entry date is known exactly, with unknown different state at the previous instant. */
		{
		    if (d->state[i-1] == d->state[i]) 
			contrib = 1; /* death-death transition has probability 1 */
		    else { 
			/*  Sum over states that state[i-1] can represent */
			contrib=0;
			for (j = 0; j < np; ++j)
			    for (k = 0; k < m->nst; ++k)
				if (k != d->state[i])
				    contrib += 
					pijt(previous[j], k,  dt, newintens, m->qvector, m->nst, 0)*
					qij(k, d->state[i], newintens, m->qvector, m->nst);
		    }
		    lik += log(contrib);
		}
	    else {
		/* Sum over states that state[i-1] can represent, and state[i] can represent */
		contrib = 0;
		for (j = 0; j < np; ++j) 
		    for (k = 0; k < nc; ++k) 
			contrib += pijt(previous[j], current[k], dt, newintens, m->qvector, m->nst, m->exacttimes);
		lik += log(contrib);
	    }
#ifdef DEBUG
 	    printf("prev=%d, curr=%d, timelag=%lf, con=%lf, lik = %lf\n", previous[0], current[0], dt, log(contrib), lik);
#endif
	}
    }
    return (-2*lik); 
}

/* Likelihood for the simple model. Data of form "from-state, to-state, time-difference" */

double liksimple_fromto(data *d, model *m)
{
    int i,j,k,np=1,nc=1;
    double lik=0, contrib;
    int *previous = (int *) S_alloc (1, sizeof(int));
    int *current = (int *) S_alloc (1, sizeof(int));
    double *newintens = (double *) S_alloc ( m->nintens , sizeof(double));
    for (i=0; i < d->nobs; ++i)
    {
	AddCovs(i, d, m, newintens);
	GetCensored(d->state[i], m, &np, &previous);
	GetCensored(d->tostate[i], m, &nc, &current);
	if ((m->ndeath > 0) && is_element(d->tostate[i], m->death, m->ndeath) && !m->exacttimes)
	{
	    if (d->state[i] == d->tostate[i]) 
		contrib = 1; /* death-death transition has probability 1 */
	    else { 
		contrib = 0;
		for (j = 0; j < np; ++j)
		    for (k = 0; k < m->nst; ++k)
			if (k != d->tostate[i])
			    contrib += 
				pijt(previous[j], k, d->time[i], newintens, m->qvector, m->nst, 0)*
				qij(k, d->tostate[i], newintens, m->qvector, m->nst);	
	    }
	    lik += log(contrib);
	}
	else {
	    contrib = 0;
	    for (j = 0; j < np; ++j) 
		for (k = 0; k < nc; ++k)
		    contrib += pijt(previous[j], current[k], d->time[i], newintens, m->qvector, m->nst, m->exacttimes);
	    lik += log(contrib);
	}
#ifdef DEBUG
 	printf("lik = %lf\n", lik);  
#endif
    }
    return (-2*lik); 
}

/* Calculates the most likely path through underlying states */ 

void Viterbi(data *d, model *m, double *fitted)
{
    int i, tru, k, kmax, obs, nc=1;
    double *newintens   = (double *) S_alloc(m->nintens, sizeof(double));
    double *newmisc     = (double *) S_alloc(m->nmisc, sizeof(double));
    double *pmat = (double *) S_alloc((m->nst)*(m->nst), sizeof(double));
    int *ptr = (int *) S_alloc((d->nobs)*(m->nst), sizeof(int));
    double *lvold = (double *) S_alloc(m->nst, sizeof(double));
    double *lvnew = (double *) S_alloc(m->nst, sizeof(double));
    int *current = (int *) S_alloc (1, sizeof(int));
    double *pout = (double *) S_alloc(m->nst, sizeof(double)); 
    double maxk, try, dt;
#ifdef DEBUG
    printf("Starting Viterbi algorithm...\n");
#endif

    for (k = 0; k < m->nst; ++k) 
	lvold[k] = 0;

    for (i = 1; i <= d->nobs; ++i)
    {
#ifdef DEBUG
	    printf("obs %d\n", i);
#endif
	if ( (i < d->nobs) && (d->subject[i] == d->subject[i-1]) )
	{
#ifdef DEBUG
	    printf("subject %d\n ", d->subject[i]);
#endif
	    dt = d->time[i] - d->time[i-1];
	    AddCovs(i-1 + m->covmatch, d, m, newintens);
	    AddMiscCovs(i-1 + m->covmatch, d, m, newmisc);
	    GetCensored(d->state[i], m, &nc, &current);
	    GetCensoredPObsTrue(pout, current, nc, newmisc, m);

	    Pmat(pmat, dt, newintens, m->qvector, m->nst, m->exacttimes);
	    /* TODO: some sort of utility function for maxima and positional maxima ? */
	    for (tru = 0; tru < m->nst; ++tru)
	    {
		kmax = 0;
		maxk = lvold[0] + log( pmat[MI(0, tru, m->nst)] );
		if (tru > 0) 
		{
		    for (k = 1; k < m->nst; ++k){
			try = lvold[k] + log( pmat[MI(k, tru, m->nst)] );
			if (try > maxk) {
			    maxk = try;
			    kmax = k;
			}
		    }
		}
		lvnew[tru] = log( pout[tru] )  +  maxk;
#ifdef DEBUG			   
		printf("obs %d, true %d, pout[%d] = %lf, lvnew = %lf, kmax=%d\n", d->state[i], tru, tru, pout[tru], lvnew[tru], kmax);
#endif
		ptr[MI(i, tru, m->nst)] = kmax;
	    }
	    for (k = 0; k < m->nst; ++k)
		lvold[k] = lvnew[k];
	}
	else 
	{
#ifdef DEBUG
	    printf("traceback for subject %d\n ", d->subject[i]);
#endif
	    /* Traceback for current patient */
	    maxk = lvold[0];
	    kmax = 0;
	    for (k=1; k < m->nst; ++k)
	    {
#ifdef DEBUG
		printf("lvold[%d] = %lf, ", k, lvold[k]);
#endif
		try = lvold[k];
		if (try > maxk) {
		    maxk = try;
		    kmax = k;
		}
	    }
#ifdef DEBUG
	    printf("kmax = %d\n", kmax);
#endif

	    obs = i-1;
	    fitted[obs] = kmax;
	    while   ( (obs > 0) && (d->subject[obs] == d->subject[obs-1]) ) 
	    {
		fitted[obs-1] = ptr[MI(obs, fitted[obs], m->nst)];
#ifdef DEBUG
		printf("ptr %d = %lf, ", obs-1, fitted[obs-1]);
#endif
		--obs;
	    }
	    fitted[obs] = 0;

	    for (k = 0; k < m->nst; ++k) 
		lvold[k] = 0;
	}
    }
}
