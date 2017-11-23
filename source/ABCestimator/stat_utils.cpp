
#include <stdio.h>
//#include <ctype.h>  //for functions like isgraph, ....
#include <time.h>

#include "stat_utils.h"
#define _GCC_
#ifdef _GCC_
   #include <stdlib.h>
   #include <unistd.h>
   #include <sys/time.h>
#else
	#include "windows.h"    //BORLAND SPECIFIC
    #include   <dos.h>      //FOR STRUCT t in get_randomSeedFromCurrentTime()
#endif

//===============================================================================
//SECTION 1) RANDOM GENERATOR
//------------------------------------------------------------------------------
#define MBIG 1000000000L
#define MSEED 161803398L
#define MZ 0
#define FAC (1.0/MBIG)

#include <fstream>
#include <iostream>
using namespace std;
#define AOS cout



double ran3(long *idum)
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if ((*idum < 0) || (iff == 0) ) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;++i) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;++k)
			for (i=1;i<=55;++i) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;

	//for debug
	if( (mj*FAC) >= 1.0){
		AOS << " ** strange ran3:" << (mj*FAC) << endl;
	}
	//end debug

	return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
//------------------------------------------------------------------------------
//This procedure will initialize ran3 with a random value
//		if seed < 0 , then seed will be used to seed ran3
//    if seed >= 0 , a random seed less than 100,000 will be used to seed ran3
#ifdef _GCC_
int my_rand(){
	//randomize();
	//return rand();
	return 232;
}
#else
int my_rand(){
	randomize();
	return rand();
}
#endif



void
initialize_ran3(const long& seed) {
   if (seed<0) {
           if (seed<-161803398) { //if seed is too big to initialize ran3 properly...
      		AOS << "initialize_ran3(const int seed)" << endl;
        	 //Then choose a random seed;
      		long lidum=-(my_rand());
      		for( ; lidum < -161803398 && lidum>=0;){
      			lidum=-(my_rand());
      		}
      		ran3(&lidum);
                AOS << "Given seed was too big, so Random generator initialized with seed = " << lidum << endl;
         	return;
           }
      	long lidum=seed;
        ran3(&lidum);

      	AOS << "Random generator initialized with seed = " << seed << endl;

      	return;
   }
   else {
      long lidum=-(my_rand());
      for( ; lidum < -161803398 && lidum>=0;){
      	lidum=-(my_rand());
      }
      AOS << "\tRandom generator initialized with lidum = " << lidum << endl;
      ran3(&lidum);
   }
};
//------------------------------------------------------------------------------
//stef_22_12_98
//returns a long integer between 0 and max_val-1
//inline
long rand_long(const long& max_val) {
   long seed=1;
	long curval= (long) (ran3(&seed)*max_val);
   if (curval==max_val) return max_val-1;
   else return curval;
}

//------------------------------------------------------------------------------
//return a seed for the gcc version

long get_randomSeedFromCurrentTime(){
#ifdef _GCC_
   long seed;
   struct timeval tv ;
   struct timezone tz;
   gettimeofday(&tv, &tz);
   seed=(long)tv.tv_usec;
   return -seed;
#else
   //time_t start=clock();
   time_t start= time(NULL);
   start=start%10000000; //nico, 16.10.02, to get a smaller number

   struct  time t;
   gettime(&t);
   start+=t.ti_hund; //nico, 16.10.02, to add hundreds of seconds, in case several jobs are launched quickly

   return (long) -start;
#endif
};

//------------------------------------------------------------------------------
// random number of a certain distribution

/* returns a variate that is uniformly distributed on the interval [a,b]*/
double UniformRandom (const double& a, const double& b){
  long lidum=1;
  return (ran3(&lidum) * (b - a) + a);
} /* UniformRandom */

/* returns a variate that is log-uniformly distributed on the interval [a,b]*/
double LogUniformRandom (const double& a, const double& b){
  long lidum=1;
  if(!a) return (1e-6 * pow(b/1e-6, ran3(&lidum)));  // if a = 0;
  return (a * pow(b/a, ran3(&lidum)));
} /* LogUniformRandom */

/* Returns a Normal random variate based on a unit variate,
   using a random generator as a source of uniform deviates.
   Adapted from the algorithm described in the book Numerical Recipes by
   Press et al.
*/
double NormalRandom (const double& dMean, const double& dStdDev){
  double w, x1, x2, dTemp3;
  long lidum=1;

   do {
      x1 = 2. * ran3(&lidum) - 1.;
      x2 = 2. * ran3(&lidum) - 1.;
      w = x1*x1 + x2*x2;
   } while (w >= 1. || w < 1E-30);

   w = sqrt((-2.*log(w))/w);
   x1 *= w;
   return (x1 * dStdDev + dMean);
} /* NormalRandom */

/* returns a variate such that the log of the variate is normally distributed.*/
double LogNormalRandom (const double& dMean, const double& dStdDev){
   if(dMean<=0) AOS << "ERROR LogNormalRandom: log dMean" << endl;
   if(dStdDev<=0) AOS << "ERROR LogNormalRandom: log dStdDev" << endl;

   return exp (NormalRandom (log (dMean), log (dStdDev)));
} /* LogNormalRandom */

//==============================================================================
/* poisson distribution */
 float poidev(float xm,  long *idum){
	static float sq,alxm,g,oldm=(-1.0);
	float em,t,y;

   //nico: changed x with xm in order to compile ??
	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran3(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran3(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran3(idum) > t);
	}
	return em;
};
//----------------------------------------------------------------

float gammln(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
//------------------------------------------------------

#ifndef _GCC_
float gamdev(int ia){
	void nrerror(char error_text[]);
	int j;
	float am,e,s,v1,v2,x,y;
        long idum = 1;

	if (ia < 1) nrerror("Error in routine gamdev");
	if (ia < 6) {
		x=1.0;
		for (j=1;j<=ia;j++) x *= ran3(&idum);
		x = -log(x);
	} else {
		do {
			do {
				do {
					v1=2.0*ran3(&idum)-1.0;
					v2=2.0*ran3(&idum)-1.0;
				} while (v1*v1+v2*v2 > 1.0);
				y=v2/v1;
				am=ia-1;
				s=sqrt(2.0*am+1.0);
				x=s*y+am;
			} while (x <= 0.0);
			e=(1.0+y*y)*exp(am*log(x/am)-s*y);
		} while (ran3(&idum) > e);
	}
	return x;
}
#endif

//-------------------------------------------------------
/* binomial distribution */
float bnldev(float pp, int n){
	int j;
   long idum = 1;
	static int nold=(-1);
	float am,em,g,angle,p,bnl,sq,t,y;
	static float pold=(-1.0),pc,plog,pclog,en,oldg;

	p=(pp <= 0.5 ? pp : 1.0-pp);
	am=n*p;
	if (n < 25) {
		bnl=0;
		for (j=1;j<=n;j++)
			if (ran3(&idum) < p) ++bnl;
	} else if (am < 1.0) {
		g=exp(-am);
		t=1.0;
		for (j=0;j<=n;j++) {
			t *= ran3(&idum);
			if (t < g) break;
		}
		bnl=(j <= n ? j : n);
	} else {
		if (n != nold) {
			en=n;
			oldg=gammln(en+1.0);
			nold=n;
		} if (p != pold) {
			pc=1.0-p;
			plog=log(p);
			pclog=log(pc);
			pold=p;
		}
		sq=sqrt(2.0*am*pc);
		do {
			do {
				angle=PI*ran3(&idum);
				y=tan(angle);
				em=sq*y+am;
			} while (em < 0.0 || em >= (en+1.0));
			em=floor(em);
			t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
				-gammln(en-em+1.0)+em*plog+(en-em)*pclog);
		} while (ran3(&idum) > t);
		bnl=em;
	}
	if (p != pp) bnl=n-bnl;
	return bnl;
}
//--------------------------------------------------------
//------------------------------------------------------------------------------
//routine to initialise the probabilities for the Gamma distributed mutation rates.
// 1) we generate n_gen random values according to the gemma density function by using
//    the routine gamma_dev. Let the vector be v_gamma_dev
// 2) we sort the resulting values and make n_loci classes
// 3) Pi = v_gama_dev[ni] where ni is the center of the i'th class

void init_gamma_weights(double* p, int n_loci, double alph){
	int n_gen=1000000;
   double *v_gama_dev;
   try {
   	v_gama_dev=new double[n_gen];
   }
   catch (...) {
   	AOS << "Memory allocation error, exiting ...\n";
      exit(0);
   }
   double *next=&v_gama_dev[0];
   for(int i=0;i<n_gen;i++, ++next){
   	//v_gama_dev[i]=(float) gamma_dev(alph);
   	*next=(double) gamma_dev(alph);
   }
   //sort the elements
	quicksort(v_gama_dev, 0,n_gen-1);

   //compute the Pi:
   int num_in_cl=n_gen/n_loci;
   int count=0;
   for(int i= num_in_cl/2; i<n_gen && count<n_loci; i+=num_in_cl){
   	p[count]=v_gama_dev[i];
      count++;
   }
   delete[]  v_gama_dev;
}

//------------------------------------------------------------------------------
void init_gamma_weights(double* p, int n_loci, double alph, const int num_rates){
	int n_gen=n_loci;
   double *v_gama_dev, *working_space;
   try {
   	v_gama_dev=new double[n_gen];
      working_space= new double[4*num_rates];
   }
   catch (...) {
   	if (v_gama_dev) delete[] v_gama_dev;
   	if (working_space) delete[] working_space;
   	AOS << "Memory allocation error, exiting ...\n";
      exit(0);
   }

   Rates4Sites (v_gama_dev, alph, num_rates , n_loci, 0, working_space);

   //sort the elements
	quicksort(v_gama_dev, 0,n_gen-1);
	for(int i=0; i<n_gen ;++i){
   	p[i]=v_gama_dev[i];
   }
   if (v_gama_dev) delete[]  v_gama_dev;
   if (working_space) delete[] working_space;
}

//------------------------------------------------------------------------------
/** sorting an array */
template<typename T>
T partition(T* a, int left, int right, int pivotIndex) {
	int storeIndex=left;
	//swap
	T temp=a[right];
	a[right]=a[pivotIndex];
	a[pivotIndex]=temp;
	for (int i=left; i < right; i++){
		if(a[i] <= a[right]){
			//swap
			temp=a[storeIndex];
			a[storeIndex]=a[i];
			a[i]=temp;
			++storeIndex;
		}
	}
	//Pivot back to final place
	temp=a[right];
	a[right]=a[storeIndex];
	a[storeIndex]=temp;
	return storeIndex;
}
template<typename T>
void quicksort(T* a, int left, int right){
	if(right > left){
		int pivotIndex=(left+right)/2;
		pivotIndex=partition(a, left, right, pivotIndex);
		quicksort(a, left, pivotIndex-1);
		quicksort(a, pivotIndex+1, right);
	}
}

double partitionP(double** a, int left, int right, int pivotIndex) {
	int storeIndex=left;
	//swap
	double* temp=a[right];
	a[right]=a[pivotIndex];
	a[pivotIndex]=temp;
	for (int i=left; i < right; i++){
		if(*a[i] <= *a[right]){
			//swap
			temp=a[storeIndex];
			a[storeIndex]=a[i];
			a[i]=temp;
			++storeIndex;
		}
	}
	//Pivot back to final place
	temp=a[right];
	a[right]=a[storeIndex];
	a[storeIndex]=temp;
	return storeIndex;
}

void quicksortP(double** a, int left, int right){
	if(right > left){
		int pivotIndex=(left+right)/2;
		pivotIndex=partitionP(a, left, right, pivotIndex);
		quicksortP(a, left, pivotIndex-1);
		quicksortP(a, pivotIndex+1, right);
	}
}


 /*
//------------------------------------------------------------------------------
// sorting an array
template<typename T>
void quicksort(int left, int right, T* array){
	int l=left-1, r=right+1, m;
	T temp;
   m=(left+right)/2;
	do{
		while(array[++l]<array[m]);
      while(array[m]<array[--r]);
		if(l <= r) {
			// swap
         temp = array[l];
			array[l] = array[r];
			array[r] = temp;
			l++;
         r--;
		}
	}while (l<=r);
	if(left<r) quicksort(left, r, array);
	if(l<right) quicksort(l, right, array);
}
 */
//------------------------------------------------------------------------------
/*
Rogers comments:

In answer to your question about my algorithm, I'm going to append my
entire gamma_dev function.  The comments at the top provide references
to the original sources.  The algorithm I use is supposedly the most
commonly used when alpha<1.

In case it is relevant, let me tell you about some of the trouble I've
run into generating gamma deviates with small values of alpha.  My
first gamma_dev function was in single precision.  It behaved very
strangely.  When alpha<0.1, the number of segregating sites went *up*
as alpha went *down*, which makes no sense at all.  I couldn't find
any error in the code, but I noticed that the code does things that
may stretch the limits of floating point arithmetic.  So I recompiled
using double precision for all variables within gamma_dev.  The
strange behavior went away.

The literature doesn't say much about the stability of these
algorithms when alpha is very small.  It seems that no one has ever
been interested in that case.  I'll bet that none of the commercial
statistical packages have tested their gamma deviate generator with
very small alpha values either.  Consequently, we can't test our
algorithms by comparing the quantiles of our generated values with
those generated by, say, SPSS.  The only sure way is to calculate
quantiles by direct integration of the density function.  I have done
this for alpha=0.1 and am about to compare the quantiles of my numbers
with these values.  I'll let you know what happens.

Alan

PS  Here's the code along with references.  */

/****************************************************************
Random deviates from standard gamma distribution with density
         a-1
        x    exp[ -x ]
f(x) = ----------------
         Gamma[a]

where a is the shape parameter.  The algorithm for integer a comes
from numerical recipes, 2nd edition, pp 292-293.  The algorithm for
a<1 uses code from p 213 of Statistical Computing, by Kennedy and
Gentle, 1980 edition.  This algorithm was originally published in:

Ahrens, J.H. and U. Dieter (1974), "Computer methods for sampling from
Gamma, Beta, Poisson, and Binomial Distributions".  COMPUTING
12:223-246.

The mean and variance of these values are both supposed to equal a.
My tests indicate that they do.

This algorithm has problems when a is small.  In single precision, the
problem  arises when a<0.1, roughly.  That is why I have declared
everything as double below.  Trouble is, I still don't know how small
a can be without causing trouble.  Mean and variance are ok at least
down to a=0.01.  f(x) doesn't seem to have a series expansion around
x=0.
****************************************************************/
double
gamma_dev(double a, double b) {
   if(b <= 0) {
      AOS << "\ngamma_dev: parameter must be positive\n";
      exit(1);
   }
   return gamma_dev(a)/b;
}


double
gamma_dev(double a) {

  int ia;
  double u, b, p, x, y=0.0, recip_a;
  long lidum=1L;

  if(a <= 0) {
    AOS << "\ngamma_dev: parameter must be positive\n";
    exit(1);
  }

  ia = (int) floor(a);  /* integer part */
  a -= ia;        /* fractional part */
  if(ia > 0) {
    y = igamma_dev(ia);  /* gamma deviate w/ integer argument ia */
    if(a==0.0) return(y);
  }

  /* get gamma deviate with fractional argument "a" */
  b = (M_E + a)/M_E;
  recip_a = 1.0/a;
  for(;;) {
    u = ran3(&lidum);
    p = b*u;
    if(p > 1) {
      x = -log((b-p)/a);
      if( ran3(&lidum) > pow(x, a-1)) continue;
      break;
    }
    else {
      x = pow(p, recip_a);
      if( ran3(&lidum) > exp(-x)) continue;
      break;
    }
  }
  return(x+y);
}
/****************************************************************
gamma deviate for integer shape argument.  Code modified from pp
292-293 of Numerical Recipes in C, 2nd edition.
****************************************************************/
double igamma_dev(int ia){
  int j;
  double am,e,s,v1,v2,x,y;
  long lidum=1L;

  if (ia < 1)
  {
    AOS << "\nError: arg of igamma_dev was <1\n";
    exit(1);
  }
  if (ia < 6)
  {
    x=1.0;
    for (j=0; j<ia; j++)
      x *= ran3(&lidum);
    x = -log(x);
  }else
  {
    do
    {
      do
      {
	do
	{                         /* next 4 lines are equivalent */
	  v1=2.0*ran3(&lidum)-1.0;       /* to y = tan(Pi * uni()).     */
	  v2=2.0*ran3(&lidum)-1.0;
	}while (v1*v1+v2*v2 > 1.0);
	y=v2/v1;
	am=ia-1;
	s=sqrt(2.0*am+1.0);
	x=s*y+am;
      }while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    }while (ran3(&lidum) > e);
  }
  return(x);
}

//------------------------------------------------------------------------------

#define FOR(i,n) for(i=0; i<n; i++)
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

//------------------------------------------------------------------------------
int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf,
    double space[])
{
/* Rates for sites from the gamma (ncatG=0) or discrete-gamma (ncatG>1).
   Rates are converted into the c.d.f. if cdf=1, which is useful for
   simulation under JC69-like models.

	. 	rates[] is a vector of size ls, and will have, in return the gamma
  		variates for sites.
	. 	ls means "length of sequence" in number of sites.
	. 	alpha is the parameter.  I fix the other parameter beta = alpha so that
  		the mean of the distribution is one.
	. 	ncatG is the number of categories in the discrete distribution.
  		ncatG = 0 means the continuous gamma distribution.  Use 8.
	. 	space[] is a vector of working space, with a minimum of size 3*ncatG
	. 	cdf stands for cumulative distribution function.

*/
   int h, ir,j, *counts=(int*)(space+2*ncatG);
   double *rK=space, *freqK=space+ncatG;

   if (alpha==0)
      for (h=0; h<ls; h++) rates[h]=1;
   else {
      if (ncatG>1) {
         DiscreteGamma (freqK, rK, alpha, alpha, ncatG, 0);
         MultiNomial (ls, ncatG, freqK, counts, space+3*ncatG);
         for (ir=0,h=0; ir<ncatG; ir++)
            for (j=0; j<counts[ir]; j++)  rates[h++]=rK[ir];
      }
      else
         for (h=0; h<ls; h++) rates[h]=rndgamma(alpha)/alpha;
      if (cdf) {
         for (h=1; h<ls; h++) rates[h]+=rates[h-1];
         abyx (1/rates[ls-1], rates, ls);
      }
   }
   return (0);
}
//------------------------------------------------------------------------------
double LnGamma (double alpha){
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673
          + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
               +.083333333333333)/x;
}
//------------------------------------------------------------------------------
double DFGamma(double x, double alpha, double beta){
   if (alpha<=0 || beta<=0) AOS << "err in DFGamma()\n";
   if (alpha>100) AOS << "large alpha in DFGamma()\n";
   return pow(beta*x,alpha)/x * exp(-beta*x - LnGamma(alpha));

}
//------------------------------------------------------------------------------
#pragma warn -8004
double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);
   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (gin);
}

//------------------------------------------------------------------------------
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
double PointChi2 (double prob, double v){
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      AOS << "\nerr IncompleteGamma\n";
      return (-1);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}
//------------------------------------------------------------------------------
/* discretization of gamma distribution with equal proportions in each
   category
*/
int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median){
   int i;
   double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if (median) {
      for (i=0; i<K; i++) rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
         freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
         freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}

//----------------------------------------------------------------------------
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
       points of the normal distribution.  26: 118-121.

*/
double PointNormal (double prob){
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);
   y = sqrt (log(1/(p1*p1)));
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}
//------------------------------------------------------------------------------
double rndgamma1 (double s);
double rndgamma2 (double s);

/* random standard gamma (Mean=Var=s,  with shape par=s, scale par=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1
*/
double rndgamma (double s){
   double r=0;
   if (s<=0)      AOS << "jgl gamma..\n";
   else if (s<1)  r=rndgamma1 (s);
   else if (s>1)  r=rndgamma2 (s);
   else           r=-log(rndu());
   return (r);
}

//------------------------------------------------------------------------------
/* random standard gamma for s<1
   switching method
*/
double rndgamma1 (double s){
   double r, x=0,small=1e-37,w;
   static double a,p,uf,ss=10,d;

   if (s!=ss) {
      a=1-s;
      p=a/(a+s*exp(-a));
      uf=p*pow(small/a,s);
      d=a*log(a);
      ss=s;
   }
   for (;;) {
      r=rndu();
		if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
      else if (r>uf)  x=a*pow(r/p,1/s), w=x;
      else            return (0);
      r=rndu ();
      if (1-r<=w && r>0)
         if (r*(w+1)>=1 || -log(r)<=w)  continue;
      break;
   }
   return (x);
}

//------------------------------------------------------------------------------
/* random standard gamma for s>1
   Best's (1978) t distribution method
*/
double rndgamma2 (double s){
   double r,d,f,g,x;
   static double b,h,ss=0;
   if (s!=ss) {
      b=s-1;
      h=sqrt(3*s-0.75);
      ss=s;
   }
   for (;;) {
      r=rndu ();
      g=r-r*r;
      f=(r-0.5)*h/sqrt(g);
      x=b+f;
      if (x <= 0) continue;
      r=rndu();
      d=64*r*r*g*g*g;
      if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
   }
   return (x);
}

//------------------------------------------------------------------------------
int xtoy (double x[], double y[], int n){
   int i;
   for (i=0; i<n; y[i]=x[i],i++);
   return(0);
}

//------------------------------------------------------------------------------
int abyx (double a, double x[], int n){
   int i;
   for (i=0; i<n; x[i]*=a,i++) ;
   return(0);
}

//------------------------------------------------------------------------------
static int z_rndu=137;
static unsigned w_rndu=13757;
//------------------------------------------------------------------------------
void SetSeed (int seed){
   z_rndu = 170*(seed%178) + 137;
   w_rndu=seed;
}

//------------------------------------------------------------------------------
/* U(0,1): AS 183: Appl. Stat. 31:188-190
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
*/
double rndu (void){
   static int x_rndu=11, y_rndu=23;
   double r;

   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
   r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
   return (r-(int)r);
}

//------------------------------------------------------------------------------
/* sample n times from a mutinomial distribution M(ncat, prob[])
   prob[] is considered cumulative prob if (space==NULL)
*/
int MultiNomial (int n, int ncat, double prob[], int nobs[], double space[]){
   int i, j, crude=(ncat>20), ncrude=5, lcrude[5];
   double r, *pcdf=(space==NULL?prob:space);

   FOR (i, ncat) nobs[i]=0;
   if (space) {
      xtoy (prob, pcdf, ncat);
      for (i=1; i<ncat; i++) pcdf[i]+=pcdf[i-1];
   }
   if (fabs(pcdf[ncat-1]-1) > 1e-5) AOS << "sum P!=1 in MultiNomial\n";
   if (crude) {
      for (j=1,lcrude[0]=i=0; j<ncrude; j++)  {
         while (pcdf[i]<(double)j/ncrude) i++;
         lcrude[j]=i-1;
      }
   }
   FOR (i, n) {
      r=rndu();
      j=0;
      if (crude) {
         for (; j<ncrude; j++) if (r<(j+1.)/ncrude) break;
         j=lcrude[j];
      }
      for (; j<ncat; j++) if (r<pcdf[j]) break;
      nobs[j] ++;
   }
   return (0);
}

/* -----------------------------------------------------------------------------
   BetaRandom

   returns a variate that is Beta distributed on the interval [a,b]
   with shape parameters alpha and beta.

   The Beta function has two shaping parameters, alpha and beta.
   Setting these parameters to 1.5 and 1.5 yields a normal-like
   distribution, but without tails. If alpha and beta are equal to
   1 it is a uniform distribution.

   If alpha and beta are less than 1, use a rejection algorithm;
   Otherwise use the fact that if x is distributed Gamma(alpha) and y
   Gamma(beta) then x/(x+y) is Beta(alpha, beta).

   The rejection algorithm first a Beta variate is found over the
   interval [0, 1] with not the most efficient algorithm.  This is then
   scaled at the end to desired range.

   It may be tempting to re-use the second number drawn as the first
   random number of the next iteration, and simply draw one more.
   *** Don't do it.  You will produce an incorrect distribution.  You
   must draw two new numbers for the rejection sampling to be correct.

   References:
   - Ripley, Stochastic Simulations, John Wiley and Sons, 1987, p 90.
   - J.H.Maindonald, Statistical Computation, John Wiley and Sons,
     1984, p 370.
*/
double BetaRandom (double alpha, double beta, double a, double b){
  if (b <= a) {
    AOS << "Error: bad shape or range for a beta variate - Exiting\n\n";
    exit (1);
  }
  return (a + BetaRandom(alpha, beta) * (b-a));   /* Scale to interval [a, b] */
}


double BetaRandom (double alpha, double beta) {
  double u1, u2, w;

  if (alpha <= 0 || beta <= 0) {
    AOS << "Error: bad shape or range for a beta variate - Exiting\n\n";
    exit (1);
  }

  if ((alpha < 1) && (beta < 1))
    /* use rejection */
    do {
      long lidum=1;
      u1 = ran3(&lidum); /* Draw two numbers */
      u2 = ran3(&lidum);

      u1 = pow(u1, 1/alpha); /* alpha and beta are > 0 */
      u2 = pow(u2, 1/beta);

      w = u1 + u2;

    } while (w > 1.0);

  else {
    /* use relation to Gamma */
    u1 = gamma_dev(alpha);
    u2 = gamma_dev(beta);
    w  = u1 + u2;
  }

  return u1/w;

} /* BetaRandom */

//------------------------------------------------------------------------------
//Returns a geometrically distributed random variate between zero and infinity
int geometric(const double &p) {
   if (p < 1.e-10 || p==1.0) return 0;
   long lidum=1L;
   return (int) floor(log(ran3(&lidum))/log(p));
}

//------------------------------------------------------------------------------
//Epanechnikov Kernel Ph�ntu 2006 (modfied from Laurent)
//"points" is a pointer to an array (of size numPoints) with normalized points
//used for the density estimation
//"position" is the normalized position for which the estimate is returned
double epanechnikovDensityEstimate(double* points, int numPoints, double position, double bandwidth, double c){
	double curDensity=0.0, delta;
	for(int i=0; i<numPoints; ++i){
		delta=fabs(position-points[i]);
		if (delta){
			delta=delta/bandwidth;
			if(fabs(delta)<=1.0) curDensity+=c*(1.0-delta*delta)/bandwidth;
		}
	}
	return curDensity;
}


//------------------------------------------------------------------------------            */
#pragma warn +8004
