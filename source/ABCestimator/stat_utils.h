#ifndef __STAT_UTIL_HPP__
#define __STAT_UTIL_HPP__

#include <math.h>

// A routine for generating good random numbers
double ran3(long *idum);
long rand_long(const long& max_val);    //random integer betwee 0 and max_val-1
void initialize_ran3(const long& seed); //initializes the random generator
float poidev(float xm, long *idum);     // poisson distribution
float bnldev(float pp, int n);          // binomial distribution
double UniformRandom (const double& a, const double& b);
double LogUniformRandom (const double& a, const double& b);
double NormalRandom (const double& dMean, const double& dStdDev);
double LogNormalRandom (const double& dMean, const double& dStdDev);

#ifndef _GCC_
float gamdev(int ia, long *idum);
#endif

//return a seed for the random number generator  (GCC version)
long get_randomSeedFromCurrentTime();

float gammln(float xx);

double partitionP(double** a, int left, int right, int pivotIndex);
void quicksortP(double** a, int left, int right);

template<typename T>
T partition(T* a, int left, int right, int pivotIndex);
template<typename T>
extern void quicksort(T* a, int left, int right);

//Density estimates
double epanechnikovDensityEstimate(double* points, int numPoints, double position, double bandwidth=0.4, double c=1.5);

// A routine to initialize the probabilities for the Gamma distributed mutation rates.
void init_gamma_weights(double* p, int n_loci, double alph);
//Discrete gamma case, using Ziheng Yang's routines
void init_gamma_weights(double* p, int n_loci, double alph, const int num_rates);
//returns gamma deviates
extern double gamma_dev(double a);
extern double gamma_dev(double a, double b);
extern double igamma_dev(int ia);
extern double BetaRandom (double alpha, double beta);
extern double BetaRandom (double alpha, double beta, double a, double b);
extern int geometric(const double &p);

//Ziheng Yang routines for generating discrete gamma distributions
//------------------------------------------------------------------------------
extern int Rates4Sites (double rates[],double alpha,int ncatG,int ls, int cdf, double space[]);
extern double LnGamma (double alpha);
extern double DFGamma(double x, double alpha, double beta);
extern double IncompleteGamma (double x, double alpha, double ln_gamma_alpha);
extern double PointChi2 (double prob, double v);
extern double PointNormal (double prob);
extern int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median);
extern double rndgamma (double s);
extern int MultiNomial (int n, int ncat, double prob[], int nobs[], double space[]);
extern double rndu (void);
extern int xtoy (double x[], double y[], int n);
extern int abyx (double a, double x[], int n);
//------------------------------------------------------------------------------
#define PI 3.141592654


#endif










