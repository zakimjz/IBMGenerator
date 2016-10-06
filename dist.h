#include "glob.h"

#define INIT_SEED -1

float ran0(long &idum);

//=============================  Distributions  =============================

class UniformDist;

class RandSeed
{
private:
  static UniformDist *ran1;
public:
  static void set_seed(long new_seed);
  static long new_seed(void);
	// Returns a random seed between 0 and MAXLONG, using
	// INIT_SEED to initialize the random sequence.
};


// Returns a random deviate between 0.0 and 1.0 (exclusive of
// the endpoint values). Call with a negative integer to
// initialize.
//
#define NTAB 32
class UniformDist
{
private:
  long idum;
  long iy;
  long iv[NTAB];
  
  float ran1(void);
public:
  UniformDist(void)
    : iy(0) { idum = RandSeed::new_seed(); };
  UniformDist(long seed)
    : iy(0), idum(seed) {};
  FLOAT operator()(void)	// returns a random number between 0 and 1
    { return ran1(); };
};


class PoissonDist
{
private:
  FLOAT lambda;
  float sq,alxm,g,oldm;
  UniformDist *ran1;

  float poidev(float xm);
  	// Returns as a floating-point number an integer value that is
	// a random deviate drawn from a Poisson distribution of mean xm.
public:
  PoissonDist(FLOAT mean)	// 
    : lambda(mean), oldm(-1.0) { ran1 = new UniformDist(); };
  LINT operator()()	// returns a random variable with Poisson dist.
    { return LINT( poidev(lambda) ); };
};


class NormalDist
{
private:
  FLOAT mu;	// mean
  FLOAT sigma;	// (std. deviation)^2
  int iset;
  float gset;
  UniformDist *ran1;

  float gasdev(void);
	// Returns a normally distributed deviate with zero mean and
	// unit variance.
public:
  NormalDist(FLOAT m, FLOAT s)	// mu, sigma
    : mu(m), sigma(s), iset(0) { ran1 = new UniformDist(); };
  FLOAT operator()()	// returns a random variable with Normal dist.
    { return gasdev() * sigma + mu; };
};


class ExpDist
{
private:
  float lambda;
  UniformDist *ran1;

  float expdev(void);
	// Returns an exponentially distributed, positive, random
	// deviate of unit mean.
public:
  ExpDist(FLOAT mean = 1.0)
    : lambda(mean) { ran1 = new UniformDist(); };
  FLOAT operator()()	// returns a random variable with an exp. distribution
    { return lambda * expdev(); };
};


// Used to choose k random numbers in the range [1..n]
//
class Choose
{
private:
  LINT *num;	// list of the positions
  FLOAT *rval;	// random value (used to get random ordering of the items)
public:
  Choose(LINT n, LINT k);
  ~Choose(void);
  LINT pos(LINT i) { return num[i]; };	// returns the i-th position
};
