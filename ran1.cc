#include "dist.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// "Minimal" random number generator of Park and Miller with
// Bays-Durham shuffle and added safeguards. Returns a uniform random
// deviate between 0.0 and 1.0 (exclusive of the endpoint values).
// Call with idum a negative integer to initialize; thereafter, do not
// alter idum between successive deviates in a sequence. RNMX should
// approximate the largest flating value that is less than 1.
float UniformDist::ran1(void)
{
  int j;
  long k;
  float temp;
  
  if (idum <= 0 || !iy) {	// Initialize
    if (-idum < 1)		// Be sure to prevent idum = 0
      idum = 1;
    else
      idum = -idum;
    for (j = NTAB+7; j >= 0; j-- )  // Load the shuffle table (after 8 warmups)
      {
	k = idum/IQ;
	idum = IA * (idum-k*IQ) - IR*k;
	if (idum < 0) idum += IM;
	if (j < NTAB) iv[j] = idum;
      }
    iy = iv[0];
  }
  k = idum/IQ;			  // Start here when not initializing.
  idum = IA * (idum-k*IQ) - IR*k; // Compute idum = (IA*idum) % IM without
  if (idum < 0) idum += IM;	  //     overflows by Schrage's method.
  j = iy/NDIV;			  // Will be in the range 0..NTAB-1
  iy = iv[j];			  // Output previously stored value and
  iv[j] = idum;			  //      and refull the shuffle table.
  if ((temp = AM*iy) > RNMX) 	  // Because users don't expect endpoint values
    return RNMX;
  else
    return temp;

}

