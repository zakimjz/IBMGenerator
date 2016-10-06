#include <math.h>
#include "dist.h"
#define PI 3.141592654

float PoissonDist::poidev(float xm)
// Returns as a floating-point number an integer value that is a
// random deviate drawn from a Poisson distribution of mean xm, using
// ran1(idum) as a source of uniform random deviates.
{
  float gammln(float xx);
  float em,t,y;
  
  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    t=1.0;
    do {
      em += 1.0;
      t *= (*ran1)();
    } while (t > g);
  } 
  else {
    if (xm != oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	y=tan(PI*(*ran1)());
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while ((*ran1)() > t);
  }
  return em;
}

#undef PI
