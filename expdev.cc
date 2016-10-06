#include <math.h>
#include "dist.h"

float ExpDist::expdev(void)
// Returns an exponentially distributed, positive, random deviate of
// unit mean, using ran1(idum) as the source of uniform deviates.
{
  float dum;
  
  do
    dum = (*ran1)();
  while (dum == 0.0);
  return -log(dum);
}
