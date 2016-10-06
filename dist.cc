#include "dist.h"
#include <limits.h>

//------------------------------- RandSeed -------------------------------

UniformDist *RandSeed::ran1 = new UniformDist(INIT_SEED);

void RandSeed::set_seed(long new_seed)
{
  delete ran1;
  ran1 = new UniformDist(new_seed);
}

long RandSeed::new_seed(void)
{
  float ans;

  ans = (*ran1)();
  ans *= (*ran1)();
//  cout << "rand_seed : " << idum << " " << long(-LONG_MAX * ans) << endl;
  return long(-LONG_MAX * ans);
}


//------------------------------- Choose -------------------------------


// allows selection of k random items from the string
//
Choose::Choose(LINT n, LINT k)
{
  UniformDist ud;
  LINT i, j, ival;
  FLOAT fval;

  num = new LINT [n];
  rval = new FLOAT [n];

  // associate a random value with each item
  // also copy item into num
  for (i = 0; i < n; i++)
    {
      rval[i] = ud();
      num[i] = i;
    }

  // sort num according to the values in rval
  for (i = 0; i < n; i++ )
    {
      ival = num[i]; fval = rval[i];
      for ( j = i; j > 0 && rval[j-1] > fval; j-- ) {
	  num[j] = num[j-1];
	  rval[j] = rval[j-1];
	}
      num[j] = ival;
      rval[j] = fval;
    }

  // resort first k num according to position
  for (i = 0; i < k; i++ )
    {
      ival = num[i];
      for ( j = i; j > 0 && num[j-1] > ival; j-- )
	num[j] = num[j-1];
      num[j] = ival;
    }
}


Choose::~Choose(void)
{
  delete [] num;
  delete [] rval;
}


