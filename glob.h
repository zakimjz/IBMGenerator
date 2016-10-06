// Some global declarations

//const int MAXINT       =  2147483647;
//const float FLOATMAX   =  3.40282346638528860e+38;
//const double MAXDOUBLE =  7.23700557733226211e+75;

#define NIL 0

#include <stdlib.h>

typedef char CHAR;
typedef unsigned char UCHAR;
typedef short int SINT;
typedef unsigned short USINT;
typedef int LINT;
typedef unsigned int ULINT;
typedef float FLOAT;
typedef double DOUBLE;
typedef double DOUBLEPRECISION;
typedef UCHAR BOOLEAN;

#define TRUE 1  //True and false were not defined in the original -- MJZaki
#define FALSE 0
//#define max(a, b) ((a) > (b) ? (a) : (b))
//#define min(a, b) ((a) < (b) ? (a) : (b))

// useful types
typedef char *Filename;
typedef LINT Item;
typedef LINT Tid;
typedef LINT Cid;
