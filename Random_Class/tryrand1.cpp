#define WANT_STREAM
#define WANT_MATH
#define WANT_TIME

#include "include.h"
#include "newran.h"

//#ifdef use_namespace
//using namespace NEWRAN;
//#endif



Real phi(Real x)                          // normal density
{ return (fabs(x)>8.0) ? 0 : 0.398942280 * exp(-x*x / 2); }

Real NORMAL10(Real x) { return 0.5*phi(0.5*(x-10.0)); }

Real UNIF(Real x) { return (x>=0.0 && x<=1.0) ? 1.0 : 0.0; }

Real TRIANG(Real x) { return (x>1.0) ? 0.0 : 1.0-x; }



void main(int n)
{
   Normal nn;
   Uniform uniform;
   cout << "Print 20 N(0,1) random numbers - should be the same as in sample output" << endl;
   { for (int i=0; i<20; i++) cout << nn.Next() << "\n" << flush; }

   return 0;
}
