#define WANT_STREAM
#define WANT_TIME

#include "include.h"
#include "newran.h"
#include "tryrand.h"

#ifdef use_namespace
using namespace NEWRAN;
#endif

int main()
{

   Random::Set(0.46875);

Normal Z;
cout << Z.Next() << endl;

   return 0;
}





