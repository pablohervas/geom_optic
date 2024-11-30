#include <time.h>

long double timeinseconds()
{
   long double seconds;
   struct timespec time;
   clock_gettime(CLOCK_MONOTONIC, &time);
// time.tv_sec is the number of seconds since last boot
// time.tv_nsec is the number of nanoseconds since time.tv_sec
// Hence, the time in seconds, with nine decimals of precision, is time.tv_sec + time.tv_nsec/1000000000
// A time in seconds with nine decimals of precision is equivalent to a time with a nanosecond of precision
   seconds=(long double)time.tv_sec + (long double)time.tv_nsec/1.e+09;
   return(seconds);
}

long long seed()
{
   long double ns;
   ns=timeinseconds()*1.e+09; // nanoseconds since epoch as a real number
   return (long long)ns; // nanoseconds since epoch as an integer number
}

double scalerandom(double rn, double a, double b)
{
   return a+(b-a)*rn;
}

