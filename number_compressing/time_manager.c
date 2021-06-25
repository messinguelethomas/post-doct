#include"time_manager.h"


void top_(timeval_t* tv, timezone_t* tz)
{
        gettimeofday(tv, tz);
}

unsigned long get_temp_residuel(void)
{
   timeval_t _t2,_t1;
   timezone_t tz;
   gettimeofday(&_t1, &tz);
   gettimeofday(&_t2, &tz);


   return   1000000L * _t2.tv_sec + _t2.tv_usec -
           (1000000L * _t1.tv_sec + _t1.tv_usec );
}

unsigned long cpu_time_(timeval_t _t1,timeval_t _t2) /* retourne des microsecondes */
{
   unsigned long _temps_residuel =  get_temp_residuel();
   return 1000000L * _t2.tv_sec + _t2.tv_usec -
           (1000000L * _t1.tv_sec + _t1.tv_usec ) - _temps_residuel;
}
