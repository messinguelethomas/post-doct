#ifndef __TIME_MANAGER_H
#define __TIME_MANAGER_H

#include <sys/time.h>
typedef struct timezone timezone_t;
typedef struct timeval timeval_t;


void top_(timeval_t* tv, timezone_t* tz);
unsigned long get_temp_residuel(void);
unsigned long cpu_time_(timeval_t _t1,timeval_t _t2);
#endif
