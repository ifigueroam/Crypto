/*  These are for time measures  */

/*
 *  Adjust this to the clock freq. of the machine where time intervals
 *  are computed. 
 *
 *  cpu MHz value in file /proc/cpuinfo
 */
#define CLOCK 2530.000

typedef struct {
        unsigned long hi;
        unsigned long lo;
} time_686;

#define time(x) \
__asm__ __volatile__("rdtsc" : "=d" (x.hi), "=a" (x.lo))

static inline double time_diff(time_686 b, time_686 a) {
        double db,da,res;
        db = (double)b.hi*(double)(1<<16)*(double)(1<<16)+(double)b.lo;
        da = (double)a.hi*(double)(1<<16)*(double)(1<<16)+(double)a.lo;
        if (db < da)
res = ((double)(1<<16)*(double)(1<<16)*(double)(1<<16)*(double)(1<<16)+db-da)
/(double)CLOCK;
	else
res = (db-da)/(double)CLOCK;
	return res;
}

