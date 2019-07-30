#ifndef TIMER_H
#define TIMER_H

#include <time.h>


class Timer {
    private:
        clock_t start;
    
    public:
        inline Timer() {
            start = clock();
        }

        inline int reset() {	
            start = clock();
        }

        inline float get() {
            return float( (long double)(clock() - start) / 
                          (long double)CLOCKS_PER_SEC );
        }

        inline float lap() {
            float ret = get();
            reset();
            return ret;
        }
};

#endif
