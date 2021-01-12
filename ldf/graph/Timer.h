/*
 * Timer.h
 *
 *  Created on: Jan 16, 2012
 *      Author: thang
 */

#ifndef TIMER_H_
#define TIMER_H_

#ifdef	_WIN32
	#define WIN32
#endif

#ifdef WIN32
	#include <windows.h>
	typedef LARGE_INTEGER TIME_TYPE;
#else
	#include <sys/time.h>
	typedef timeval TIME_TYPE;
#endif


class Timer
{
private:

TIME_TYPE frequency, startCount, endCount;

public:
    /**
     * Default constructor. The counter is started automatically.
     */
    Timer() {
#ifdef WIN32
		QueryPerformanceFrequency(&frequency);
		startCount.QuadPart = 0;
		endCount.QuadPart = 0;
#else
		startCount.tv_sec = startCount.tv_usec = 0;
		endCount.tv_sec = endCount.tv_usec = 0;
#endif
    	start();
    }

    /**
     * Start/reset the time counter
     */
    void   start() {
#ifdef WIN32
    	QueryPerformanceCounter(&startCount);
#else
    	gettimeofday(&startCount, NULL);
#endif
    }

    /**
     * @return	Elapsed time in second(s).
     */
    double stop() {
#ifdef WIN32
    	QueryPerformanceCounter(&endCount);
    	return ( -startCount.QuadPart * (1.0 / frequency.QuadPart) +
    	     endCount.QuadPart * (1.0 / frequency.QuadPart));
#else
    	gettimeofday(&endCount, NULL);
    	double x = startCount.tv_sec + startCount.tv_usec*0.000001;
    	double y = endCount.tv_sec  + endCount.tv_usec*0.000001;
    	return y - x;
#endif
    }
};

#endif /* TIMER_H_ */
