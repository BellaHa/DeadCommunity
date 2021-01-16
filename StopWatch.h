#pragma once

#include <chrono>

class StopWatch {
private:
    chrono::steady_clock::time_point mStart, mStop;
    double milliseconds = 0;

public:
    StopWatch() = default;

    ~StopWatch() = default;

    void start() {
        mStart = chrono::steady_clock::now();
    }

    void stop() {
        mStop = chrono::steady_clock::now();
    }

    void reset() {
        mStart = mStop = chrono::steady_clock::now();
    }

    double getSeconds() {
        getMilliSeconds();
        return milliseconds / 1000;
    }

    double getMilliSeconds() {
        milliseconds = (double) chrono::duration_cast<chrono::milliseconds>(mStop - mStart).count();
        reset();
        return milliseconds;
    }
};
