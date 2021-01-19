#pragma once

class Constant {
public:
    Constant();

    ~Constant();

    static int K;
    static bool IS_WEIGHTED;
    static int COMMUNITY_POPULATION;
    static const double PERCENTAGE_THRESHOLD;
    static const double EPSILON;
    static const double DELTA;
    static const int NUM_THREAD;
    static bool IS_BOUNDED_THRESHOLD;
    static const bool MODEL; // true : LT ; false : IC
    static const bool GCS;
    static const double EIG_TIME;
    static  int NUMBER_OF_COMMS;
};
