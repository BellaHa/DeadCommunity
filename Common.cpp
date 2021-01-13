#include "Common.h"
#include <algorithm>

Common *Common::instance = nullptr;

Common::Common() {
    seed = new int[100];
    for (int i = 0; i < 100; i++) {
        seed[i] = rand();
    }
}

Common::~Common() {
}

Common *Common::getInstance() {
    if (instance == nullptr)
        instance = new Common();
    return instance;
}

// long Common::nChoosek(long N, long K) {
//     // This function gets the total number of unique combinations based upon N and K.
//     // N is the total number of items.
//     // K is the size of the group.
//     // Total number of unique combinations = N! / ( K! (N - K)! ).
//     // This function is less efficient, but is more likely to not overflow when N and K are large.
//     // Taken from:  http://blog.plover.com/math/choose.html
//     //
//     if (K > N) return 0;
//     long r = 1;
//     long d;
//     for (d = 1; d <= K; d++) {
//         r *= N--;
//         r /= d;
//     }
//     return r;
// }

unsigned Common::nChoosek(unsigned n, unsigned k) {
    if (k > n) return 0;
    if (k == 0) return 1;
    unsigned re = n;
    for (int i = 2; i <= k; i++) {
        re *= (n - i + 1);
        re /= i;
    }
    return re;
}

bool Common::isIntersected(vector<int> *set1, vector<int> *set2) {
    int j = 0;
    for (int i = 0; i < set1->size(); i++) {
        while (j < set2->size() && (*set2)[j] < (*set1)[i]) {
            j++;
        }
        if (j < set2->size() && (*set2)[j] == (*set1)[i])
            return true;
        else if (j == set2->size())
            return false;
    }
    return false;
}

vector<int> Common::setDifference(vector<int> *set1, vector<int> *set2) {
    vector<int> re;
    for (int i = 0; i < set1->size(); i++) {
        if (find(set2->begin(), set2->end(), set1->at(i)) == set2->end()) {
            re.push_back(set1->at(i));
        }
    }
    return re;
}

unsigned Common::randomInThread() {
    unsigned tmp = seed[omp_get_thread_num() % 100];
    tmp = tmp * 17931 + 7391;
    seed[omp_get_thread_num() % 100] = tmp;
    return tmp;
}

double Common::lognCk(unsigned int n, unsigned int k) {
    // double _n = 0, _k = 0, n_k = 0;
    // for (int i = 1; i <= n; i++) _n += log(1.0 * i); // n!
    // for (int i = 1; i <= k; i++) _k += log(1.0 * i); // k!
    // for (int i = 1; i <= n - k; i++) n_k += log(1.0 * i); // n-k!
    // return _n - _k - n_k;

    if (k > n) return 0.;
    if (k == 0) return 1.;
    double up = 0., down = 0.;
    for (int d = 1; d <= k; d++) {
        up += log(1. * n);
        n--;
        down += log(1. * d);
    }
    return up - down;
}
