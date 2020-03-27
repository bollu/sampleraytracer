#pragma once
#include <math.h>  // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdio.h>
#include <functional>
#include <iostream>
#include <vector>

static const int MAXRANDS = 1e4;

template <typename T>
struct Trace {
    double *rands;
    unsigned short Xi[3];
    float score = 0.0;
    int nrands = 0;
    int ix = 0;

    Trace(double *mem, unsigned short Xi[3]) : rands(mem) {
        for (int i = 0; i < 3; ++i) {
            this->Xi[i] = Xi[i];
        }
        score = 0;
        nrands = ix = 0;
    };

    Trace(const Trace &other) = delete;

    double rand() {
        assert(nrands < MAXRANDS);
        if (ix >= nrands) {
            const double r = erand48(Xi);
            rands[ix] = r;
            nrands = ix = ix + 1;
            return r;
        } else {
            return rands[ix++];
        }
    }
};

template <typename V>
std::ostream &operator<<(std::ostream &o, const Trace<V> &t) {
    return o << "Trace[score=" << t.score << "|nrands=" << t.nrands << "]";
}

// We need to copy the old code, because we nee to chain together these
// steps, so that we can sample "across" hastings calls. Otherwise we are
// not using the "warmed up" chain.
// Assumes that T.Xi has been setup already.
template <typename V, typename F, typename... Args>
V metropolisStep(Trace<V> &t, const int nmoves_per_sample, int &naccept, F f,
                 Args... args) {
    V prevv = f(t);
    double prevscore = t.score + log(t.nrands);
    int prev_nrands = t.nrands;

    for (int i = 1; i < nmoves_per_sample; i++) {
        // 1. perturb
        const unsigned int rix = t.nrands > 0 ? (nrand48(t.Xi) % t.nrands) : 0;
        const double prev_rand_at_rix = t.rands[rix];
        t.rands[rix] = erand48(t.Xi);

        // 2. sample
        t.score = 0;
        t.ix = 0;
        const V curv = f(t);
        const double curscore = t.score + log(t.nrands);
        const double acceptr = log(erand48(t.Xi));
        // TODO: double-check that this is indeed the correct
        // sampling criteria
        // 3. accept
        if (acceptr < curscore - prevscore) {
            prevv = curv;
            prevscore = curscore;
            prev_nrands = t.nrands;
            naccept++;
        } else {
            // use t = t (old trace)
            t.nrands = prev_nrands;
            t.rands[rix] = prev_rand_at_rix;
        }
    }
    return prevv;
};

template <typename V, typename F, typename... Args>
void sampleMH(int nsamples, int nmoves_per_sample, unsigned short Xi[3],
              double *randmem, V *out, F f, Args... args) {
    Trace<V> t(randmem, Xi);

    int naccept = 0;
    for (int i = 0; i < nsamples; ++i) {
        out[i] = metropolisStep(t, nmoves_per_sample, naccept, f, args...);
    }
}

// template <typename V, typename F, typename ...Args>
// V metropolisStep(Trace &t, const int nsamps, int &nmoves_per_sample,
//                  F f, Args... args) {
//                  // std::function<V(Trace &)> f) {
//
//     if (t.nrands > 0) {
//         const unsigned int rix = nrand48(t.Xi) % t.nrands;
//         t.rands[rix] = erand48(t.Xi);
//     }
//     // don't reset the trace. use the warmed up trace.
//     bool shouldquit = false;
//     V prevv(f(t, shouldquit));
//     if (shouldquit) { return prevv; }
//
//     double prevscore = 0;
//
//     for (int i = 1; i < nsamps; i++) {
//         // 1. perturb
//         const int prev_nrands = t.nrands;
//         const unsigned int rix = t.nrands > 0 ? nrand48(t.Xi) % t.nrands :
//         -1; const double prev_rand_at_rix = t.nrands  > 0?  t.rands[rix] :
//         -42;
//
//         t.rands[rix] = erand48(t.Xi);
//         // 2. sample
//         t.reset();
//         V curv = f(t, shouldquit,  args...);
//         if (shouldquit) { return curv; }
//         const double curscore = t.score + t.nrands;
//
//         const double acceptr = log(erand48(t.Xi));
//         // TODO: double-check that this is indeed the correct
//         // sampling criteria
//         const bool accept = acceptr < curscore - prevscore;
//         // 3. accept
//         if (accept) {
//             prevscore = curscore;
//             prevv = curv;
//         } else {
//             // use t = t (old trace)
//             t.nrands = prev_nrands;
//             if (rix > 0) { t.rands[rix] = prev_rand_at_rix; }
//         }
//     }
//     return prevv;
// };
