#include <string.h>
#include "stdlib.h"
#include "assert.h"
#include "trace.h"
#include "plot.h"
#include <iostream>
using namespace std;


double gaussian(Trace<double> &trace) {
    double r = 0.5 - trace.rand();
    r *= 5;
    trace.score = -r*r;
    return r;
};

static const int MAXHEIGHT = 15;
static const int NSAMPLES = 1e5;
static const int NMOVES_PER_SAMPLE = 20;
static const int NBINS = 40;

int plot_gaussian() {
    std::cout << "##Gaussian:##\n";
    double *randmem = new double[MAXRANDS];
    double gs[NSAMPLES];

    unsigned short Xi[3] = {0, 0, 1};
    sampleMH<double>(NSAMPLES, NMOVES_PER_SAMPLE, Xi, randmem, gs, gaussian);

    double gshist[NBINS];
    histogram<double>(gs, NSAMPLES, gshist, NBINS);

    plot(std::cout, gshist, NBINS, MAXHEIGHT);
    return 0;
}

int plot_uniform() {
    std::cout << "##Uniform:##\n";
    double *randmem = new double[MAXRANDS];
    double gs[NSAMPLES];
    unsigned short Xi[3] = {0, 0, 1};
    sampleMH<double>(NSAMPLES, NMOVES_PER_SAMPLE, Xi, randmem, gs, [](Trace<double> &trace) {
            trace.score = 0; return trace.rand();
            });

    double gshist[NBINS];
    histogram<double>(gs, NSAMPLES, gshist, NBINS);
    plot(std::cout, gshist, NBINS, MAXHEIGHT);
    return 0;
}

int plot_modsin() {
    std::cout << "##|sin|:##\n";
    double *randmem = new double[MAXRANDS];
    double gs[NSAMPLES];

    unsigned short Xi[3] = {0, 0, 1};
    sampleMH<double>(NSAMPLES, NMOVES_PER_SAMPLE, Xi, randmem, gs, 
            [](Trace<double> &trace) { double r = 20 * trace.rand(); trace.score = fabs(sin(r)); return r; });

    double gshist[NBINS];
    histogram<double>(gs, NSAMPLES, gshist, NBINS);

    plot(std::cout, gshist, NBINS, MAXHEIGHT);
    return 0;
}

int plot_floor() {
    std::cout << "##ceil:##\n";
    double *randmem = new double[MAXRANDS];
    double gs[NSAMPLES];

    unsigned short Xi[3] = {0, 0, 1};
    sampleMH<double>(NSAMPLES, NMOVES_PER_SAMPLE, Xi, randmem, gs, 
            [](Trace<double> &trace) { double r = 100 * trace.rand(); assert(r >= 0); trace.score = fabs(ceil(r)); return r; });

    double gshist[NBINS];
    histogram<double>(gs, NSAMPLES, gshist, NBINS);

    plot(std::cout, gshist, NBINS, MAXHEIGHT);
    return 0;
}


int main() {
    plot_gaussian();
    plot_uniform();
    plot_modsin();
    plot_floor();
}
