#include <assert.h>
#include <iostream>

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

// write the plot:
// - into the output stream o
// - plot first n values from arr
// - with maximum value of height MAXHEIGHT
template<typename T>
void plot(std::ostream &o, T *arr, int n, const int MAXHEIGHT) {
    int maxv = 0;
    for(int i = 0; i < n; ++i) { assert(arr[i] >= 0); maxv = max(arr[i], maxv); }

    // (0, 0) --> +ve x
    // |
    // v
    // +ve y


    for(int y = 0; y < MAXHEIGHT; ++y) {
        for(int x = 0; x < n; ++x) {
            double rescaled = double(arr[x]) / maxv;
            int h = int(rescaled * MAXHEIGHT + 0.5);

            if (y == MAXHEIGHT-1-h) { o << "_ ";  }
            else if (y >= MAXHEIGHT-1-h) { o  << "| "; } else { o << "  "; }
        }
        o << "\n";
    }
}

// compute the historgram
// - of the first n values of arr
// - store the result in hist
// - with number of bins nbins.
// - ASSUMPTION: hist[0..nbins-1] has been PRE-ALLOCATED.
template<typename T>
void histogram (const T *arr, int n,
        T *hist, int nbins) {
    for(int i = 0; i < nbins; ++i) hist[i] = 0;

    assert(n > 0);
    T minv = arr[0]; T maxv = arr[0];
    for(int i = 1; i < n; ++i) {
        minv = min(minv, arr[i]);
        maxv = max(maxv, arr[i]);
    }


    // map (min, max) -> (0, nbins-1)
    // slope = y2-y1/x2-x1 = (nbins-1)/(max - min);
    // slope(x - min) = y - 0
    T slope = 0;
    if (maxv != minv) {
        slope = (nbins - 1) / (maxv - minv);
    }
    
    // std::cerr << "maxv: " << maxv << 
    //     " | minv: " << minv << 
    //     " | slope: " << slope << "\n";

    
    
    for(int i = 0; i < n; ++i) {
        int ix = int(slope * (arr[i] - minv));
        assert (ix >= 0);
        assert (ix < nbins);
        hist[ix] += 1;
    }
}
