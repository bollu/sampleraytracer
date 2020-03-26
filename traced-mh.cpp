#include <assert.h>
#include <math.h>    // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdio.h>   //        Remove "-fopenmp" for g++ version < 4.2
#include <stdlib.h>  // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <string.h>  // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include "trace.h"

static const bool LOG = false;
const int NTRACES = 4;
double *tracemem;

template <typename T>
T max(T a, T b) {
    return a > b ? a : b;
}
template <typename T>
T max3(T a, T b, T c) {
    return max(a, max(b, c));
}

template <typename T>
T min(T a, T b) {
    return a < b ? a : b;
}
template <typename T>
T min3(T a, T b, T c) {
    return min(a, min(b, c));
}

struct Vec {         // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z;  // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) {
        x = x_;
        y = y_;
        z = z_;
    }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const {
        return x * b.x + y * b.y + z * b.z;
    }  // cross:
    Vec operator%(Vec &b) {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};
struct Ray {
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
    double rad;   // radius
    Vec p, e, c;  // position, emission, color
    Refl_t refl;  // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_)
        : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const {  // returns distance, 0 if nohit
        Vec op = p - r.o;  // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t, eps = 1e-4, b = op.dot(r.d),
                  det = b * b - op.dot(op) + rad * rad;
        if (det < 0)
            return 0;
        else
            det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

const Vec zero = Vec();
Sphere spheres[] = {
    // Scene: radius, position, emission, color, material
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), zero, Vec(.75, .25, .25),
           DIFF),  // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), zero, Vec(.25, .25, .75),
           DIFF),                                                     // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), zero, Vec(.75, .75, .75), DIFF),  // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), zero, zero, DIFF),         // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), zero, Vec(.75, .75, .75), DIFF),  // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), zero, Vec(.75, .75, .75),
           DIFF),                                                      // Top
    Sphere(16.5, Vec(27, 16.5, 47), zero, Vec(1, 1, 1) * .999, SPEC),  // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), zero, Vec(1, 1, 1) * .999, SPEC),  // Glas
    Sphere(500, Vec(50, 81.6 + 500 - .27, 81.6), Vec(20, 20, 20), zero,
           DIFF)  // Lite
};
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray &r, double &t, int &id) {
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t) {
            t = d;
            id = i;
        }
    return t < inf;
}

// https://www.rapidtables.com/convert/color/rgb-to-hsv.html
Vec rgb2hsv(Vec rgb) {
    float cmax = max3(rgb.x, rgb.y, rgb.z);
    float cmin = min3(rgb.x, rgb.y, rgb.z);
    float d = cmax - cmin;
    return Vec(-42, cmax == 0 ? 0 : d / cmax, cmax);
}

Vec radiance(const Ray &r, int depth, Trace &rands) {
    double t;    // distance to intersection
    int id = 0;  // id of intersected object
    if (!intersect(r, t, id)) {
        // rands.score += log(1e-3);
        return Vec();  // if miss, return black
    }

    const Sphere &obj = spheres[id];  // the hit object
    Vec x = r.o + r.d * t, n = (x - obj.p).norm(),
        nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
    double p =
        f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;  // max refl
    if (++depth > 5) {
        if (rands.rand() < p) {
            f = f * (1 / p);
        } else {
            // rands.score += log1p(0);
            // rands.score += log(10 * rgb2hsv(obj.e).z);
            return obj.e; /*r.r*/
        }
    }
    if (obj.refl == DIFF) {  // Ideal DIFFUSE reflection
        // rands.score += log(1.1);
        double r1 = 2 * M_PI * rands.rand(), r2 = rands.rand(), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(),
            v = w % u;
        Vec d =
            (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, rands));
    } else if (obj.refl == SPEC) {  // Ideal SPECULAR reflection
        //rands.score += log(1.1);
        // send more towards SPECULAR
        return obj.e +
               f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, rands));
    } else {
        assert(obj.refl == REFR);
        // rands.score += log(1.1);
        Ray reflRay(x,
                    r.d - n * 2 * n.dot(r.d));  // Ideal dielectric REFRACTION
        bool into = n.dot(nl) > 0;              // Ray from outside going in?
        double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc,
               ddn = r.d.dot(nl), cos2t;
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) <
            0)  // Total internal reflection
            return obj.e + f.mult(radiance(reflRay, depth, rands));
        Vec tdir =
            (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t))))
                .norm();
        double a = nt - nc, b = nt + nc, R0 = a * a / (b * b),
               c = 1 - (into ? -ddn : tdir.dot(n));
        double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re,
               P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
        return obj.e +
               f.mult(depth > 2
                          ? (rands.rand() < P
                                 ?  // Russian roulette
                                 radiance(reflRay, depth, rands) * RP
                                 : radiance(Ray(x, tdir), depth, rands) * TP)
                          : radiance(reflRay, depth, rands) * Re +
                                radiance(Ray(x, tdir), depth, rands) * Tr);
    }
}



int main(int argc, char *argv[]) {
    int w = 512, h = 512,
        samps = argc == 2 ? atoi(argv[1]) : 1;              // # samples
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm());  // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r,
        *c = new Vec[w * h];

    tracemem = new double[MAXRANDS];

#pragma omp parallel for schedule(dynamic, 1) private(r)  // OpenMP
    for (int y = 0; y < h; y++) {  // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4,
                100. * y / (h - 1));
        // cols: why does he use same initialization?
        for (unsigned short x = 0; x < w; x++) {
            // 2x2 subpixel rows
            for (int sy = 0; sy < 1; sy++) {
                const int i = (h - y - 1) * w + x;
                // 2x2 subpixel cols
                for (int sx = 0; sx < 1; sx++) {
                    // Vec r = Vec();

                    for (int s = 0; s < samps/NTRACES; s++) {
                        int naccept;
                        const Vec r = metropolisStep<Vec>(
                            tracemem,
                            NTRACES, naccept, [&](Trace &trace) {
                                double r1 = 2 * trace.rand(),
                                       dx = r1 < 1 ? sqrt(r1) - 1
                                                   : 1 - sqrt(2 - r1);
                                double r2 = 2 * trace.rand(),
                                       dy = r2 < 1 ? sqrt(r2) - 1
                                                   : 1 - sqrt(2 - r2);
                                Vec d =
                                    cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                    cy * (((sy + .5 + dy) / 2 + y) / h - .5) +
                                    cam.d;
                                Vec rad =
                                    radiance(Ray(cam.o + d * 140, d.norm()), 0,
                                             trace) *
                                    (1. / (samps / NTRACES));
                                const Vec hsv = rgb2hsv(rad);
                                // not sure this is correct. Usually. you still
                                // need to return the sample on an incorrect
                                // step.
                                return rad;
                            });

                        c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));

                    }
                    
                    // Camera rays are pushed ^^^^^ forward to start in
                    // interior
                    // c[i] = c[i] + Vec(clamp(r.x * .25), clamp(r.y *.25), clamp(r.z * .25));
                }  // sx
            }      // sy
        }          // x

    }                                       // y
    FILE *f = fopen("traced-mh.ppm", "w");  // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}

