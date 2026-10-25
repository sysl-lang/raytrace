/* The showcase scene of `raytrace`, ported to C for the speed comparison in README.md.
 *
 * It is the sysl program's algorithm step for step -- the same PCG32 streams keyed by pixel and
 * sample, the same order of random draws, the same BVH rule, the same materials behind a table of
 * function pointers (the C spelling of a trait object) -- so the two images agree to within the
 * ties a differently ordered sort can break, and the timings compare like with like.
 *
 *   cc -O2 -o rt-c bench/rt.c -lm
 *   ./rt-c [width height spp] [--rays]     (single-threaded; writes c.ppm)
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct { double x, y, z; } V;
static inline V v3(double x, double y, double z) { V r = {x, y, z}; return r; }
static inline V add(V a, V b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline V sub(V a, V b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline V mul(V a, V b) { return v3(a.x * b.x, a.y * b.y, a.z * b.z); }
static inline V scl(V a, double k) { return v3(a.x * k, a.y * k, a.z * k); }
static inline V dvd(V a, double k) { return scl(a, 1.0 / k); }
static inline V neg(V a) { return v3(-a.x, -a.y, -a.z); }
static inline double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline V cross(V a, V b) { return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
static inline double len2(V a) { return dot(a, a); }
static inline V unit(V a) { return dvd(a, sqrt(len2(a))); }
static inline double axis(V a, int n) { return n == 0 ? a.x : n == 1 ? a.y : a.z; }
static inline V vmin(V a, V b) { return v3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z)); }
static inline V vmax(V a, V b) { return v3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z)); }
static inline V reflect(V v, V n) { return sub(v, scl(n, 2.0 * dot(v, n))); }
static inline V refract(V v, V n, double eta) {
    double c = fmin(-dot(v, n), 1.0);
    V perp = scl(add(v, scl(n, c)), eta);
    V para = scl(n, -sqrt(fabs(1.0 - len2(perp))));
    return add(perp, para);
}

/* PCG32, as sysl.rand has it */
typedef struct { uint64_t state, inc; } Rng;
static inline uint32_t next_u32(Rng *g) {
    uint64_t old = g->state;
    g->state = old * 6364136223846793005ULL + g->inc;
    uint32_t xs = (uint32_t)(((old >> 18) ^ old) >> 27);
    uint32_t rot = (uint32_t)(old >> 59);
    return (xs >> rot) | (xs << ((0u - rot) & 31));
}
static inline uint64_t next_u64(Rng *g) { uint64_t hi = next_u32(g); return hi << 32 | next_u32(g); }
static inline double unit01(Rng *g) { return (double)(next_u64(g) >> 11) / 9007199254740992.0; }
static Rng rng(uint64_t seed, uint64_t seq) {
    Rng g = {0, seq << 1 | 1};
    next_u32(&g); g.state += seed; next_u32(&g);
    return g;
}
static uint64_t below(Rng *g, uint64_t n) {
    uint64_t t = (0 - n) % n, r = next_u64(g);
    while (r < t) r = next_u64(g);
    return r % n;
}
static inline double uniform(Rng *g, double lo, double hi) { return lo + (hi - lo) * unit01(g); }
static V random_unit_vector(Rng *g) {
    for (;;) {
        double x = uniform(g, -1, 1), y = uniform(g, -1, 1), z = uniform(g, -1, 1);
        V p = v3(x, y, z);
        double l = len2(p);
        if (l > 1e-160 && l <= 1.0) return dvd(p, sqrt(l));
    }
}
static V random_in_unit_disk(Rng *g) {
    for (;;) {
        double x = uniform(g, -1, 1), y = uniform(g, -1, 1);
        V p = v3(x, y, 0);
        if (len2(p) < 1.0) return p;
    }
}

typedef struct { V o, d; double time; } Ray;
static inline V at(Ray r, double t) { return add(r.o, scl(r.d, t)); }
typedef struct { double lo, hi; } Iv;

/* textures */
typedef struct Tex Tex;
struct Tex { V (*value)(const Tex *, double u, double v, V p); };
typedef struct { Tex t; V c; } Solid;
static V solid_value(const Tex *t, double u, double v, V p) { return ((const Solid *)t)->c; }
static Tex *solid(V c) { Solid *s = malloc(sizeof *s); s->t.value = solid_value; s->c = c; return &s->t; }
typedef struct { Tex t; double inv; Tex *even, *odd; } Checker;
static V checker_value(const Tex *t, double u, double v, V p) {
    const Checker *c = (const Checker *)t;
    int x = (int)floor(c->inv * p.x), y = (int)floor(c->inv * p.y), z = (int)floor(c->inv * p.z);
    return (x + y + z) % 2 == 0 ? c->even->value(c->even, u, v, p) : c->odd->value(c->odd, u, v, p);
}
static Tex *checker(double scale, V a, V b) {
    Checker *c = malloc(sizeof *c); c->t.value = checker_value; c->inv = 1.0 / scale;
    c->even = solid(a); c->odd = solid(b); return &c->t;
}
/* the planet: assets/make_assets.py's 128x64 image, generated here with the same formulas
   (eight-bit sRGB, row-major from the top), sampled by (u, v) and converted to linear light */
#define PW 128
#define PH 64
typedef struct { Tex t; int w, h; unsigned char rgb[PW * PH * 3]; } Image;
static double pymod(double a, double b) { double r = fmod(a, b); return (r != 0 && (r < 0) != (b < 0)) ? r + b : r; }
static double to_linear(unsigned char b) { double s = (double)b / 255.0; return s * s; }
static V image_value(const Tex *t, double u, double v, V p) {
    const Image *m = (const Image *)t;
    double uu = u < 0 ? 0 : u > 1 ? 1 : u;
    double vv = 1.0 - (v < 0 ? 0 : v > 1 ? 1 : v);
    int i = (int)(uu * (double)m->w), j = (int)(vv * (double)m->h);
    if (i > m->w - 1) i = m->w - 1;
    if (j > m->h - 1) j = m->h - 1;
    const unsigned char *px = m->rgb + (j * m->w + i) * 3;
    return v3(to_linear(px[0]), to_linear(px[1]), to_linear(px[2]));
}
static Tex *planet(void) {
    Image *m = malloc(sizeof *m); m->t.value = image_value; m->w = PW; m->h = PH;
    const double pi = 3.141592653589793;
    for (int y = 0; y < PH; y++) {
        double lat = 90 - (y + 0.5) * 180 / PH;
        for (int x = 0; x < PW; x++) {
            double lon = (x + 0.5) * 360 / PW - 180;
            int land = sin(lon * (pi / 180) * 3) * cos(lat * (pi / 180) * 2) > 0.2;
            int r = land ? 60 : 20, g = land ? 140 : 60, b = land ? 50 : 170;
            if (fabs(lat) > 70) { r = 235; g = 240; b = 245; }
            if (fabs(pymod(lon, 30)) < 1.5 || fabs(pymod(lat, 30)) < 1.5) { r = 250; g = 220; b = 80; }
            unsigned char *px = m->rgb + (y * PW + x) * 3;
            px[0] = (unsigned char)r; px[1] = (unsigned char)g; px[2] = (unsigned char)b;
        }
    }
    return &m->t;
}
typedef struct { Tex t; V vecs[256]; int px[256], py[256], pz[256]; double scale; } Marble;
static void perm(Rng *g, int *p) {
    for (int i = 0; i < 256; i++) p[i] = i;
    for (size_t i = 256; i > 1;) { i--; size_t j = below(g, i + 1); int t = p[i]; p[i] = p[j]; p[j] = t; }
}
static double noise(const Marble *m, V p) {
    double fx = floor(p.x), fy = floor(p.y), fz = floor(p.z);
    double u = p.x - fx, v = p.y - fy, w = p.z - fz;
    int i = (int)fx, j = (int)fy, k = (int)fz;
    double uu = u * u * (3 - 2 * u), vv = v * v * (3 - 2 * v), ww = w * w * (3 - 2 * w), acc = 0;
    for (int di = 0; di < 2; di++)
        for (int dj = 0; dj < 2; dj++)
            for (int dk = 0; dk < 2; dk++) {
                V g = m->vecs[m->px[(i + di) & 255] ^ m->py[(j + dj) & 255] ^ m->pz[(k + dk) & 255]];
                double a = di, b = dj, c = dk;
                acc += (a * uu + (1 - a) * (1 - uu)) * (b * vv + (1 - b) * (1 - vv)) *
                    (c * ww + (1 - c) * (1 - ww)) * dot(g, v3(u - a, v - b, w - c));
            }
    return acc;
}
static double turb(const Marble *m, V p, int depth) {
    double acc = 0, w = 1; V q = p;
    for (int i = 0; i < depth; i++) { acc += w * noise(m, q); w *= 0.5; q = scl(q, 2.0); }
    return fabs(acc);
}
static V marble_value(const Tex *t, double u, double v, V p) {
    const Marble *m = (const Marble *)t;
    return scl(v3(0.5, 0.5, 0.5), 1.0 + sin(m->scale * p.z + 10.0 * turb(m, p, 7)));
}
static Tex *marble(uint64_t seed, double scale) {
    Marble *m = malloc(sizeof *m); m->t.value = marble_value; m->scale = scale;
    Rng g = rng(seed, 7);
    for (int i = 0; i < 256; i++) {
        double x = uniform(&g, -1, 1), y = uniform(&g, -1, 1), z = uniform(&g, -1, 1);
        m->vecs[i] = unit(v3(x, y, z));
    }
    perm(&g, m->px); perm(&g, m->py); perm(&g, m->pz);
    return &m->t;
}

/* materials */
typedef struct Mat Mat;
typedef struct { V p, n; const Mat *mat; double t, u, v; int front; } Hit;
typedef struct { V att; Ray ray; } Scatter;
struct Mat {
    int (*scatter)(const Mat *, Ray in, const Hit *h, Rng *g, Scatter *out);
    V (*emitted)(const Mat *, double u, double v, V p);
};
static V no_emission(const Mat *m, double u, double v, V p) { return v3(0, 0, 0); }
typedef struct { Mat m; Tex *tex; } Lambertian;
static int lambert_scatter(const Mat *m, Ray in, const Hit *h, Rng *g, Scatter *out) {
    V d = add(h->n, random_unit_vector(g));
    if (fabs(d.x) < 1e-8 && fabs(d.y) < 1e-8 && fabs(d.z) < 1e-8) d = h->n;
    const Lambertian *l = (const Lambertian *)m;
    out->att = l->tex->value(l->tex, h->u, h->v, h->p);
    out->ray = (Ray){h->p, d, in.time};
    return 1;
}
static Mat *lambertian(Tex *t) { Lambertian *l = malloc(sizeof *l); l->m.scatter = lambert_scatter; l->m.emitted = no_emission; l->tex = t; return &l->m; }
typedef struct { Mat m; V albedo; double fuzz; } Metal;
static int metal_scatter(const Mat *m, Ray in, const Hit *h, Rng *g, Scatter *out) {
    const Metal *me = (const Metal *)m;
    V r = add(unit(reflect(in.d, h->n)), scl(random_unit_vector(g), me->fuzz));
    if (dot(r, h->n) <= 0) return 0;
    out->att = me->albedo; out->ray = (Ray){h->p, r, in.time};
    return 1;
}
static Mat *metal(V a, double f) { Metal *m = malloc(sizeof *m); m->m.scatter = metal_scatter; m->m.emitted = no_emission; m->albedo = a; m->fuzz = f; return &m->m; }
typedef struct { Mat m; double ior; } Dielectric;
static double schlick(double c, double ri) { double r = (1 - ri) / (1 + ri), r0 = r * r; return r0 + (1 - r0) * pow(1 - c, 5.0); }
static int glass_scatter(const Mat *m, Ray in, const Hit *h, Rng *g, Scatter *out) {
    double ior = ((const Dielectric *)m)->ior, ri = h->front ? 1.0 / ior : ior;
    V u = unit(in.d);
    double c = fmin(-dot(u, h->n), 1.0), s = sqrt(1 - c * c);
    int total = ri * s > 1.0;
    V d = (total || schlick(c, ri) > unit01(g)) ? reflect(u, h->n) : refract(u, h->n, ri);
    out->att = v3(1, 1, 1); out->ray = (Ray){h->p, d, in.time};
    return 1;
}
static Mat *dielectric(double ior) { Dielectric *d = malloc(sizeof *d); d->m.scatter = glass_scatter; d->m.emitted = no_emission; d->ior = ior; return &d->m; }
typedef struct { Mat m; Tex *tex; } Light;
static int light_scatter(const Mat *m, Ray in, const Hit *h, Rng *g, Scatter *out) { return 0; }
static V light_emitted(const Mat *m, double u, double v, V p) { const Light *l = (const Light *)m; return l->tex->value(l->tex, u, v, p); }
static Mat *light(V c) { Light *l = malloc(sizeof *l); l->m.scatter = light_scatter; l->m.emitted = light_emitted; l->tex = solid(c); return &l->m; }

/* geometry */
typedef struct { Iv x, y, z; } Box;
static Iv iv_ordered(double a, double b) { return a <= b ? (Iv){a, b} : (Iv){b, a}; }
static Iv iv_pad(Iv i) { return i.hi - i.lo < 0.0001 ? (Iv){i.lo - 0.00005, i.hi + 0.00005} : i; }
static Box corners(V a, V b) {
    Box r = {iv_pad(iv_ordered(a.x, b.x)), iv_pad(iv_ordered(a.y, b.y)), iv_pad(iv_ordered(a.z, b.z))};
    return r;
}
static Iv iv_hull(Iv a, Iv b) { return (Iv){fmin(a.lo, b.lo), fmax(a.hi, b.hi)}; }
static Box hull(Box a, Box b) { return (Box){iv_hull(a.x, b.x), iv_hull(a.y, b.y), iv_hull(a.z, b.z)}; }
static Iv box_axis(Box b, int n) { return n == 1 ? b.y : n == 2 ? b.z : b.x; }
static int box_hit(const Box *b, Ray r, Iv t) {
    double lo = t.lo, hi = t.hi;
    for (int a = 0; a < 3; a++) {
        Iv ax = box_axis(*b, a);
        double inv = 1.0 / axis(r.d, a), o = axis(r.o, a);
        double t0 = (ax.lo - o) * inv, t1 = (ax.hi - o) * inv;
        if (t0 < t1) { if (t0 > lo) lo = t0; if (t1 < hi) hi = t1; }
        else { if (t1 > lo) lo = t1; if (t0 < hi) hi = t0; }
        if (hi <= lo) return 0;
    }
    return 1;
}

typedef struct Obj Obj;
struct Obj { int (*hit)(const Obj *, Ray r, Iv t, Hit *h); Box bbox; };
static void facing(Hit *h, Ray r, double t, V p, V out, const Mat *m, double u, double v) {
    int front = dot(r.d, out) < 0;
    h->p = p; h->n = front ? out : neg(out); h->mat = m; h->t = t; h->u = u; h->v = v; h->front = front;
}
typedef struct { Obj o; Ray c; double r; const Mat *m; } Sphere;
static int sphere_hit(const Obj *ob, Ray r, Iv ti, Hit *h) {
    const Sphere *s = (const Sphere *)ob;
    V c = at(s->c, r.time), oc = sub(c, r.o);
    double a = len2(r.d), hh = dot(r.d, oc), disc = hh * hh - a * (len2(oc) - s->r * s->r);
    if (disc < 0) return 0;
    double sq = sqrt(disc), root = (hh - sq) / a;
    if (!(ti.lo < root && root < ti.hi)) { root = (hh + sq) / a; if (!(ti.lo < root && root < ti.hi)) return 0; }
    V p = at(r, root), out = dvd(sub(p, c), s->r);
    double theta = acos(-out.y), phi = atan2(-out.z, out.x) + 3.141592653589793;
    facing(h, r, root, p, out, s->m, phi / (2 * 3.141592653589793), theta / 3.141592653589793);
    return 1;
}
static Obj *sphere(V c0, V c1, double rad, const Mat *m) {
    Sphere *s = malloc(sizeof *s); s->o.hit = sphere_hit; s->c = (Ray){c0, sub(c1, c0), 0}; s->r = rad; s->m = m;
    V rv = v3(rad, rad, rad);
    s->o.bbox = hull(corners(sub(c0, rv), add(c0, rv)), corners(sub(c1, rv), add(c1, rv)));
    return &s->o;
}
typedef struct { Obj o; V q, u, v, n, w; double d; const Mat *m; } Quad;
static int quad_hit(const Obj *ob, Ray r, Iv ti, Hit *h) {
    const Quad *q = (const Quad *)ob;
    double den = dot(q->n, r.d);
    if (fabs(den) < 1e-8) return 0;
    double t = (q->d - dot(q->n, r.o)) / den;
    if (!(ti.lo <= t && t <= ti.hi)) return 0;
    V p = at(r, t), pl = sub(p, q->q);
    double al = dot(q->w, cross(pl, q->v)), be = dot(q->w, cross(q->u, pl));
    if (al < 0 || al > 1 || be < 0 || be > 1) return 0;
    facing(h, r, t, p, q->n, q->m, al, be);
    return 1;
}
static Obj *quad(V Q, V u, V v, const Mat *m) {
    Quad *q = malloc(sizeof *q); q->o.hit = quad_hit; q->q = Q; q->u = u; q->v = v; q->m = m;
    V n = cross(u, v); q->n = unit(n); q->d = dot(q->n, Q); q->w = dvd(n, dot(n, n));
    q->o.bbox = hull(corners(Q, add(add(Q, u), v)), corners(add(Q, u), add(Q, v)));
    return &q->o;
}
typedef struct { Obj o; V a, e1, e2, n; const Mat *m; } Tri;
static int tri_hit(const Obj *ob, Ray r, Iv ti, Hit *h) {
    const Tri *t = (const Tri *)ob;
    V pv = cross(r.d, t->e2);
    double det = dot(t->e1, pv);
    if (fabs(det) < 1e-12) return 0;
    double inv = 1.0 / det;
    V tv = sub(r.o, t->a);
    double u = dot(tv, pv) * inv;
    if (u < 0 || u > 1) return 0;
    V qv = cross(tv, t->e1);
    double v = dot(r.d, qv) * inv;
    if (v < 0 || u + v > 1) return 0;
    double tt = dot(t->e2, qv) * inv;
    if (!(ti.lo < tt && tt < ti.hi)) return 0;
    facing(h, r, tt, at(r, tt), t->n, t->m, u, v);
    return 1;
}
static Obj *tri(V a, V b, V c, const Mat *m) {
    Tri *t = malloc(sizeof *t); t->o.hit = tri_hit; t->a = a; t->e1 = sub(b, a); t->e2 = sub(c, a);
    t->n = unit(cross(t->e1, t->e2)); t->m = m;
    t->o.bbox = hull(corners(a, add(a, t->e1)), corners(a, add(a, t->e2)));
    return &t->o;
}
typedef struct { Obj o; Obj *l, *r; } Node;
static int node_hit(const Obj *ob, Ray r, Iv ti, Hit *h) {
    const Node *n = (const Node *)ob;
    if (!box_hit(&n->o.bbox, r, ti)) return 0;
    int hl = n->l->hit(n->l, r, ti, h);
    Iv t2 = {ti.lo, hl ? h->t : ti.hi};
    int hr = n->r->hit(n->r, r, t2, h);
    return hl || hr;
}
static int sort_axis;
static int by_axis(const void *a, const void *b) {
    double x = box_axis((*(Obj *const *)a)->bbox, sort_axis).lo, y = box_axis((*(Obj *const *)b)->bbox, sort_axis).lo;
    return x < y ? -1 : x > y ? 1 : 0;
}
static Obj *build(Obj **xs, size_t n) {
    if (n == 1) return xs[0];
    Box b = {{INFINITY, -INFINITY}, {INFINITY, -INFINITY}, {INFINITY, -INFINITY}};
    for (size_t i = 0; i < n; i++) b = hull(b, xs[i]->bbox);
    Node *nd = malloc(sizeof *nd); nd->o.hit = node_hit; nd->o.bbox = b;
    if (n == 2) { nd->l = xs[0]; nd->r = xs[1]; return &nd->o; }
    double sx = b.x.hi - b.x.lo, sy = b.y.hi - b.y.lo, sz = b.z.hi - b.z.lo;
    sort_axis = sx > sy ? (sx > sz ? 0 : 2) : (sy > sz ? 1 : 2);
    qsort(xs, n, sizeof *xs, by_axis);
    nd->l = build(xs, n / 2); nd->r = build(xs + n / 2, n - n / 2);
    return &nd->o;
}

/* the scene: `showcase` in scenes.sysl, item for item */
static Obj *scene(void) {
    static Obj *items[16]; size_t n = 0;
    items[n++] = sphere(v3(0, -1000, 0), v3(0, -1000, 0), 1000, lambertian(checker(0.5, v3(0.2, 0.3, 0.1), v3(0.9, 0.9, 0.9))));
    items[n++] = sphere(v3(-2.2, 1, 0), v3(-2.2, 1, 0), 1, lambertian(solid(v3(0.7, 0.3, 0.2))));
    items[n++] = sphere(v3(0, 1, 0), v3(0, 1, 0), 1, dielectric(1.5));
    items[n++] = sphere(v3(0, 1, 0), v3(0, 1, 0), -0.8, dielectric(1.5));
    items[n++] = sphere(v3(2.2, 1, 0), v3(2.2, 1, 0), 1, metal(v3(0.8, 0.8, 0.85), 0.05));
    items[n++] = sphere(v3(-1.4, 0.5, 2), v3(-1.4, 0.5, 2), 0.5, lambertian(marble(42, 4.0)));
    items[n++] = sphere(v3(1.3, 0.6, 2.1), v3(1.3, 0.6, 2.1), 0.6, lambertian(planet()));
    items[n++] = sphere(v3(0.1, 0.35, 2.4), v3(0.1, 0.6, 2.4), 0.35, lambertian(solid(v3(0.2, 0.4, 0.9))));

    Mat *gold = metal(v3(0.9, 0.7, 0.3), 0.2);
    FILE *f = fopen("assets/icosphere.obj", "r");
    if (!f) { perror("assets/icosphere.obj"); exit(1); }
    V vs[64]; size_t nv = 0; Obj **tris = malloc(128 * sizeof *tris); size_t nt = 0;
    char line[256];
    while (fgets(line, sizeof line, f)) {
        double x, y, z; int a, b, c;
        if (sscanf(line, "v %lf %lf %lf", &x, &y, &z) == 3) vs[nv++] = add(scl(v3(x, y, z), 0.55), v3(-0.5, 0.55, -2.0));
        else if (sscanf(line, "f %d %d %d", &a, &b, &c) == 3) tris[nt++] = tri(vs[a - 1], vs[b - 1], vs[c - 1], gold);
    }
    fclose(f);
    items[n++] = build(tris, nt);
    items[n++] = quad(v3(-1.5, 4.5, -1), v3(3, 0, 0), v3(0, 0, 2), light(v3(6, 6, 5)));
    return build(items, n);
}

typedef struct { int w, h, spp, depth; V bg, center, p00, du, dv, du_disk, dv_disk; double defocus; uint64_t seed; } Cam;
static Cam camera(int w, int h, int spp) {
    Cam c = {w, h, spp, 20, v3(0.5, 0.6, 0.75)};
    V from = v3(0, 2.2, 7.5), to = v3(0, 0.8, 0), vup = v3(0, 1, 0);
    double vfov = 35, defocus = 0.6, focus = 7.6, pi = 3.141592653589793;
    double hh = tan(vfov * (pi / 180) / 2), vh = 2 * hh * focus, vw = vh * (double)w / (double)h;
    V W = unit(sub(from, to)), U = unit(cross(vup, W)), Vv = cross(W, U);
    c.du = scl(U, vw / w); c.dv = scl(neg(Vv), vh / h);
    V ul = add(sub(sub(from, scl(W, focus)), scl(U, vw / 2)), scl(Vv, vh / 2));
    c.p00 = add(ul, scl(add(c.du, c.dv), 0.5));
    double rad = focus * tan(defocus / 2 * (pi / 180));
    c.center = from; c.du_disk = scl(U, rad); c.dv_disk = scl(Vv, rad); c.defocus = defocus; c.seed = 1;
    return c;
}
static V trace(const Cam *c, Ray r, const Obj *w, Rng *g) {
    V color = v3(0, 0, 0), thr = v3(1, 1, 1);
    for (int d = 0; d < c->depth; d++) {
        Hit h;
        if (!w->hit(w, r, (Iv){0.001, INFINITY}, &h)) return add(color, mul(thr, c->bg));
        color = add(color, mul(thr, h.mat->emitted(h.mat, h.u, h.v, h.p)));
        Scatter s;
        if (!h.mat->scatter(h.mat, r, &h, g, &s)) return color;
        thr = mul(thr, s.att); r = s.ray;
    }
    return color;
}
static V pixel(const Cam *c, int i, int j, const Obj *w) {
    V sum = v3(0, 0, 0);
    for (int n = 0; n < c->spp; n++) {
        Rng g = rng(c->seed + (uint64_t)n, (uint64_t)(j * c->w + i));
        double ox = unit01(&g) - 0.5, oy = unit01(&g) - 0.5;
        V s = add(add(c->p00, scl(c->du, i + ox)), scl(c->dv, j + oy));
        V o = c->center;
        if (c->defocus > 0) { V d = random_in_unit_disk(&g); o = add(add(c->center, scl(c->du_disk, d.x)), scl(c->dv_disk, d.y)); }
        Ray r = {o, sub(s, o), unit01(&g)};
        sum = add(sum, trace(c, r, w, &g));
    }
    return dvd(sum, c->spp);
}
static unsigned char to_byte(double x) {
    if (!(x > 0)) return 0;
    double s = sqrt(x);
    return (unsigned char)(256 * (s < 0 ? 0 : s > 0.999 ? 0.999 : s));
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }

int main(int argc, char **argv) {
    int w = 400, h = 225, spp = 64, rays = 0;
    if (argc >= 4) { w = atoi(argv[1]); h = atoi(argv[2]); spp = atoi(argv[3]); }
    for (int i = 1; i < argc; i++) if (!strcmp(argv[i], "--rays")) rays = 1;
    Obj *world = scene();
    Cam c = camera(w, h, spp);
    double t0 = now();
    if (rays) {
        int hits = 0;
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++) {
                V p = add(add(c.p00, scl(c.du, i)), scl(c.dv, j));
                Hit hh;
                hits += world->hit(world, (Ray){c.center, sub(p, c.center), 0}, (Iv){0.001, INFINITY}, &hh);
            }
        double s = now() - t0;
        printf("%d primary rays, %d hits, %.0f ms, %.3f Mrays/s\n", w * h, hits, s * 1000, w * h / s / 1e6);
        return 0;
    }
    unsigned char *img = malloc((size_t)w * h * 3);
    for (int j = 0; j < h; j++)
        for (int i = 0; i < w; i++) {
            V p = pixel(&c, i, j, world);
            unsigned char *o = img + ((size_t)j * w + i) * 3;
            o[0] = to_byte(p.x); o[1] = to_byte(p.y); o[2] = to_byte(p.z);
        }
    double s = now() - t0;
    FILE *f = fopen("c.ppm", "wb");
    fprintf(f, "P6\n%d %d\n255\n", w, h); fwrite(img, 1, (size_t)w * h * 3, f); fclose(f);
    printf("showcase %dx%d at %d spp on 1 thread: %.0f ms -> c.ppm\n", w, h, spp, s * 1000);
    return 0;
}
