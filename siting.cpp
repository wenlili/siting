////////////////////////////////////////////////////////////////////////////////
// utilities
////////////////////////////////////////////////////////////////////////////////

#include <limits.h>  // INT_MIN
#include <math.h>    // sqrt
#include <stdlib.h>  // abs, atoi
#include <sys/time.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

typedef unsigned short int usint;
typedef unsigned long long int ullint;

#include <omp.h>

void die(const char *msg) {
  cout << "ERROR: " << msg << endl;
  exit(1);
}

#define RANDOM_MAX 0x7fffffff  // 2^31 - 1

inline int random(const int seed) {
  // glibc: m = 2^31; a = 1103515245; c = 12345
  return (int)((1103515245U * ((unsigned)seed & 0x7fffffffU) + 12345U) &
               0x7fffffffU);
}

template <class C>
inline C sign(const C x) {
  return x < 0 ? -1 : (x == 0 ? 0 : 1);
}

template <class C>
inline C square(const C x) {
  return x * x;
}

inline double calc_time(struct timeval &begin, struct timeval &end) {
  return ((end.tv_sec - begin.tv_sec) * 1000000u +
          end.tv_usec - begin.tv_usec) / 1e6;
}

////////////////////////////////////////////////////////////////////////////////
// vix
////////////////////////////////////////////////////////////////////////////////

int test_one_target(const int nrows, const usint *elevs,
                    const int ox, const int oy, const int oz,
                    const int tx, const int ty, const int tz) {
  //const float oz = elevs[ox * nrows + oy] + oht;
  //const float tz = elevs[tx * nrows + ty] + tht;
  int d, id, stride = 1;
  float z, drecip, slope, zslope;
  const int iteroverx = (abs(tx - ox) >= abs(ty - oy));
  if (abs(ox - tx) <= 1 && abs(oy - ty) <= 1) return 1;

  //      iteroverx, s, sign(slope):
  //
  //             \0,1,-1|      /
  //              \     |0,1,1/
  //               \    |    /
  //    1,-1,-1     \   |   /  1,1,1
  //                 \  |  /
  //                  \3|2/
  //                  4\|/1
  //  ------------------O-------------------
  //                  5/|\8
  //    1,-1,1        /6|7\    1,1,-1
  //                 /  |  \
  //                /   |   \
  //               /    |    \
  //              /     |     \
  //             /0,-1,1|      \
  //            /       |0,-1,-1\
  //
  // zslope is the slope of the LOS from the observer to the target.
  // If the slope from the observer to any intermediate point on the
  // LOS is larger, then the target is hidden.

  if (iteroverx) {
    d = tx - ox;
    if (d > 0) {
      // if (d <= 1) die("Illegal d, 1st and 8th octants");
      drecip = 1.0 / d;
      slope =  (ty - oy) * drecip;  // may be -ve.
      zslope = (tz - oz) * drecip;
      for (id = 1; id < d; id += stride, stride <<= 1) {
        z = elevs[(ox + id) * nrows + (oy + (int)round(slope * id))];
        if (z > oz + id * zslope) return 0;
        // if ((z - oz) / id > zslope) return 0;
      }
    } else {
      // if (d >= -1) die("Illegal d, 4th and 5th octants");
      drecip = 1.0 / d;
      slope =  (ty - oy) * drecip;  // may be -ve.
      zslope = (tz - oz) * drecip;
      for (id = -1; id > d; id -= stride, stride <<= 1) {
        z = elevs[(ox + id) * nrows + (oy + (int)round(slope * id))];
        if (z > oz + id * zslope) return 0;
        // if ((z - oz) / id < zslope) return 0;
      }
    }
  } else {  // Iterating over y
    d = ty - oy;
    if (d > 0) {
      // if (d <= 1) die("Illegal d in 2nd and 3rd octants");
      drecip = 1.0 / d;
      slope =  (tx - ox) * drecip;  // may be -ve.
      zslope = (tz - oz) * drecip;
      for (id = 1; id < d; id += stride, stride <<= 1) {
        z = elevs[(ox + (int)round(slope * id)) * nrows + (oy + id)];
        if (z > oz + id * zslope) return 0;
        // if ((z - oz) / id) > zslope) return 0;
      }
    } else {
      // if (d >= -1) die("Illegal d in 6th and 7th octants");
      drecip = 1.0 / d;
      slope =  (tx - ox) * drecip;  // may be -ve.
      zslope = (tz - oz) * drecip;
      for (id = -1; id > d; id -= stride, stride <<= 1) {
        z = elevs[(ox + (int)round(slope * id)) * nrows + (oy + id)];
        if (z > oz + id * zslope) return 0;
        // if ((z - oz) / id) < zslope) return 0;
      }
    }
  }
  return 1;
}

void calc_one_vix(const int tid,
                  const int nrows, const usint *elevs,
                  const int roi, const int oht, const int tht,
                  const int ntests, unsigned char *vix) {
  const int ox = tid / nrows;
  const int oy = tid % nrows;
  const int oz = elevs[ox * nrows + oy] + oht;
  int ntarget = 0;
  int nvis = 0;
  const int vsxmax = min(nrows - 1, ox + roi);  // viewshed bounds
  const int vsymax = min(nrows - 1, oy + roi);
  const int vsxmin = max(0, ox - roi);
  const int vsymin = max(0, oy - roi);
  int r = tid;
  // for (int i = 0; i < 10 * ntests; i++) {  // iterate over random targets
  while (ntarget < ntests) {
    int tx, ty, tz;
    int visq;
    r = random(r);
    //tx = (int)(r * (vsxmax - vsxmin + 0.99999f) / RANDOM_MAX) + vsxmin;
    tx = (int)((2*roi+0.99999f)*r/RANDOM_MAX) + (ox-roi);
    r = random(r);
    //ty = (int)(r * (vsymax - vsymin + 0.99999f) / RANDOM_MAX) + vsymin;
    ty = (int)((2*roi+0.99999f)*r/RANDOM_MAX) + (oy-roi);
    // if (tx == 0 && ty == 0) continue;
    if (square(tx - ox) + square(ty - oy) > square(roi)) {
      tx = ox + (tx - ox)/3; // golden
      ty = oy + (ty - oy)/3;
      //continue; // consistent
    }
    if (tx >=0 && tx < nrows && ty >= 0 && ty < nrows) {
        tz = elevs[tx * nrows + ty] + tht;
        visq = test_one_target(nrows, elevs, ox, oy, oz, tx, ty, tz);
        // cout << "test_one_target(" << ox << ',' << oy << ',' << tx << ',' << ty << ")=" << visq << endl;
    } else {
        visq = 0;
    }
    ntarget++;
    if (visq) nvis++;

    // Stopping rule: do at least 10 points.
    // Then, continue until vix >= .5 or <= .1.
    // This could be improved by using the variance, and by looking at
    // other observers, to select the best.

    // if (ntarget >= ntests) break;
    // if (ntarget < 10) continue;
    // if (v >= 0.5 || v <= 0.1) break;
  }
  float v = (float)nvis / ntarget;
  // cout << "obs at (" << ox << ',' << oy << "), z=" << elevs[ox * nrows + oy] 
  //      << ", vix=" << nvis << '/' << ntarget << '=' << v << endl;
  vix[ox * nrows + oy] = (unsigned char)min(255, (int)(v * 255.999f));
}

void calc_vix(const int nrows, const usint *elevs,
              const int roi, const int oht, const int tht,
              const int ntests, unsigned char *vix) {
  #pragma omp parallel for schedule(guided)
  for (int tid = 0; tid < square(nrows); tid++)
    calc_one_vix(tid, nrows, elevs, roi, oht, tht, ntests, vix);
}

////////////////////////////////////////////////////////////////////////////////
// findmax
////////////////////////////////////////////////////////////////////////////////

void process_one_block(const int bx, const int by, const int nrows,
                       unsigned char *vix, const float blocksize,
                       const int nwantedperblock, int *obs) {
  const int xmin = (int)(blocksize * bx);
  const int xmax = min((int)(blocksize * (bx + 1)), nrows);
  const int ymin = (int)(blocksize * by);
  const int ymax = min((int)(blocksize * (by + 1)), nrows);
  const int width = ymax - ymin;
  const int npoints = (xmax - xmin) * width;
  for (int i = 0; i < nwantedperblock; i++) {
    int p1x = xmin;
    int p1y = ymin;
    unsigned char v1 = vix[p1x * nrows + p1y];
    int h1 = p1x * (p1x + p1y) * 010101010101;
    for (int j = 1; j < npoints; j++) {
      int p2x = xmin + j / width;
      int p2y = ymin + j % width;
      unsigned char v2 = vix[p2x * nrows + p2y];
      int h2 = p2x * (p2x + p2y) * 010101010101;
      if (v1 < v2 || (v1 == v2 && h1 < h2)) {
        p1x = p2x;
        p1y = p2y;
        v1 = v2;
        h1 = h2;
      }
    }
    obs[i * 2] = p1x;
    obs[i * 2 + 1] = p1y;
    vix[p1x * nrows + p1y] = 0;  // -1 is better
  }
}

void find_max(const int nrows, unsigned char *vix, const float blocksize, const int /*nwanted*/,
              const int nblocks, const int nwantedperblock, int *obs) {
  #pragma omp parallel for schedule(guided)
  for (int i = 0; i < square(nblocks); i++)
    process_one_block(i/nblocks, i%nblocks, nrows, vix, blocksize,
                      nwantedperblock, &obs[i*nwantedperblock*2]);
}

////////////////////////////////////////////////////////////////////////////////
// viewshed
////////////////////////////////////////////////////////////////////////////////

void set_vis(const int nwpr, const int row, const int col,
             ullint *shed) {
  shed[row * nwpr + col / 64] |= 1ULL << (63 - col % 64);
}

void calc_one_shed(const int tid,
                   const int nrows, const usint *elevs,
                   const int roi, const int oht, const int tht,
                   const int nsheds, const int *obs,
                   ullint *sheds) {
  const int observer[2] = {obs[2 * tid], obs[2 * tid + 1]};
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;
  ullint *thisshed = &sheds[tid * nr * nwpr];

  int delta[2];
  int target[2];
  int p[2];  // Current point
  int sig, pelev;
  float horizon_slope, slope, horizon_alt;
  int observer_alt;
  int inciny;

  // Clipping xmin etc at 0, nrows-1 makes the viewshed depend on the roi, so don't.
  const int xmin = observer[0] - roi;
  const int ymin = observer[1] - roi;
  const int xmax = observer[0] + roi;
  const int ymax = observer[1] + roi;
  const int xwidth = xmax - xmin;
  const int ywidth = ymax - ymin;
  const int perimeter = 2 * (xwidth + ywidth);  // This formula is subtle

  set_vis(nwpr, roi, roi, thisshed);  // Observer is visible from itself.

  // Observer distance above sea level, including distance above ground.
  observer_alt = elevs[observer[0] * nrows + observer[1]] + oht;

  // The target is in turn every point along the smaller of the border or a box
  // of side 2 * radius around the observer.

  // xmax etc are coords of pixels, not of the edges between the pixels.  I.e.,
  // xmin=5, xmax=7 means 3 pixels.
  // A 3x3 regions has a perimeter of 8.

  for (int ip = 0; ip < perimeter; ip++) {
    if (ip < xwidth) {
      target[0] = xmin + ip;
      target[1] = ymin;
    } else if (ip < 2 * xwidth) {
      target[0] = 1 + xmin - xwidth + ip;
      target[1] = ymax;
    } else if (ip < 2 * xwidth + ywidth) {
      target[0] = xmin;
      target[1] = 1 + ymin - 2 * xwidth + ip;
    } else {
      target[0] = xmax;
      target[1] = ymin - 2 * xwidth - ywidth + ip;
    }

    // don't clip
    // if (target[0] < 0) target[0] = 0;
    // if (target[0] >= nrows) target[0] = nrows - 1;
    // if (target[1] < 0) target[1] = 0;
    // if (target[1] >= nrows) target[1] = nrows - 1;

    // This occurs only when observer is on the edge of the region.
    if (observer[0] == target[0] && observer[1] == target[1]) continue;

    // Run a line of sight out from obs to target.
    delta[0] = target[0] - observer[0];
    delta[1] = target[1] - observer[1];
    inciny = (abs(delta[0]) < abs(delta[1]));  // outer parens reqd?

    // Step along the coord (X or Y) that varies the most from the observer to
    // the target.  Inciny says which coord that is.  Slope is how fast the
    // other coord varies.
    slope = (float)delta[1 - inciny] / (float)delta[inciny];
    sig = (delta[inciny] > 0 ? 1 : -1);
    horizon_slope = -99999.f;  // Slope (in vertical plane) to horizon so far.

    // i=0 would be the observer, which is always visible.
    for (int i = sig; i != delta[inciny] + sig; i += sig) {
      p[inciny] = observer[inciny] + i;
      p[1 - inciny] = observer[1 - inciny] + (int)round(i * slope);

      // Have we reached the edge of the area?
      if (p[0] < 0 || p[0] >= nrows || p[1] < 0 || p[1] >= nrows) break;
      if ((square(p[0] - observer[0]) + square(p[1] - observer[1]) > square(roi))) break;

      pelev = elevs[p[0] * nrows + p[1]];
      // Slope from the observer, incl the observer_ht, to this point, at ground
      // level.  The slope is projected into the plane XZ or YZ, depending on
      // whether X or Y is varying faster, and thus being iterated thru.
      float s = (float)(pelev - observer_alt) / (float)abs(p[inciny] - observer[inciny]);
      if (horizon_slope < s) horizon_slope = s;
      horizon_alt = observer_alt + horizon_slope * abs(p[inciny] - observer[inciny]);
      if (pelev + tht >= horizon_alt)
        set_vis(nwpr, p[0] - observer[0] + roi, p[1] - observer[1] + roi, thisshed);
    }
  }
}

void calc_sheds(const int nrows, const usint *elevs,
                const int roi, const int oht, const int tht,
                const int nsheds, const int *obs,
                ullint *sheds) {
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;
  for (int i = 0; i < nsheds * nr * nwpr; i++)
    sheds[i] = 0ULL;
  #pragma omp parallel for schedule(guided)
  for (int tid = 0; tid < nsheds; tid++)
    calc_one_shed(tid, nrows, elevs, roi, oht, tht, nsheds, obs, sheds);
}

////////////////////////////////////////////////////////////////////////////////
// site
////////////////////////////////////////////////////////////////////////////////

const int BIT_COUNT[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};

inline int popcll(ullint i) {
  int sum = 0;
  const unsigned char *c;
  c = reinterpret_cast<const unsigned char *>(&i);
  for (int j = 0; j < 8; j++)
    sum += BIT_COUNT[c[j]];
  return sum;
}

void calc_union(const int nrows, const int roi,
                const int *observer, const ullint *shed,
                ullint *cumshed) {
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;
  const int cumnwpr = (nrows + 63) / 64;

  int firstword = (observer[1] - roi) / 64;
  int firstbit = (observer[1] - roi) % 64;
  if (firstbit < 0) {
    firstword--;
    firstbit += 64;
  }

  #pragma omp parallel for schedule(guided)
  for (int row = 0; row < nr; row++) {  // each row of shed
    int cumrow = observer[0] - roi + row;
    if (cumrow >= 0 && cumrow < nrows)  // row inside terrain
      for (int word = 0; word < nwpr; word++) {  // each word of row
        int cumword = firstword + word;
        if (cumword >= 0 && cumword < cumnwpr)  // word inside terrain
          cumshed[cumrow * cumnwpr + cumword] |= shed[row * nwpr + word] >> firstbit;
        if (firstbit != 0 && cumword + 1 >= 0 && cumword + 1 < cumnwpr)  // firstbit != 0 and word + 1 inside terrain
          cumshed[cumrow * cumnwpr + cumword + 1] |= shed[row * nwpr + word] << (64 - firstbit);
      }
  }
}

int is_obs_vis(const int cumnwpr, const ullint *cumshed,
               const int *observer) {
  int i = observer[0] * cumnwpr * 64 + observer[1];
  if ((cumshed[i / 64] & 1ULL << (63 - i % 64)) != 0)
    return 1;
  else
    return 0;
}

void union_area(const int bid, const int nrows, const int roi,
                const int intervis, const int nsheds, const int *obs,
                const ullint *sheds,
                const int nr, const int nwpr, const int cumnwpr,
                const int nusedsheds, const char *usedq,
                const int lastobsx, const int lastobsy,
                const ullint *cumshed, int *testshedarea) {
  //if (usedq[bid]) return;

  int observer[2];
  observer[0] = obs[2 * bid];
  observer[1] = obs[2 * bid + 1];

  //if (nusedsheds > 0 && square(observer[0] - lastobsx) + square(observer[1] - lastobsy) > square(2 * roi)) return;

  int *areas = 0;
  areas = new int[nr];
  if (!areas) die("out of memory");
  for (int i = 0; i < nr; i++)
    areas[i] = 0;

  int visible;
  if (intervis && nusedsheds > 0)
    visible = is_obs_vis(cumnwpr, cumshed, observer);

  //#pragma omp parallel for schedule(guided)
  for (int row = 0; row < nr; row++)
    if (!intervis || nusedsheds == 0 || visible) {  // calculate one row of extra area
      const int cumrow = observer[0] - roi + row;
      if (cumrow >= 0 && cumrow < nrows) {
        int firstword = (observer[1] - roi) / 64;
        int firstbit = (observer[1] - roi) % 64;
        if (firstbit < 0) {
          firstword--;
          firstbit += 64;
        }
        int lastword = (observer[1] + roi) / 64;

        ullint prevvalue = 0ULL;
        ullint value, cumvalue, tempvalue;
        int sum = 0;

        for (int cumword = firstword; cumword <= lastword; cumword++)
          if (cumword >= 0 && cumword < cumnwpr) {
            int word = cumword - firstword;  // definition out of loop?
            if (cumword == 0 && word > 0) prevvalue = sheds[bid * nr * nwpr + row * nwpr + word - 1];
            if (word < nwpr)
              value = sheds[bid * nr * nwpr + row * nwpr + word];
            else
              value = 0ULL;
            cumvalue = cumshed[cumrow * cumnwpr + cumword];
            tempvalue = cumvalue;
            tempvalue |= value >> firstbit;
            if (firstbit != 0) tempvalue |= prevvalue << (64 - firstbit);
            tempvalue ^= cumvalue;
            sum += popcll(tempvalue);
            prevvalue = value;
          }

        areas[row] = sum;
      }
    }

  int sum = 0;
  for (int i = 0; i < nr; i++)
    sum += areas[i];
  testshedarea[bid] = sum;

  delete[] areas;
}

void site_it(const int nrows, const int roi, const int intervis,
             const int nsheds, const int *obs, const ullint *sheds, char *selected) {
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;
  const int cumnwpr = (nrows + 63) / 64;

  int *usedsheds;                     // list of sheds used so far.
  char *usedq;                        // whether each particular shed has been used.
  ullint *cumshed;    // Cumulative viewshed.
  int *testshedarea;

  usedsheds = 0;
  usedq = 0;
  cumshed = 0;
  testshedarea = 0;
  int *areas = 0;
  usedsheds = new int[nsheds];
  usedq = new char[nsheds];
  cumshed = new ullint[nrows * cumnwpr];
  testshedarea = new int[nsheds];
  areas = new int[nsheds];
  if (!usedsheds || !usedq || !cumshed || !testshedarea || !areas)
    die("Memory exhausted. Program terminates.");

  for (int i = 0; i < nsheds; i++) usedq[i] = 0;
  for (int i = 0; i < nrows * cumnwpr; i++) cumshed[i] = 0ULL;
  //cout << "Total area=" << square(nrows) << endl;
  //cout << "#nusedsheds newshed obsx obsy area extraarea newcumarea areapercentage alreadyvisible" << endl;

  int nusedsheds = 0;
  int lastobsx = 0;
  int lastobsy = 0;
  int cumarea = 0;
  for (int i = 0; i < nsheds; i++)
    testshedarea[i] = 0;

  while (1) {
    vector<int> updatelist;
    #pragma omp parallel for schedule(guided)
    for (int bid = 0; bid < nsheds; bid++)
      if (!usedq[bid])
        if (nusedsheds == 0 || square(obs[bid*2]-lastobsx) + square(obs[bid*2+1]-lastobsy) <= square(2*roi))
          #pragma omp critical
          updatelist.push_back(bid);
    #pragma omp parallel for schedule(guided)
    for (int i = 0; i < updatelist.size(); i++)
      union_area(updatelist[i], nrows, roi, intervis, nsheds, obs, sheds,
                 nr, nwpr, cumnwpr, nusedsheds, usedq, lastobsx, lastobsy,
                 cumshed, testshedarea);

    int newshed = -1;
    int extraarea = 0;
#ifdef SEQ
    // Sequential version.
    for (int i = 0; i < nsheds; i++) {
      if (usedq[i])
        continue;
      if (testshedarea[i] > extraarea) {
        extraarea = testshedarea[i];
        newshed = i;
      }
    }
#else
    // OpenMP version.
    #pragma omp parallel
    {
      int nthreads = omp_get_num_threads();
      int size = (nsheds + nthreads - 1) / nthreads;
      int id = omp_get_thread_num();
      int start = id * size;
      int end = (id + 1) * size;
      if (end > nsheds) end = nsheds;
      int m_newshed = -1;
      int m_extraarea = 0;
      for (int i = start; i < end; i++) {
        if (usedq[i]) continue;
        if (testshedarea[i] > m_extraarea) {
          m_extraarea = testshedarea[i];
          m_newshed = i;
        }
      }
      #pragma omp critical
      if (m_extraarea > extraarea) {
        extraarea = m_extraarea;
        newshed = m_newshed;
      }
    }
#endif
    if (extraarea == 0) {
      //cout << "No more new observers that will add new area." << endl;
      break;
    }

    //int alreadyvisible = is_obs_vis(cumnwpr, cumshed, &obs[2 * newshed]);  // Must calc before updating cumshed.
    calc_union(nrows, roi, &obs[2 * newshed], &sheds[newshed * nr * nwpr], cumshed);
    usedsheds[nusedsheds++] = newshed;
    usedq[newshed] = 1;
    lastobsx = obs[newshed * 2];
    lastobsy = obs[newshed * 2 + 1];

    cumarea += extraarea;
    double areapercentage = 100.0 * cumarea / square(nrows);
    if (areapercentage > 95)
        break;

    if (nusedsheds == 1)
      for (int i = 0; i < nsheds; i++)
        areas[i] = testshedarea[i];

    /*
    cout << setw(6) << nusedsheds << setw(6) << newshed
         << setw(6) << obs[newshed * 2] << setw(6) << obs[newshed * 2 + 1]
         << setw(8) << areas[newshed] << setw(8) << extraarea
         << setw(10) << cumarea << setw(8) << areapercentage << setw(2) << alreadyvisible << endl;
    */
  }
  cout << " nusedsheds:" << nusedsheds << " coverage:" << 100.0*cumarea/square(nrows);

  for (int i = 0; i < nsheds; i++)
    selected[i] = usedq[i];
  delete[] usedsheds;
  delete[] usedq;
  delete[] cumshed;
  delete[] testshedarea;
  delete[] areas;
}

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
#ifndef SEQ
  omp_set_dynamic(0);
  omp_set_num_threads(50);
#endif

  double elapsed_secs;
  struct timeval begin, middle, end;
  gettimeofday(&begin, NULL);
  middle = begin;

  int nrows;                      // # rows, cols in this cell
  int roi;                        // radius of interest
  int oht;                        // ht of observer above terrain
  int tht;                        // target ht above terrain
  int ntests;                     // # of targets tested per observer
  int blocksize0;                 // Requested number of rows in one block of the input visibility index array.
  float blocksize;               // Perturbed blocksize0, to remove sliver blocks at the end.
  int nwanted0;                   // Desired number of output observers.
  int nwanted;                    // Modified because blocksize was perturbed.
  int nblocks;                    // Number of blocks on one side of it.
  int nwantedperblock;            // Number of observers to find per block, before culling.
  int intervis;                   // Should new observers be visible to existing ones?

  usint *elevs;      // terrain elevation (input)
  unsigned char *vix;             // visibility index * 256 (output)
  int *obs;                       // observers
  ullint *sheds;  // viewsheds
  char *selected;

  //cout << "[SITE, compiled from " << __FILE__ << " on " << __DATE__ << ", " << __TIME__ << ']' << endl;
  if (argc != 10) {
    cout << "argc=" << argc << endl;
    die("SITE requires 9 arguments: nrows, roi, oht/tht, ntests, blocksize, nwanted, intervis, infile, outfile");
  }

  nrows = atoi(argv[1]);
  roi = atoi(argv[2]);
  oht = tht = atoi(argv[3]);
  ntests = atoi(argv[4]);
  blocksize0 = atoi(argv[5]);
  nwanted0 = atoi(argv[6]);
  intervis = atoi(argv[7]);

  //cout << "nrows=" << nrows << ", roi=" << roi << ", oht=" << oht << ", tht=" << tht
  //     << ", ntests=" << ntests << "\nblocksize0=" << blocksize0 << ", nwanted0=" << nwanted0
  //     << ", intervis=" << intervis << endl;

  if (nrows <= 0 || nrows > 20000) die("Unreasonable value for nrows");
  if (roi < 1 || roi > 10000) die("Unreasonable value for roi.");
  if (tht < 0 || tht > 1000000) die("Unreasonable value for tht.");
  if (ntests < 1 || ntests > 1000) die("Unreasonable value for ntests.");
  if (blocksize0 < 10 || blocksize0 > 2000) die("Unreasonable value for blocksize0.");
  if (nwanted0 < 100 || nwanted0 > 2000000) die("Unreasonable value for nwanted0.");
  if (intervis != 0 && intervis != 1) die("Unreasonable value for intervis.");

  // Perturb blocksize so that the last block won't be really small.
  blocksize = (float)nrows / (int)((float)nrows / blocksize0 + 0.5f);  // floating point block size
  nblocks = (int)(nrows / blocksize + 0.5f);                             // number of blocks
  nwantedperblock = (int)((float)nwanted0 / square(nblocks) + 0.99999f);  // number of wanted per block
  nwanted = nwantedperblock * square(nblocks);                          // number of wanted
  int lastsize = nrows - (int)(blocksize * (nblocks - 1));              // size of the last block
  if (square(lastsize) < nwantedperblock)                               // too small
    die("The last block is too small for nwantedperblock.");

  //cout << "blocksize=" << blocksize << ", nblocks=" << nblocks
  //     << ", nwantedperblock=" << nwantedperblock << ", nwanted=" << nwanted << endl;

  // number of rows per shed and number of words per row
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;

  elevs = 0;
  vix = 0;
  obs = 0;
  sheds = 0;
  selected = 0;
  elevs = new usint[square(nrows)];
  vix = new unsigned char[square(nrows)];
  obs = new int[2 * nwanted];
  sheds = new ullint[nwanted * nr * nwpr];
  selected = new char[nwanted];
  if (!elevs || !vix || !obs || !sheds || !selected)
    die("Memory exhausted. Program terminates.");
  /*
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nrows; j++) {
      cin.read(reinterpret_cast<char *>(&elevs[i * nrows + j]), sizeof(usint));
      if (cin.fail()) {
        cout << "Error: at i=" << i << ", j=" << j << endl;
        die("Input failed");
      }
    }
  */
  ifstream ifs(argv[8]);
  ifs.read((char *)elevs, square(nrows)*sizeof(usint));
  if (ifs.fail())
    die("Input failed");
  ifs.close();

  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cout << "input:" << elapsed_secs;

  calc_vix(nrows, elevs, roi, oht, tht, ntests, vix);
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cout << " vix:" << elapsed_secs;

  find_max(nrows, vix, blocksize, nwanted, nblocks, nwantedperblock, obs);
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cout << " findmax:" << elapsed_secs;
  
  calc_sheds(nrows, elevs, roi, oht, tht, nwanted, obs, sheds);
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cout << " viewshed:" << elapsed_secs;

  site_it(nrows, roi, intervis, nwanted, obs, sheds, selected);
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cout << " site:" << elapsed_secs;

  ofstream ofs(argv[9]);
  for (int i = 0; i < nwanted; i++)
    if (selected[i])
      ofs << obs[2*i] << ',' << obs[2*i+1] << '\n';
  ofs.close();

  delete[] elevs;
  delete[] vix;
  delete[] obs;
  delete[] sheds;
  delete[] selected;

  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cout << " output:" << elapsed_secs;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(begin, end);
  cout << " total:" << elapsed_secs << endl;
}
