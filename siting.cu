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

using namespace std;

typedef unsigned short int usint;
typedef unsigned long long int ullint;

void die(const char *msg) {
  cerr << "ERROR: " << msg << endl;
  exit(1);
}

#define RANDOM_MAX 0x7fffffff  // 2^31 - 1

__device__ int random(const int seed) {
  // glibc: m = 2^31; a = 1103515245; c = 12345
  return (int)((1103515245U * ((unsigned)seed & 0x7fffffffU) + 12345U) &
               0x7fffffffU);
}

template <class C>
__device__ C sign(const C x) {
  return x < 0 ? -1 : (x == 0 ? 0 : 1);
}

template <class C>
__host__ __device__ C square(const C x) {
  return x * x;
}

inline double calc_time(struct timeval &begin, struct timeval &end) {
  return ((end.tv_sec - begin.tv_sec) * 1000000u +
          end.tv_usec - begin.tv_usec) / 1e6;
}

////////////////////////////////////////////////////////////////////////////////
// vix
////////////////////////////////////////////////////////////////////////////////

__device__ int test_one_target(const int nrows, const usint *elevs,
                               const int ox, const int oy, const int oz,
                               const int tx, const int ty, const int tz) {
  if (abs(ox - tx) <= 1 && abs(oy - ty) <= 1) return 1;
  int dx = tx - ox, dy = ty - oy;
  int px, py, pz;  // Current point
  int inciny = abs(dx) < abs(dy);
  int sign;
  float slope, zslope;
  sign = (inciny*dy + (1-inciny)*dx) > 0 ? 1 : -1;
  slope = (float)(inciny*dx + (1-inciny)*dy) / (inciny*dy + (1-inciny)*dx);
  zslope = (float)(tz - oz) / (inciny ? dy : dx);
  const int limit = inciny ? dy : dx;
  int stride = 1;
  for (int i = sign; abs(i) < abs(limit); i += stride*sign, stride <<= 1) {  // *= 1.9, 2.2, 2.5
    int j = round(i * slope);
    px = ox + (inciny*j + (1-inciny)*i);
    py = oy + (inciny*i + (1-inciny)*j);
    pz = elevs[px * nrows + py];
    if (pz > oz + i * zslope) return 0;
  }
  return 1;
}

__global__ void calc_one_vix(const int nrows, const usint *elevs,
                             const int roi, const int oht, const int tht,
                             const int ntests, unsigned char *vix) {
  const int bid = blockIdx.y * gridDim.x + blockIdx.x;
  const int tid = bid * blockDim.x + threadIdx.x;
  if (tid >= square(nrows)) return;

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
      //continue; // too slow
    }
    if (tx >=0 && tx < nrows && ty >= 0 && ty < nrows) {
        tz = elevs[tx * nrows + ty] + tht;
        visq = test_one_target(nrows, elevs, ox, oy, oz, tx, ty, tz);
        // cerr << "test_one_target(" << ox << ',' << oy << ',' << tx << ',' << ty << ")=" << visq << endl;
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
  // cerr << "obs at (" << ox << ',' << oy << "), z=" << elevs[ox * nrows + oy] 
  //      << ", vix=" << nvis << '/' << ntarget << '=' << v << endl;
  vix[ox * nrows + oy] = (unsigned char)min(255, (int)(v * 255.999f));
}

void calc_vix(const int nrows, const usint *h_elevs,
              const int roi, const int oht, const int tht,
              const int ntests, unsigned char *h_vix) {
  usint *d_elevs;
  unsigned char *d_vix;

  if (cudaMalloc((void **)&d_elevs, square(nrows) * sizeof(usint)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_vix, square(nrows) * sizeof(unsigned char)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMemcpy(d_elevs, h_elevs, square(nrows) * sizeof(usint), cudaMemcpyHostToDevice) != cudaSuccess)
    die("cudaMemcpy failed");

  const size_t dimblock = 128;
  // const size_t dimgrid = square(nrows) / dimblock + (square(nrows) % dimblock ? 1 : 0);
  int s = (int)sqrt(square(nrows) / dimblock);
  if (square(s) * dimblock < square(nrows)) s++;
  const dim3 dimgrid(s, s);
  calc_one_vix<<<dimgrid, dimblock>>>(nrows, d_elevs, roi, oht, tht, ntests, d_vix);

  if (cudaMemcpy(h_vix, d_vix, square(nrows) * sizeof(unsigned char), cudaMemcpyDeviceToHost) != cudaSuccess)
    die("cudaMemcpy failed");

  if (cudaFree(d_elevs) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_vix) != cudaSuccess) die("cudaFree failed");
}

////////////////////////////////////////////////////////////////////////////////
// findmax
////////////////////////////////////////////////////////////////////////////////

__global__ void process_one_block(const int nrows, unsigned char *vix,
                                  const float blocksize, const int nblockrows,
                                  const int nwantedperblock, int *obs) {
  // no need for if (blockIdx.x >= square(nblockrows)) return;
  extern __shared__ int results[];
  __shared__ int xmin, xmax, ymin, ymax;
  if (threadIdx.x == 0) {
    int bx = blockIdx.x / nblockrows;
    int by = blockIdx.x % nblockrows;
    xmin = (int)(blocksize * bx);
    xmax = min((int)(blocksize * (bx + 1)), nrows);
    ymin = (int)(blocksize * by);
    ymax = min((int)(blocksize * (by + 1)), nrows);
  }
  __syncthreads();

  const int width = ymax - ymin;
  const int npoints = (xmax - xmin) * (ymax - ymin);
  const int npointsperthread = npoints / blockDim.x + (npoints % blockDim.x ? 1 : 0);
  for (int i = 0; i < nwantedperblock; i++) {
    int t = threadIdx.x * npointsperthread;  // the first point, probably used
    int p1x = xmin + t / width;
    int p1y = ymin + t % width;
    unsigned char v1 = vix[p1x * nrows + p1y];
    int h1 = p1x * (p1x + p1y) * 010101010101;
    for (int j = t + 1; j < t + npointsperthread && j < npoints; j++) {
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
    results[threadIdx.x * 4] = p1x;
    results[threadIdx.x * 4 + 1] = p1y;
    results[threadIdx.x * 4 + 2] = v1;
    results[threadIdx.x * 4 + 3] = h1;
    if (threadIdx.x == 0) {
      // p1x... is the result of thread 0
      for (int j = 1; j < blockDim.x; j++) {
        int p2x = results[j * 4];
        int p2y = results[j * 4 + 1];
        int v2 = results[j * 4 + 2];
        int h2 = results[j * 4 + 3];
        if (v1 < v2 || (v1 == v2 && h1 < h2)) {
          p1x = p2x;
          p1y = p2y;
          v1 = v2;
          h1 = h2;
        }
      }
      obs[blockIdx.x * nwantedperblock * 2 + i * 2] = p1x;
      obs[blockIdx.x * nwantedperblock * 2 + i * 2 + 1] = p1y;
      vix[p1x * nrows + p1y] = 0;  // used
    }
    __syncthreads();
  }
}

void find_max(const int nrows, const unsigned char *h_vix,
              const float blocksize, const int nwanted, const int nblockrows,
              const int nwantedperblock, int *h_obs) {
  unsigned char *d_vix;
  int *d_obs;
  if (cudaMalloc((void **)&d_vix, square(nrows) * sizeof(unsigned char)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_obs, nwanted * 2 * sizeof(int)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMemcpy(d_vix, h_vix, square(nrows) * sizeof(unsigned char), cudaMemcpyHostToDevice) != cudaSuccess)
    die("cudaMemcpy failed");

  const size_t dimgrid = square(nblockrows);  // two dimensional
  const size_t dimblock = 256;
  process_one_block<<<dimgrid, dimblock, dimblock * 4 * sizeof(int)>>>(
      nrows, d_vix, blocksize, nblockrows, nwantedperblock, d_obs);

  if (cudaMemcpy(h_obs, d_obs, nwanted * 2 * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
    die("cudaMemcpy failed");

  if (cudaFree(d_vix) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_obs) != cudaSuccess) die("cudaFree failed");
}

////////////////////////////////////////////////////////////////////////////////
// viewshed
////////////////////////////////////////////////////////////////////////////////

__device__ void set_vis(const int nwpr, const int row, const int col,
                        ullint *shed) {
  atomicOr(&shed[row * nwpr + col / 64], 1ULL << (63 - col % 64));
}

__global__ void calc_one_shed(const int nrows, const usint *elevs,
                              const int roi, const int oht, const int tht,
                              const int nsheds, const int *obs,
                              ullint *sheds) {
  if (blockIdx.x >= nsheds) return;
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;
  ullint * const thisshed = &sheds[blockIdx.x * nr * nwpr];
  __shared__ int ox, oy, oz;
  if (threadIdx.x == 0) {
    ox = obs[2 * blockIdx.x];
    oy = obs[2 * blockIdx.x + 1];
    oz = elevs[ox * nrows + oy] + oht;
    set_vis(nwpr, roi, roi, thisshed);
  }
  __syncthreads();

  // Clipping xmin etc at 0, nrows-1 makes the viewshed depend on the roi, so don't.
  const int xmin = ox - roi;
  const int ymin = oy - roi;
  const int xmax = ox + roi;
  const int ymax = oy + roi;
  const int xwidth = xmax - xmin;
  const int ywidth = ymax - ymin;
  const int perimeter = 2 * (xwidth + ywidth);  // This formula is subtle
  const int ntpt = perimeter / blockDim.x + (perimeter % blockDim.x ? 1 : 0);
  int dx, dy;
  int tx, ty;
  int px, py;  // Current point
  int inciny;
  int sign;
  float slope, zslope;

  const int sector = threadIdx.x;
  for (int ip = sector * ntpt; ip < (sector + 1) * ntpt && ip < perimeter; ip++) {
    if (ip < xwidth) {
      tx = xmin + ip;
      ty = ymin;
    } else if (ip < 2 * xwidth) {
      tx = 1 + xmin - xwidth + ip;
      ty = ymax;
    } else if (ip < 2 * xwidth + ywidth) {
      tx = xmin;
      ty = 1 + ymin - 2 * xwidth + ip;
    } else {
      tx = xmax;
      ty = ymin - 2 * xwidth - ywidth + ip;
    }

    // Run a line of sight out from obs to target.
    dx = tx - ox;
    dy = ty - oy;
    inciny = abs(dx) < abs(dy);
    sign = (inciny*dy + (1-inciny)*dx) > 0 ? 1 : -1;
    slope = (float)(inciny*dx + (1-inciny)*dy) / (inciny*dy + (1-inciny)*dx);
    zslope = -99999.f;

    // i=0 would be the observer, which is always visible.
    for (int i = sign; i != (inciny ? dy : dx) + sign; i += sign) {
      int j = round(i * slope);
      px = ox + (inciny*j + (1-inciny)*i);
      py = oy + (inciny*i + (1-inciny)*j);

      // Have we reached the edge of the area?
      if (px < 0 || px >= nrows || py < 0 || py >= nrows) break;
      if (square(px - ox) + square(py - oy) > square(roi)) break;

      int pelev = elevs[px * nrows + py];
      float s = (float)(pelev - oz) / abs(i);
      if (zslope < s) zslope = s;
      float hz = oz + zslope * abs(i);
      if (pelev + tht >= hz)
        set_vis(nwpr, px - ox + roi, py - oy + roi, thisshed);
    }
  }
}

void calc_sheds(const int nrows, const usint *h_elevs,
                const int roi, const int oht, const int tht,
                const int nsheds, const int *h_obs,
                ullint *h_sheds) {
  usint *d_elevs;
  int *d_obs;
  ullint *d_sheds;
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;

  if (cudaMalloc((void **)&d_elevs, square(nrows) * sizeof(usint)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_obs, 2 * nsheds * sizeof(int)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_sheds, nsheds * nr * nwpr * sizeof(ullint)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMemcpy(d_elevs, h_elevs, square(nrows) * sizeof(usint), cudaMemcpyHostToDevice) != cudaSuccess)
    die("cudaMemcpy failed");
  if (cudaMemcpy(d_obs, h_obs, 2 * nsheds * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess)
    die("cudaMemcpy failed");
  if (cudaMemset(d_sheds, 0, nsheds * nr * nwpr * sizeof(ullint)) != cudaSuccess)
    die("cudaMemset failed");

  const size_t dimblock = 256;  // must be a multiple of 32
  const size_t dimgrid = nsheds;
  calc_one_shed<<<dimgrid, dimblock>>>(nrows, d_elevs, roi, oht, tht, nsheds, d_obs, d_sheds);

  if (cudaMemcpy(h_sheds, d_sheds, nsheds * nr * nwpr * sizeof(ullint), cudaMemcpyDeviceToHost)
      != cudaSuccess)
    die("cudaMemcpy failed");

  if (cudaFree(d_elevs) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_obs) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_sheds) != cudaSuccess) die("cudaFree failed");
}

////////////////////////////////////////////////////////////////////////////////
// site
////////////////////////////////////////////////////////////////////////////////

__host__ int is_obs_vis(const int cumnwpr,
                        const ullint *cumshed,
                        const int *observer) {
  int i = observer[0] * cumnwpr * 64 + observer[1];
  if ((cumshed[i / 64] & 1ULL << (63 - i % 64)) != 0)
    return 1;
  else
    return 0;
}

__device__ int is_obs_vis(const int cumnwpr,
                          const ullint *cumshed,
                          const int obsx, const int obsy) {
  int i = obsx * cumnwpr * 64 + obsy;
  if ((cumshed[i / 64] & 1ULL << (63 - i % 64)) != 0)
    return 1;
  else
    return 0;
}

__global__ void calc_extra_area(const int nrows, const int roi,
                                const int *obs, const int *updatelist,
                                const ullint *sheds,
                                const int nr, const int nwpr, const int cumnwpr,
                                const ullint *cumshed,
                                int *testshedarea) {
  const int oi = updatelist[blockIdx.x];
  const int obsx = obs[oi*2];
  const int obsy = obs[oi*2+1];
  const ullint *shed = sheds + oi*nr*nwpr;
  int *extraarea = testshedarea + oi;

  extern __shared__ int areas[];
  areas[threadIdx.x] = 0;

  // calculate nrpt rows of extra area
  const int nrpt = nr / blockDim.x + (nr % blockDim.x ? 1 : 0);
  int sum = 0;

  for (int row = threadIdx.x * nrpt; row < (threadIdx.x + 1) * nrpt && row < nr; row++) {
    const int cumrow = obsx - roi + row;
    if (cumrow >= 0 && cumrow < nrows) {
      int firstword = (obsy - roi) / 64;
      int firstbit = (obsy - roi) % 64;
      if (firstbit < 0) {
        firstword--;
        firstbit += 64;
      }
      int lastword = (obsy + roi) / 64;

      ullint prevvalue = 0ULL;
      ullint value, cumvalue, tempvalue;

      for (int cumword = firstword; cumword <= lastword; cumword++)
        if (cumword >= 0 && cumword < cumnwpr) {
          int word = cumword - firstword;  // definition out of loop?
          if (cumword == 0 && word > 0) prevvalue = shed[row * nwpr + word - 1];
          if (word < nwpr)
            value = shed[row * nwpr + word];
          else
            value = 0ULL;
          cumvalue = cumshed[cumrow * cumnwpr + cumword];
          tempvalue = cumvalue;
          tempvalue |= value >> firstbit;
          if (firstbit != 0) tempvalue |= prevvalue << (64 - firstbit);
          tempvalue ^= cumvalue;
          sum += __popcll(tempvalue);
          prevvalue = value;
        }
    }
  }

  areas[threadIdx.x] = sum;
  __syncthreads();  // wait for all threads

  if (threadIdx.x == 0) {
    int sum = 0;
    for (int i = 0; i < blockDim.x; i++)
      sum += areas[i];
    *extraarea = sum;
  }
}

__global__ void union_area(const int nrows, const int roi, const int intervis,
                           const int nsheds, const int *obs,
                           const ullint *sheds,
                           const int nr, const int nwpr, const int cumnwpr,
                           const int nusedsheds, const char *usedq,
                           const int lastobs,
                           const int lastobsx, const int lastobsy,
                           const ullint *cumshed,
                           int *testshedarea, int *updatelist) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int obsx, obsy, valid;

  // calculate multiple extra areas
  for (int oi = tid; oi < nsheds; oi += gridDim.x * blockDim.x) {
    valid = 1;
    if (usedq[oi]) {
      valid = 0;
    } else {
      obsx = obs[oi*2];
      obsy = obs[oi*2 + 1];
      if (nusedsheds > 0) {  // lastobs >= 0
        if (square(obsx-lastobsx) + square(obsy-lastobsy) > square(2*roi))
          valid = 0;
        else if (intervis && !is_obs_vis(cumnwpr, cumshed, obsx, obsy))  // if intervis, reset all after the first
          valid = 0;  // testshedarea[oi] = 0;  // reset invisible ones to zero
      }
    }

    if (valid) {
      /*
      cudaStream_t s;
      cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
      calc_extra_area<<<1, nr, nr*sizeof(int), s>>>(nrows, roi, obsx, obsy,
                                                    &sheds[oi*nr*nwpr],
                                                    nr, nwpr, cumnwpr,
                                                    cumshed, &testshedarea[oi]);
      cudaStreamDestroy(s);
      // calc_extra_area<<<1, nr, nr*sizeof(int)>>>(nrows, roi, obsx, obsy, &sheds[oi*nr*nwpr],
      //                                            nr, nwpr, cumnwpr, cumshed, &testshedarea[oi]);
      */
      int index = atomicAdd(updatelist+499999, 1);
      updatelist[index] = oi;
    }
  }
}

__device__ void swap(int *x, int *y) {
  int z = *x;
  *x = *y;
  *y = z;
}

__global__ void findtopobs(const int nsheds, const char *usedq,
                           const int *testshedarea, int *top100) {
  __shared__ int obs[256];
  __shared__ int areas[256];
  obs[threadIdx.x] = 0;
  areas[threadIdx.x] = 0;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int size = gridDim.x * blockDim.x;
  for (int i = tid; i < nsheds; i += size)
    if (!usedq[i]) {
      int extraarea = testshedarea[i];
      if (extraarea > areas[threadIdx.x]) {
        obs[threadIdx.x] = i;
        areas[threadIdx.x] = extraarea;
      }
    }
  __syncthreads();

  if (threadIdx.x < 128 && areas[threadIdx.x + 128] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 128]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 128]);
  }
  if (threadIdx.x < 64 && areas[threadIdx.x + 64] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 64]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 64]);
  }
  if (threadIdx.x < 32 && areas[threadIdx.x + 32] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 32]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 32]);
  }
  if (threadIdx.x < 16 && areas[threadIdx.x + 16] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 16]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 16]);
  }
  if (threadIdx.x < 8 && areas[threadIdx.x + 8] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 8]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 8]);
  }
  if (threadIdx.x < 4 && areas[threadIdx.x + 4] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 4]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 4]);
  }
  if (threadIdx.x < 2 && areas[threadIdx.x + 2] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 2]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 2]);
  }
  if (threadIdx.x < 1 && areas[threadIdx.x + 1] > areas[threadIdx.x]) {
    swap(&obs[threadIdx.x], &obs[threadIdx.x + 1]);
    swap(&areas[threadIdx.x], &areas[threadIdx.x + 1]);
  }

  if (threadIdx.x == 0) {
    top100[blockIdx.x] = obs[0];
    top100[blockIdx.x + 100] = areas[0];
  }
}

__global__ void calc_union(const int nrows, const int roi,
                           const int nsheds, const int *obs,
                           const ullint *sheds,
                           const int nr, const int nwpr, const int cumnwpr,
                           const int lastobs, char *usedq,
                           ullint *cumshed) {
  // lastobs >= 0
  __shared__ int lastobsx, lastobsy;
  __shared__ const ullint *shed;
  if (threadIdx.x == 0) {
    usedq[lastobs] = 1;  // set gridDim.x times
    lastobsx = obs[lastobs*2];
    lastobsy = obs[lastobs*2 + 1];
    shed = &sheds[lastobs*nr*nwpr];
  }
  __syncthreads();  // wait for thread 0

  int firstword = (lastobsy - roi) / 64;
  int firstbit = (lastobsy - roi) % 64;
  if (firstbit < 0) {
    firstword--;
    firstbit += 64;
  }

  int row = blockIdx.x*blockDim.x + threadIdx.x;
  if (row < nr) {  // a row of shed
    int cumrow = lastobsx - roi + row;
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

void site_it(const int nrows, const int roi, const int intervis,
             const int nsheds, const int *h_obs, const ullint *h_sheds, char *selected) {
  const int nr = 2 * roi + 1;
  const int nwpr = (nr + 63) / 64;
  const int cumnwpr = (nrows + 63) / 64;

  int *usedsheds;     // list of sheds used so far.
  char *h_usedq;      // whether each particular shed has been used.
  int h_top100[200];  // top 100 tentative observers and extra areas

  usedsheds = 0;
  h_usedq = 0;
  int *areas = 0;
  usedsheds = new int[nsheds];
  h_usedq = new char[nsheds];
  areas = new int[nsheds];
  if (!usedsheds || !h_usedq || !areas)
    die("Memory exhausted. Program terminates.");
  for (int i = 0; i < nsheds; i++) h_usedq[i] = 0;

  int *d_obs;
  ullint *d_sheds;
  char *d_usedq;
  ullint *d_cumshed;
  int *d_testshedarea;
  int *d_top100;
  if (cudaMalloc((void **)&d_obs, 2 * nsheds * sizeof(int)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_sheds, nsheds * nr * nwpr * sizeof(ullint)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_usedq, nsheds * sizeof(char)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_cumshed, nrows * cumnwpr * sizeof(ullint)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_testshedarea, nsheds * sizeof(int)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMalloc((void **)&d_top100, 200 * sizeof(int)) != cudaSuccess)
    die("cudaMalloc failed");
  if (cudaMemcpy(d_obs, h_obs, 2 * nsheds * sizeof(int), cudaMemcpyHostToDevice) != cudaSuccess)
    die("cudaMemcpy failed");
  if (cudaMemcpy(d_sheds, h_sheds, nsheds * nr * nwpr * sizeof(ullint), cudaMemcpyHostToDevice)
      != cudaSuccess)
    die("cudaMemcpy failed");
  if (cudaMemset(d_usedq, 0, nsheds * sizeof(char)) != cudaSuccess)
    die("cudaMemset failed");
  if (cudaMemset(d_cumshed, 0, nrows*cumnwpr*sizeof(ullint)) != cudaSuccess)
    die("cudaMemset failed");

  int *d_updatelist;
  if (cudaMalloc((void **)&d_updatelist, 500000*sizeof(int)) != cudaSuccess) die("cudaMalloc failed");

  int nusedsheds = 0;
  int lastobs = -1;
  int lastobsx = 0;
  int lastobsy = 0;
  int cumarea = 0;
  if (cudaMemset(d_testshedarea, 0, nsheds * sizeof(int)) != cudaSuccess)
    die("cudaMemset failed");
  size_t dimgrid = (nsheds+255)/256;
  size_t dimblock = 256;
  // size_t dimgrid = (nsheds + dimblock - 1) / dimblock;

  //cout << "Total area=" << square(nrows) << endl;
  //cout << "#nusedsheds newshed obsx obsy area extraarea newcumarea areapercentage" << endl;

  while (1) {
    if (cudaMemset(d_updatelist, 0, 500000*sizeof(int)) != cudaSuccess) die("cudaMemset failed");

    union_area<<<dimgrid, dimblock>>>(
        nrows, roi, intervis, nsheds, d_obs, d_sheds, nr, nwpr, cumnwpr,
        nusedsheds, d_usedq, lastobs, lastobsx, lastobsy, d_cumshed, d_testshedarea, d_updatelist);

    int updatelistsize;
    if (cudaMemcpy(&updatelistsize, d_updatelist+499999, sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess) die("cudaMemcpy failed");
      //cout << "\nsize = " << updatelistsize;

    calc_extra_area<<<updatelistsize, nr, nr*sizeof(int)>>>(
        nrows, roi, d_obs, d_updatelist, d_sheds, nr, nwpr, cumnwpr, d_cumshed, d_testshedarea);

    cudaDeviceSynchronize();

    // find top 100 observers
    findtopobs<<<100, 256>>>(nsheds, d_usedq, d_testshedarea, d_top100);
    if (cudaMemcpy(h_top100, d_top100, 200 * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
      die("cudaMemcpy failed");

    int newshed = 0;
    int extraarea = 0;
    for (int i = 0; i < 100; i++) {
      if (h_top100[i + 100] > extraarea) {
        extraarea = h_top100[i + 100];
        newshed = h_top100[i];
      }
    }
    if (extraarea == 0) {
      //cout << "No more new observers that will add new area." << endl;
      break;
    }

    usedsheds[nusedsheds++] = newshed;
    h_usedq[newshed] = 1;
    lastobs = newshed;
    lastobsx = h_obs[newshed * 2];
    lastobsy = h_obs[newshed * 2 + 1];

    // set usedq and calculate cumshed
    calc_union<<<nr, 1>>>(nrows, roi, nsheds, d_obs, d_sheds, nr, nwpr, cumnwpr, lastobs, d_usedq, d_cumshed);
    cudaDeviceSynchronize();

    cumarea += extraarea;
    double areapercentage = 100.0 * cumarea / square(nrows);
    if (areapercentage > 95)
        break;

    if (nusedsheds == 1) {
      if (cudaMemcpy(areas, d_testshedarea, nsheds * sizeof(int), cudaMemcpyDeviceToHost) != cudaSuccess)
        die("cudaMemcpy failed");
      if (intervis)
        if (cudaMemset(d_testshedarea, 0, nsheds * sizeof(int)) != cudaSuccess)
          die("cudaMemset failed");
    }

    /*
    cout << setw(6) << nusedsheds << setw(6) << newshed
         << setw(6) << h_obs[newshed * 2] << setw(6) << h_obs[newshed * 2 + 1]
         << setw(8) << areas[newshed] << setw(8) << extraarea
         << setw(10) << cumarea << setw(8) << areapercentage << endl;
    */
  }
  cout << " nusedsheds:" << nusedsheds << " coverage:" << 100.0*cumarea/square(nrows);

  if (cudaFree(d_obs) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_sheds) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_usedq) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_cumshed) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_testshedarea) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_top100) != cudaSuccess) die("cudaFree failed");
  if (cudaFree(d_updatelist) != cudaSuccess) die("cudaFree failed");

  for (int i = 0; i < nsheds; i++)
    selected[i] = h_usedq[i];
  delete[] usedsheds;
  delete[] h_usedq;
  delete[] areas;
}

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
  // clock_t begin, middle, end;
  double elapsed_secs;
  // begin = middle = clock();
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
  unsigned char *vix;             // visibility index * 256 (output); vim
  int *obs;                       // observers
  ullint *sheds;  // viewsheds
  char *selected;

  //cerr << "[SITE, compiled from " << __FILE__ << " on " << __DATE__ << ", " << __TIME__ << ']' << endl;
  if (argc != 10) {
    cerr << "argc=" << argc << endl;
    die("SITE requires 9 arguments: nrows, roi, oht/tht, ntests, blocksize, nwanted, intervis, infile, outfile");
  }

  nrows = atoi(argv[1]);
  roi = atoi(argv[2]);
  oht = tht = atoi(argv[3]);
  ntests = atoi(argv[4]);
  blocksize0 = atoi(argv[5]);
  nwanted0 = atoi(argv[6]);
  intervis = atoi(argv[7]);

  //cerr << "nrows=" << nrows << ", roi=" << roi << ", oht=" << oht << ", tht=" << tht
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

  //cerr << "blocksize=" << blocksize << ", nblocks=" << nblocks
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
  if (cudaMallocHost((void **)&elevs, square(nrows) * sizeof(usint)) != cudaSuccess)
    die("cudaMallocHost failed");
  if (cudaMallocHost((void **)&vix, square(nrows) * sizeof(unsigned char)) != cudaSuccess)
    die("cudaMallocHost failed");
  if (cudaMallocHost((void **)&obs, 2 * nwanted * sizeof(int)) != cudaSuccess)
    die("cudaMallocHost failed");
  if (cudaMallocHost((void **)&sheds, nwanted * nr * nwpr * sizeof(ullint)) != cudaSuccess)
    die("cudaMallocHost failed");
*/
  /*
  for (int i = 0; i < nrows; i++)
    for (int j = 0; j < nrows; j++) {
      cin.read(reinterpret_cast<char *>(&elevs[i * nrows + j]), sizeof(usint));
      if (cin.fail()) {
        cerr << "Error: at i=" << i << ", j=" << j << endl;
        die("Input failed");
      }
    }
  */
  ifstream ifs(argv[8]);
  ifs.read((char *)elevs, square(nrows)*sizeof(usint));
  if (ifs.fail())
    die("Input failed");
  ifs.close();

  // end = clock();
  // elapsed_secs = float(end - middle) / CLOCKS_PER_SEC;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cerr << "input:" << elapsed_secs;

  calc_vix(nrows, elevs, roi, oht, tht, ntests, vix);
  cudaDeviceSynchronize();
  // end = clock();
  // elapsed_secs = float(end - middle) / CLOCKS_PER_SEC;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cerr << " vix:" << elapsed_secs;
  /*
  int x;
  ullint sum = 0;
  ifstream ifs("vim.bin");
  for (int i = 0; i < square(nrows); i++) {
    ifs >> x;
    sum += square(vix[i]-x);
  }
  cerr << " RMS VIM error:" << sqrt((double)sum/square(nrows));
  ifs.close();
  */

  find_max(nrows, vix, blocksize, nwanted, nblocks, nwantedperblock, obs);
  cudaDeviceSynchronize();
  // end = clock();
  // elapsed_secs = float(end - middle) / CLOCKS_PER_SEC;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cerr << " findmax:" << elapsed_secs;
  
  calc_sheds(nrows, elevs, roi, oht, tht, nwanted, obs, sheds);
  cudaDeviceSynchronize();
  // end = clock();
  // elapsed_secs = float(end - middle) / CLOCKS_PER_SEC;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cerr << " viewshed:" << elapsed_secs;

  site_it(nrows, roi, intervis, nwanted, obs, sheds, selected);
  cudaDeviceSynchronize();
  // end = clock();
  // elapsed_secs = float(end - middle) / CLOCKS_PER_SEC;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cerr << " site:" << elapsed_secs;

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
/*
  if (cudaFreeHost(elevs) != cudaSuccess) die("cudaFreeHost failed");
  if (cudaFreeHost(vix) != cudaSuccess) die("cudaFreeHost failed");
  if (cudaFreeHost(obs) != cudaSuccess) die("cudaFreeHost failed");
  if (cudaFreeHost(sheds) != cudaSuccess) die("cudaFreeHost failed");
*/
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(middle, end);
  middle = end;
  cerr << " output:" << elapsed_secs;
  // end = clock();
  // elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  gettimeofday(&end, NULL);
  elapsed_secs = calc_time(begin, end);
  cerr << " total:" << elapsed_secs << endl;
}
