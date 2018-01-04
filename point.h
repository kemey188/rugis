#ifndef _H_RUGIS_POINT_H_
#define _H_RUGIS_POINT_H_

#include <stdint.h>
#include <algorithm>

template <typename T>
struct Point2D {
  T x, y;
  Point2D():x(0), y(0) {}
  Point2D(T x, T y): x(x), y(y) {}
};

template <typename T>
struct Point3D {
  T x, y, z;
  Point3D(): x(0), y(0), z(0){}
  Point3D(T x, T y, T z): x(x), y(y), z(z){}
};

typedef Point2D<uint32_t> XYPoint;
typedef Point2D<double> LatLngPoint;

template <typename T>
inline T Clip(T max, T min, T v) {
    return std::min<T>(max, std::max<int>(min, v));
}

#endif
