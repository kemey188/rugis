#ifndef _H_RUGIS_GEOMETRY_H_
#define _H_RUGIS_GEOMETRY_H_

namespace geo {
  typedef unsigned long u64;
  typedef unsigned int u32;
  typedef int i32;

  double Distance(double lat1, double lng1, double lat2, double lng2);
  u64 LatLng2Id(double lat, double lng, int level);
  u64 LatLng2Id(i32 lat, i32 lng, int level);
  void Id2LatLng(u64 id, float *lat, float *lng, int *level);
  int GPS2GCJ(double wg_lat, double wg_lng, double *gcj_lat, double *gcj_lng);
  int DistanceOfId(u64 c1, u64 c2);
  bool LatLng2UTM(double lat, double lng, double &x, double &y);
  bool UTM2LatLng(double x, double y, double &lat, double &lng);
};

#endif
