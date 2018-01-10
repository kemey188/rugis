#ifndef _H_RUGIS_COMMON_COORDSYS_H_
#define _H_RUGIS_COMMON_COORDSYS_H_

#include <math.h>
#include "point.h"

const int kGroundLevel = 19;
const int kMaxMapSize = (1<<(kGroundLevel+8)) - 1;
const double kMinLatitude = -85.05112878;
const double kMaxLatitude = 85.05112878;
const double kMinLngitude = -180;
const double kMaxLngitude = 180;

// 以瓦片形式表示的地图上的一个点
// 最大精度为kGroundLevel
// X,Y = tile[level-bit]|inner[8-bit]|padding[19-level-bit] 
struct MapPoint: XYPoint {

  MapPoint(int tx, int ty, int x, int y):
    XYPoint(((tx<<8) + x),
        ((ty<<8) + y)) { }

  MapPoint(int level, int tx, int ty, int x, int y):
    XYPoint(((tx<<8) + x)<<(kGroundLevel - level),
        ((ty<<8) + y)<<(kGroundLevel - level)) { }

  MapPoint(double lat, double lng) {
    lat = Clip(lat, kMinLatitude, kMaxLatitude);
    lng = Clip(lng, kMinLngitude, kMaxLngitude);
    double x = (lng + 180) / 360;
    double sinLatitude = sin(lat * M_PI / 180);
    double y = 0.5 - log((1 + sinLatitude) / (1 - sinLatitude)) / (4 * M_PI);
    size_t mapSize = 256UL<<kGroundLevel;
    this->x = Clip<int>(x * mapSize + 0.5, 0, mapSize - 1);
    this->y = kMaxMapSize - Clip<int>(y * mapSize + 0.5, 0, mapSize - 1);
  }

  // 取瓦片坐标
  XYPoint tile(int level) {
    if (level >= kGroundLevel)
      level = kGroundLevel;
    int delta = level-kGroundLevel;
    return XYPoint(x>>(8+delta), y>>(8+delta));
  }

  // 瓦片内的相对坐标
  XYPoint inner(int level) {
    if (level >= kGroundLevel)
      level = kGroundLevel;
    int delta = level-kGroundLevel;
    return XYPoint((x>>delta)&0xff, (y>>delta)&0xff);
  }
};

// 球面坐标: WGS84, GCJ02, BD09
// 平面坐标: 瓦片坐标

#endif
