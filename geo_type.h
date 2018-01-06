#ifndef _H_GEO_TYPE_H_
#define _H_GEO_TYPE_H_

#include <string>
#include "google/s2/s2loop.h"

enum AOIType {
  kUnknown = 0,
  kAdmin = 1, // adcode
  kBArea = 2, // aoi
  kTArea = 3, // business_area
  kBuild = 4  // building
};

struct AOI {
  int id;
  int adcode; // adcode
  int clazz;  // category
  int area;   // area
  int circum; // circum
  int latE6;  // lat,lng of center
  int lngE6;  // 
  std::string name;
  S2Loop   *bounds;

  AOI():
    id(0), adcode(0), clazz(0), area(0), circum(0),
    latE6(0), lngE6(0),
    bounds(NULL) { }

  AOI(int id):
    id(id), adcode(0), clazz(0), area(0), circum(0),
    latE6(0), lngE6(0),
    bounds(NULL) { }

  AOI(int id, int adcode, int clazz, int area, int circum, int latE6, int lngE6, std::string name):
    id(id), adcode(adcode), clazz(clazz), area(area), circum(circum),
    latE6(latE6), lngE6(lngE6), name(name),
    bounds(NULL) { }

  AOI(const AOI &other) {
    id     = other.id;
    adcode = other.adcode;
    clazz  = other.clazz;
    area   = other.area;
    circum = other.circum;
    latE6  = other.latE6;
    lngE6  = other.lngE6;
    name   = other.name;
    if (other.bounds) 
      bounds = other.bounds->Clone();
    else
      bounds = NULL;
  }

  ~AOI() {
    if(bounds) delete bounds;
    bounds = NULL;
  }

  bool operator==(const AOI &other) const {
    return id < other.id;
  }

  bool operator<(const AOI &other) const {
    return id < other.id;
  }
};

struct POI {
  int id;
  int adcode;
  int clazz;
  int latE6;
  int lngE6;
  POI():
    id(0), adcode(0), clazz(0), latE6(0), lngE6(0) {
    }
};

#endif
