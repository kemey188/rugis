#ifndef _H_RUGIS_AOI_INDEX_H_
#define _H_RUGIS_AOI_INDEX_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <map>
#include "geo_type.h"

struct AreaIndex {
  typedef boost::geometry::model::d2::point_xy<double> point_t;
  typedef boost::geometry::model::box<boost::geometry::model::d2::point_xy<double>>  box_t;
  typedef std::pair<box_t, int> value_t;
  typedef boost::geometry::index::rtree< value_t, boost::geometry::index::rstar<16,4> > rtree_t;

  bool Load(const char *prefix);

  // 返回包含Point(lat, lng)的区域ID
  // 如果有多个点,返回面积最小的区域
  int  Find(double lat, double lng) const;
  void Find(double lat, double lng, std::vector<int>& res) const;

  inline std::string Name(int id) const {
    const AOI *aoi = GetEntry(id);
    return aoi ? aoi->name : std::string();
  }

  const AOI *GetEntry(int id) const {
    std::vector<AOI>::const_iterator it = std::lower_bound(
        entries.begin(), entries.end(), AOI(id));
    if (it != entries.end() && it->id == id) {
      return &(*it);
    } else return NULL;
  }

  AOI *GetEntry(int id) {
    std::vector<AOI>::iterator it = std::lower_bound(
        entries.begin(), entries.end(), AOI(id));
    if (it != entries.end() && it->id == id) {
      return &(*it);
    } else return NULL;
  }

  rtree_t                     index;
  std::vector<AOI>            entries;

};

#endif
