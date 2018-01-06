#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <boost/algorithm/string/split.hpp>

#include "google/s2/s2loop.h"
#include "google/util/coding/coder.h"
#include "aoi_index.h"

inline bool FileExists(const std::string fpath) {
  struct stat st;
  return 0 == stat(fpath.c_str(), &st);
}

#define CheckFile(path) do { \
  if(!FileExists(path)) { \
    std::cerr<<path<<" not exists\n"; \
    return false; \
  } \
} while(0)

bool AreaIndex::Load(const char *prefix) {

  std::string index_path(prefix);
  std::string bound_path(prefix);

  index_path.append(".idx");
  bound_path.append(".dat");

  CheckFile(index_path);
  CheckFile(bound_path);
  
  // build index
  int id, adcode, clazz, area, circum, latE6, lngE6;
  std::string name;
  double lat_lo, lng_lo, lat_hi, lng_hi;

  std::vector<value_t> boxes;

  std::ifstream ifs(index_path.c_str());

  if (!ifs) return false;

  std::string line;
  std::vector<std::string> parts; 

  while(std::getline(ifs, line)) {

    parts.clear();
    boost::split(parts, line, boost::is_any_of("\t,"));

    if (parts.size() != 12) {
      std::cerr<<"bad line "<<line;
      return false;
    }

    id     = atoi(parts[0].c_str());
    adcode = atoi(parts[1].c_str());
    clazz  = atoi(parts[2].c_str());
    area   = atoi(parts[3].c_str());
    circum = atoi(parts[4].c_str());
    latE6  = atoi(parts[5].c_str());
    lngE6  = atoi(parts[6].c_str());
    lat_lo = strtod(parts[7].c_str(),  NULL);
    lng_lo = strtod(parts[8].c_str(),  NULL);
    lat_hi = strtod(parts[9].c_str(),  NULL);
    lng_hi = strtod(parts[10].c_str(), NULL);
    name   = parts[11];

    boxes.push_back(std::make_pair(
          box_t(
            point_t(lat_lo, lng_lo),
            point_t(lat_hi, lng_hi)),
          id));
    entries.push_back(AOI(id, adcode, clazz, area, circum, latE6, lngE6, name));
  }
  ifs.close();

  index.insert(boxes.begin(), boxes.end());
  std::sort(entries.begin(), entries.end());

  // read coordinates
  ifs.open(bound_path.c_str(), std::ios::binary);
  if (!ifs) return false;

  uint64_t hdr = 0;
  std::string coords;
  Decoder coder;
  while(ifs.read((char *)&hdr, 8)) {
    int len = hdr & 0xffffffff;
    int rid = hdr >> 32;
    coords.resize(len);
    ifs.read((char*)(coords.c_str()), len);
    AOI *p = GetEntry(rid);
    if (p) {
      p->bounds = new S2Loop();
      coder.reset(coords.c_str(), coords.length());
      p->bounds->Decode(&coder);
    }
  }
  ifs.close();

  return true;
}

void AreaIndex::Find(double lat, double lng, std::vector<int>& res) const {
  typedef boost::geometry::model::d2::point_xy<double> point_t;

  res.clear();

  std::vector<value_t> cands;

  point_t pt(lat, lng);
  index.query(boost::geometry::index::contains(pt), std::back_inserter(cands));

  if (cands.empty()) {
    double eps = 1e-5;
    box_t target(point_t(lat-eps, lng-eps), point_t(lat+eps, lng+eps));
    index.query(boost::geometry::index::intersects(target), std::back_inserter(cands));
  }

  // std::cout<<"find "<<cands.size()<<" candidates\n";
  S2Point point = S2LatLng::FromDegrees(lat, lng).ToPoint();
  for(const value_t &c: cands) {
    const AOI *aoi = GetEntry(c.second);
    if (aoi && aoi->bounds && aoi->bounds->Contains(point))
      res.push_back(c.second);
  }

  if (res.empty()) {
    S2Point point = S2LatLng::FromDegrees(lat+1e-4, lng+1e-4).ToPoint();
    for(const value_t &c: cands) {
      const AOI *aoi = GetEntry(c.second);
      if (aoi && aoi->bounds && aoi->bounds->Contains(point))
        res.push_back(c.second);
    }
  }
}

const double kMaxArea = 4*M_PI;

int AreaIndex::Find(double lat, double lng) const {
  std::vector<int> res;
  Find(lat, lng, res);
  if (res.empty()) return 0;

  size_t i,j;
  int min_area = 1000000000;
  for(i = 0, j = 0; i < res.size(); ++ i) {
    const AOI *aoi = GetEntry(res[i]);
    if (aoi->area < min_area) {
      min_area = aoi->area;
      j = i;
    }
  }
  return res[j];
}
