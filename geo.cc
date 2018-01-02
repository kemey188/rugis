#include <algorithm>
#include <cmath>
#include "s2/s2cellid.h"
#include "s2/s2latlng.h"
#include "geo.h"
#include "coordtrans.h"

namespace geo {
  static const double kEarthRadius = 6378137;
  static const double kEarthCircumferenceMeters = 2 * M_PI * kEarthRadius; //40075017;

  inline double EarthMetersToRadians(double meters) {
    return meters / kEarthRadius;
  }

  inline double RadiansToEarthMeters(double radians) {
    return radians * kEarthRadius;
  }

  double Distance(double latitude1, double lngitude1, double latitude2, double lngitude2) {
    double dx = sin((latitude1 - latitude2) * M_PI / 360);
    double dy = sin((lngitude1 - lngitude2) * M_PI / 360);
    return 2 * kEarthRadius * asin(sqrt(dx*dx + cos(latitude1*M_PI/180)*cos(latitude2*M_PI/180)*dy*dy));
  }

  int DistanceOfId(uint64_t c1, uint64_t c2) {
    return RadiansToEarthMeters(
        S2CellId(c1).ToLatLng().GetDistance(
          S2CellId(c2).ToLatLng()).radians());
  }

  u64 LatLng2Id(double lat, double lng, int level) {
    return S2CellId::FromLatLng(S2LatLng::FromDegrees(lat, lng)).parent(level).id();
  }

  u64 LatLng2Id(i32 lat, i32 lng, int level) {
    return S2CellId::FromLatLng(S2LatLng::FromE6(lat, lng)).parent(level).id();
  }

  void Id2LatLng(u64 id, float *lat, float *lng, int *level) {
    S2CellId sid(id);
    S2LatLng ll = sid.ToLatLng();
    if (lat) *lat = ll.lat().degrees();
    if (lng) *lng = ll.lng().degrees();
    if (level) *level = sid.level();
  }

  // BMap.Projection.convertLL2MC(new BMap.Point(113.768261,23.036282))
  // => 12664762.69, 2619536

  #define CHECK_LATLNG(lat, lng) \
  if (lat >= 90 || lat <= -90 || lng >= 180 || lng <= -180) return false

  bool LatLng2UTM(double lat, double lng, double &x, double &y) {
    CHECK_LATLNG(lat, lng);
    int band  = int(std::floor(fabs(lat)));
    if (band >= 75) band = 0;
    else band = 5 - band/15;

    const double *C = LL2MC[band];
    const double alpha = fabs(lat) / C[9];

    x = C[0] + C[1] * fabs(lng);
    y = C[2] + C[3] * alpha
      + C[4] * std::pow(alpha, 2)
      + C[5] * std::pow(alpha, 3)
      + C[6] * std::pow(alpha, 4)
      + C[7] * std::pow(alpha, 5)
      + C[8] * std::pow(alpha, 6);

    if (lng < 0) x *= -1;
    if (lat < 0) y *= -1;

    return true;
  }

  bool UTM2LatLng(double x, double y, double &lat, double &lng) {
    int band = 0;
    for(int i = 0; i < 6; ++ i) {
      if (fabs(y) >= MCBAND[i]) { band = i; break; }
    }
    const double *C = MC2LL[band];
    const double alpha = fabs(y) / C[9];

    lng = C[0] + C[1] * fabs(x);
    lat = C[2] + C[3] * alpha
      + C[4] * std::pow(alpha, 2)
      + C[5] * std::pow(alpha, 3)
      + C[6] * std::pow(alpha, 4)
      + C[7] * std::pow(alpha, 5)
      + C[8] * std::pow(alpha, 6);
    if (x < 0) lng *= -1;
    if (y < 0) lat *= -1;
    return true;
  }

  int GPS2GCJ(double wg_lat, double wg_lng, double *gcj_lat, double *gcj_lng) {
	  *gcj_lat = 0;
	  *gcj_lng = 0;
	  double g_lat = 0;
	  double g_lng = 0;
	  int r_val = -1;
	  r_val = GPS_to_GCJ(wg_lng, wg_lat, g_lng , g_lat);

	  if(r_val==0) {
		  r_val = 1;
		  *gcj_lat = g_lat;
		  *gcj_lng = g_lng;
	  } else {
		  r_val = 0;
	  }

	  return r_val;
  }

  int GCJ2GPS(double gcj_lat, double gcj_lng, double *wg_lat, double *wg_lng) {
	  *wg_lat = 0;
	  *wg_lng = 0;
	  double wg84_lat = 0;
	  double wg84_lng = 0;
	  int r_val = -1;
	  r_val = GCJ_to_GPS(gcj_lng, gcj_lat, wg84_lng, wg84_lat);
	  if(r_val==0) {
		  r_val = 1;
		  *wg_lat = wg84_lat;
		  *wg_lng = wg84_lng;
	  } else {
		  r_val = 0;
	  }
	  return r_val; 
  }

};
