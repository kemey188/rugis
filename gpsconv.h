#ifndef __H_RUGIS_GPSCONVERTER_H__
#define __H_RUGIS_GPSCONVERTER_H__

// wgs84 -> gcj02
int GPS2GCJ(double wg_lat, double wg_lng, double *gcj_lat, double *gcj_lng);
// wgs84 -> mars
bool wgs84_to_mars_simple(double wgLat, double wgLon, double &mgLat, double &mgLon);
// wgs84 -> mars
bool wgs84_to_mars(const double wgLat, const double wgLon, double& mgLat, double& mgLon, int altitude = 0);
// bd <-> gcj02
void bd_encrypt(double gg_lat, double gg_lon, double &bd_lat, double &bd_lon);
void bd_decrypt(double bd_lat, double bd_lon, double &gg_lat, double &gg_lon);

#endif
