#ifndef _H_RUGIS_DBSCAN_FILTER_H__
#define _H_RUGIS_DBSCAN_FILTER_H__

#include <vector>
#include "google/s2/s2cellid.h"
#include "google/s2/s2latlng.h"
#include "google/s2/s2cellunion.h"

#ifdef __cplusplus
extern "C" {
#endif

int dbscan(double *lat, double *lng, short *cluster, int n, double dc, int minPts);
int dpfilt(double *lat, double *lng, int n, double threshold, double dc);
int dpfilt_slow(double *lat, double *lng, int n, double threshold, double dc);


#ifdef __cplusplus
}
#endif


inline static S2CellId LatLng2CellId(double lat, double lng) {
    return S2CellId::FromLatLng(S2LatLng::FromDegrees(lat, lng));
}

struct DPoint {
    S2CellId id;
    double lat;
    double lng;
    int  samples;
    bool is_core;
    bool is_noise;
    int  cluster;
    int  density;

    DPoint(double lat, double lng):
        id(LatLng2CellId(lat, lng)),
        lat(lat),
        lng(lng),
        samples(1),
        is_core(false),
        is_noise(false),
        cluster(-1),
        density(0)
    {
    }

    DPoint(double lat, double lng, int cnt):
        id(LatLng2CellId(lat, lng)),
        lat(lat),
        lng(lng),
        samples(cnt),
        is_core(false),
        is_noise(false),
        cluster(-1),
        density(0)
    {
    }

    DPoint(S2CellId id):
        id(id),
        lat(0),
        lng(0),
        samples(1),
        is_core(false),
        is_noise(false),
        cluster(-1),
        density(0)
    {
        S2LatLng ll = id.ToLatLng();
        lat = ll.lat().degrees();
        lng = ll.lat().degrees();
    }

    DPoint(S2CellId id, int cnt):
        id(id),
        lat(0),
        lng(0),
        samples(cnt),
        is_core(false),
        is_noise(false),
        cluster(-1),
        density(0)
    {
        S2LatLng ll = id.ToLatLng();
        lat = ll.lat().degrees();
        lng = ll.lat().degrees();
    }
};

struct Cluster {
    int    cnt;
    double lat;
    double lng;
    Cluster():
        cnt(0),
        lat(0),
        lng(0)
    {
    }
};

struct MobileDetectOption {
    int     min_pts;
    double  max_expand_radius;
    double  min_expand_radius;
    int     dbscan_min_pts;
    double  dbscan_min_radius;
    double  core_scale;
    double  noise_level;
    MobileDetectOption():
        min_pts(20),
        max_expand_radius(10000),
        min_expand_radius(1000),
        dbscan_min_pts(0),
        dbscan_min_radius(25),
        core_scale(3),
        noise_level(0.05)
    {
    }
};

struct PruneOption {
    int     dbscan_min_pts;
    double  dbscan_min_radius;
    double  noise_level;
    PruneOption():
        dbscan_min_pts(1),
        dbscan_min_radius(25),
        noise_level(0.05)
    {
    }
};

inline static bool operator<(const DPoint &r, const DPoint &l) {
    return r.id < l.id;
}

inline static bool operator==(const DPoint &r, const DPoint &l) {
    return r.id == l.id;
}

typedef std::vector<DPoint> DPointVec;
int dbscan(DPointVec &points, double dc, int minPts);
int quick_union(const DPointVec &points, std::vector<Cluster> &unions, double dc);
int density_peak_cluster_fast(DPointVec &points, double dc);
S2LatLng density_peak_centroid(DPointVec &points, double dc);
S2LatLng density_peak_centroid_slow(DPointVec &points, double dc);
int density_peak_filter(DPointVec &points, double threshold = 60, double dc = 10);
int density_peak_filter_slow(DPointVec &points, double threshold=60, double dc=10);

// high-level routines
//
// static double kDefaultUnionRadius    = 1000;
// static int    kDefaultMinPts         = 0;
// static double kDefaultNeiborDistance = 45;
// static double kDefaultCoreScale      = 3;

bool is_mobile(DPointVec &points, MobileDetectOption opt);
void prune(DPointVec &points, PruneOption opt);
void merge_cluster_0(DPointVec &points, double distance);
void merge_cluster_1(DPointVec &points, double distance);

#endif
