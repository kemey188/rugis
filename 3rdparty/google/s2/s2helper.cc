#include "s2.h"
#include "s2cap.h"
#include "s2cellid.h"
#include "s2latlng.h"
#include "s2regioncoverer.h"
#include "strings/strutil.h"

const double kEarthCircumferenceMeters = 1000 * 40075.017;

double EarthMetersToRadians(double meters) {
    return (2 * M_PI) * (meters / kEarthCircumferenceMeters);
}

double RadiansToEarthMeters(double radians) {
    return (radians * kEarthCircumferenceMeters) / (2 * M_PI);
}

// Returns the cell level with a side most closely matching the
// specified number of meters.
int GetClosestLevel(double meters) {
    return S2::kAvgEdge.GetClosestLevel(EarthMetersToRadians(meters));
}

// Returns the distance between two locations.
double DistanceBetweenLocations(
        double lat1, double lng1,
        double lat2, double lng2) {
    S2LatLng latlng1 = S2LatLng::FromDegrees(lat1, lng1);
    S2LatLng latlng2 = S2LatLng::FromDegrees(lat2, lng2);
    return RadiansToEarthMeters(latlng1.GetDistance(latlng2).radians());
}

// Generates a list of cells covering lat,lng between the target
// s2 cell levels (min_level, max_level).
int IndexCells(double lat, double lng,
        int min_level, int max_level,
        std::vector<S2CellId> &cells
        ) {
    S2CellId id(S2CellId::FromLatLng(S2LatLng::FromDegrees(lat, lng)));
    do {
        if (id.level() >= min_level && id.level() <= max_level) {
            cells.push_back(id);
        }
        id = id.parent();
    } while (id.level() > 0);
    return cells.size();
}

// Generates a list of cells at the target s2 cell levels which cover
// a cap of radius 'radius_meters' with center at lat & lng.
int NeighborCells(double lat, double lng, double radius_in_meters,
        int min_level, int max_level,
        std::vector<S2CellId>& covering) {
    covering.clear();
    const double radius_radians = EarthMetersToRadians(radius_in_meters);
    // Note: h = 1 - cos(theta) ~ 2 * sin(theta/2) * sin(theta/2)
    //       sin(theta) ~ theta, theta->0
    // But it's not accuracy, why???
    // const S2Cap region = S2Cap::FromAxisHeight(
    //         S2LatLng::FromDegrees(lat, lng).Normalized().ToPoint(),
    //         (radius_radians * radius_radians) / 2);
    const S2Cap region = S2Cap::FromAxisAngle(
            S2LatLng::FromDegrees(lat, lng).Normalized().ToPoint(),
            S1Angle::Radians(radius_radians));
    S2RegionCoverer coverer;
    coverer.set_min_level(min_level);
    coverer.set_max_level(max_level);

    coverer.GetCovering(region, &covering);
    return covering.size();
}

