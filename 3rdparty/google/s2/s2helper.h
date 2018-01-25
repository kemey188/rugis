#ifndef _H_QMAP_S2HELPER_H_
#define _H_QMAP_S2HELPER_H_

#include "s2cellid.h"
#include <vector>

// Convert between distance of earth and radians 
double EarthMetersToRadians(double meters);
double RadiansToEarthMeters(double radians);

// Returns the cell level with a side most closely matching the
// specified number of meters.
int GetClosestLevel(double meters);

// Returns the distance between two locations.
double DistanceBetweenLocations(
        double lat1, double lng1,
        double lat2, double lng2);

// Generates a list of cells covering lat,lng between the target
// s2 cell levels (min_level, max_level).
int IndexCells(double lat, double lng,
        int min_level, int max_level,
        std::vector<S2CellId> &cells
        );


// Generates a list of cells at the target s2 cell levels which cover
// a cap of radius 'radius_meters' with center at lat & lng.
int NeighborCells(double lat, double lng,
        double radius_in_meters,
        int min_level, int max_level,
        std::vector<S2CellId>& covering);

#endif
