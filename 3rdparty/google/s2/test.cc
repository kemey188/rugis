#include <stdlib.h>
#include <stdio.h>
#include "s2.h"
#include "s2latlng.h"
#include "s2cellid.h"
#include "s2cell.h"

static inline double area(S2CellId sid)
{
    return S2Cell(sid).ExactArea() * 40680631590769;
}

int main(int , char **argv)
{
    //Build Sequence: lat,lng -> x,y,z -> f,u,v -> f,s,t -> cellid
    double lat = strtod(argv[1], NULL);
    double lng = strtod(argv[2], NULL);
    double level = strtoul(argv[3], NULL, 10);

    S2LatLng s2ll = S2LatLng::FromDegrees(lat, lng);
    S2CellId cid = S2CellId::FromLatLng(s2ll);
    const S2Point &point = s2ll.ToPoint();

    printf("%.6f,%.6f => %.6f,%.6f => %.6f,%.6f,%.6f => %llx\t%g\n",
            lat, lng,
            s2ll.lat().radians(),
            s2ll.lng().radians(),
            point.x(),
            point.y(),
            point.z(),
            cid.id(),
            area(cid)
            );

    S2CellId cid2 = cid.parent(level);

    const S2Point &point2 = cid2.ToPoint();
    S2LatLng s2ll2 = cid2.ToLatLng();

    printf("%.6f,%.6f => %f,%f => %.6f,%.6f,%.6f => %llx\t%g\n",
            lat, lng,
            s2ll2.lat().degrees(),
            s2ll2.lng().degrees(),
            point2.x(),
            point2.y(),
            point2.z(),
            cid2.id(),
            area(cid2)
            );

    return 0;
}
