# rugis
A python library for geography science

***FUNCTIONS***:

    Adjacent(...)
        Ajacent(id) -> (cell_id, ..)
        return 8 adjacent cells
    AdvanceId(...)
        AdvanceId(cell_id, step)
    BoundBox(...)
        BBox((lat1, lng1), ...) -> lat_lo, lng_lo, lat_hi, lng_hi
    ClosestLevel(...)
        GetClosestLevel(distance) -> level
    Distance(...)
        Distance(lat1, lng1, lat2, lng2)
    DistanceOfId(...)
        Distance(cell_id1, cell_id2)
    EncodeArea(...)
        EncodeLoop((lat1, lng1), ...)
    GPS2GCJ(...)
    Id2LatLng(...)
        Id2LatLng(id) -> (lat, lng, level)
    Id2PosFaceLevel(...)
        Id2PosFaceLevel(id) -> (pos, face, level)
    Id2QuadKey(...)
        Id2QuadKey(id) -> (face, ..)
    Indexes(...)
        Indexes(lat, lng, max_level, min_level) -> (cell_id, ..)
    LatLng2Id(...)
        LatLng2Id(lat, lng, level) -> id
    LatLng2UTM(...)
        LatLng2UTM(lat, lng) -> x,y
    LatLng2XY(...)
        LatLng2XY(lat, lng) -> i,j
    Neighbors(...)
        Neighbors(lat, lng, radius, min_level, max_level) -> (cell_id, ..)
    Next(...)
        Next(id) -> next cell id in hillbert curve space
    Parent(...)
        Parent(id) -> parent id at level
    Prev(...)
        Prev(id) -> next cell id in hillbert curve space
    Simplify(...)
        PySimplify((lat1, lng1), ...)
    UTM2LatLng(...)
        UTM2LatLng(x, y) -> lat,lng
    dbscan(...)

# Install	
	python setup.py build
	
After build, you will find rugis.so in **./build/lib.linux-x86_64-2.7/**

Bingo, start your projection :)
