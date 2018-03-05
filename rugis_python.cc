#include <Python.h>
#include <structmember.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <limits.h>
#include "google/s2/s2helper.h"
#include "google/s2/s2latlng.h"
#include "google/s2/s2loop.h"
#include "google/util/coding/coder.h"
#include "radar.h"
#include "coord.h"
#include "geo.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/point.hpp>

static PyObject* PyIsPointInBox(PyObject* self, PyObject* args) {
  (void)self;
  double box_upleft_lat, box_upleft_lng, box_downright_lat, box_downright_lng, lat, lng;

  if (!PyArg_ParseTuple(args, "dddddd", &box_upleft_lat, &box_upleft_lng, &box_downright_lat, &box_downright_lng, &lat, &lng)) {
    Py_RETURN_NONE;
  }
  int result = geo::isInBox(box_upleft_lat, box_upleft_lng, box_downright_lat, box_downright_lng, lat, lng);
  return PyInt_FromLong(result);
}

static PyObject* PyGetClosestLevel(PyObject* self, PyObject* args) {
  (void)self;
  double meters;
  if (!PyArg_ParseTuple(args, "d", &meters, &meters)) {
    Py_RETURN_NONE;
  }
  int result = GetClosestLevel(meters);
  return PyInt_FromLong(result);
}

static PyObject* PyDistanceOfLL(PyObject* self, PyObject* args) {
  (void)self;
  double lat1, lng1, lat2, lng2;
  if (!PyArg_ParseTuple(args, "dddd", &lat1, &lng1, &lat2, &lng2)) {
    Py_RETURN_NONE;
  }
  double result = DistanceBetweenLocations(lat1, lng1, lat2, lng2);
  return PyFloat_FromDouble(result);
}

static PyObject* PyDistanceOfId(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t c1 = 0, c2 = 0;
  if (!PyArg_ParseTuple(args, "ll", &c1, &c2)) {
    Py_RETURN_NONE;
  }
  double result = RadiansToEarthMeters(
      S2CellId(c1).ToLatLng().GetDistance(
        S2CellId(c2).ToLatLng()).radians());
  return PyFloat_FromDouble(result);
}

static PyObject* PyLatLng2Id(PyObject* self, PyObject* args) {
  (void)self;
  double lat, lng;
  int level = 25;
  if (!PyArg_ParseTuple(args, "ddi", &lat, &lng, &level)) {
    Py_RETURN_NONE;
  }
  S2CellId id = S2CellId::FromLatLng(S2LatLng::FromDegrees(lat, lng));
  return PyInt_FromLong(id.parent(level).id());
}

static PyObject* PyId2LatLng(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  if (!PyArg_ParseTuple(args, "l", &id)) {
    return NULL;
  }
  S2CellId cid(id);
  S2LatLng ll = cid.ToLatLng();
  double lat = ll.lat().degrees();
  double lng = ll.lng().degrees();
  PyObject* list = PyList_New(3);
  PyList_SET_ITEM(list, 0, PyFloat_FromDouble(lat));
  PyList_SET_ITEM(list, 1, PyFloat_FromDouble(lng));
  PyList_SET_ITEM(list, 2, PyInt_FromLong(cid.level()));
  return list;
}

static PyObject* PyLatLng2UTM(PyObject* self, PyObject* args) {
  (void)self;
  double lat,lng;
  if (!PyArg_ParseTuple(args, "dd", &lat, &lng)) {
    Py_RETURN_NONE;
  }
  double x,y;
  geo::LatLng2UTM(lat, lng, x, y);

  PyObject* list = PyList_New(2);
  PyList_SET_ITEM(list, 0, PyFloat_FromDouble(x));
  PyList_SET_ITEM(list, 1, PyFloat_FromDouble(y));
  return list;
}

static PyObject* PyUTM2LatLng(PyObject* self, PyObject* args) {
  (void)self;
  double lat, lng;
  double x = 0,y = 0;
  if (!PyArg_ParseTuple(args, "dd", &x, &y)) {
    Py_RETURN_NONE;
  }
  geo::LatLng2UTM(x, y, lat, lng);

  PyObject* list = PyList_New(2);
  PyList_SET_ITEM(list, 0, PyFloat_FromDouble(lat));
  PyList_SET_ITEM(list, 1, PyFloat_FromDouble(lng));
  return list;
}

static PyObject* PyAdvanceId(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id = 0;
  int step = 0;
  if (!PyArg_ParseTuple(args, "li", &id, &step)) {
    Py_RETURN_NONE;
  }
  S2CellId cid(id);
  cid = cid.advance_wrap(step);
  return PyInt_FromLong(cid.id());
}

static PyObject* PyIndexCells(PyObject* self, PyObject* args) {
  (void)self;
  double lat, lng;
  int min_level, max_level;
  if (!PyArg_ParseTuple(args, "ddii", &lat, &lng, &min_level, &max_level)) {
    Py_RETURN_NONE;
  }
  std::vector<S2CellId> cells;
  IndexCells(lat, lng, min_level, max_level, cells);
  PyObject* list = PyList_New(cells.size());
  for (int i = 0, l = cells.size(); i < l; i++) {
    PyList_SET_ITEM(list, i, PyInt_FromLong(cells[i].id()));
  }
  return list;
}
// lukaimin

static PyObject* PyNeighborCells(PyObject* self, PyObject* args) {
  (void)self;
  double lat, lng, radius;
  int min_level, max_level;

  if (!PyArg_ParseTuple(args, "dddii", &lat, &lng, &radius, &min_level, &max_level)) {
    Py_RETURN_NONE;
  }

  std::vector<S2CellId> cells;
  NeighborCells(lat, lng, radius, min_level, max_level, cells);

  PyObject* list = PyList_New(cells.size());
  for (int i = 0, l = cells.size(); i < l; i++) {
    PyList_SET_ITEM(list, i, PyInt_FromLong(cells[i].id()));
  }
  return list;
}

static PyObject* PyAdjacentCells(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  if (!PyArg_ParseTuple(args, "l", &id)) {
    Py_RETURN_NONE;
  }

  S2CellId cid(id);
  std::vector<S2CellId> cells;
  cid.AppendAllNeighbors(cid.level(), &cells);

  PyObject* list = PyList_New(cells.size());
  for (int i = 0, l = cells.size(); i < l; i++) {
    PyList_SET_ITEM(list, i, PyInt_FromLong(cells[i].id()));
  }
  return list;
}

static PyObject* PyId2QuadKey(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  if (!PyArg_ParseTuple(args, "l", &id)) {
    Py_RETURN_NONE;
  }
  S2CellId cid(id);
  int      face = cid.face();
  uint64_t pos  = cid.pos();
  PyObject* list = PyList_New(cid.level()+1);
  PyList_SET_ITEM(list, 0, PyInt_FromLong(face));
  for(int i = 0; i < cid.level(); ++ i) {
    int k = (pos >> (58 - i*2)) & 0x3;
    PyList_SET_ITEM(list, i+1, PyInt_FromLong(k));
  }
  return list;
}

static PyObject* PyNextId(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  if (!PyArg_ParseTuple(args, "l", &id)) {
    return NULL;
  }
  S2CellId cid(id);
  return PyInt_FromLong(cid.next().id());
}

static PyObject* PyPrevId(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  if (!PyArg_ParseTuple(args, "l", &id)) {
    return NULL;
  }
  S2CellId cid(id);
  return PyInt_FromLong(cid.prev().id());
}

static PyObject* PyParentId(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  int level;
  if (!PyArg_ParseTuple(args, "li", &id, &level)) {
    return NULL;
  }
  S2CellId cid(id);
  return PyInt_FromLong(cid.parent(level).id());
}

static PyObject* PyId2PosFaceLevel(PyObject* self, PyObject* args) {
  (void)self;
  uint64_t id;
  if (!PyArg_ParseTuple(args, "l", &id)) {
    Py_RETURN_NONE;
  }
  S2CellId cid(id);
  PyObject* list = PyList_New(3);
  PyList_SET_ITEM(list, 0, PyInt_FromLong(cid.pos()));
  PyList_SET_ITEM(list, 1, PyInt_FromLong(cid.face()));
  PyList_SET_ITEM(list, 2, PyInt_FromLong(cid.level()));
  return list;
}

// lukaimin
static PyObject *PyGPS2GCJ(PyObject *self, PyObject *args) {
  (void)self;
  double lat, lng;
  if (!PyArg_ParseTuple(args, "dd", &lat, &lng)) {
    Py_RETURN_NONE;
  }
  double xlat = 0, xlng = 0;
  geo::GPS2GCJ(lat, lng, &xlat, &xlng);
  PyObject* list = PyTuple_New(2);
  PyTuple_SET_ITEM(list, 0, PyFloat_FromDouble(xlat));
  PyTuple_SET_ITEM(list, 1, PyFloat_FromDouble(xlng));
  return list;
}

//geometric 
const long double _x_PI_ = 3.14159265358979324 * 3000.0 / 180.0;

// kemey.RU
static PyObject *PyBD092GCJ(PyObject *self, PyObject *args) {
  (void)self;
  double lat, lng;
  if (!PyArg_ParseTuple(args, "dd", &lat, &lng)) {
    Py_RETURN_NONE;
  }
  double gcj_lat = 0, gcj_lng = 0;
  double x = lng - 0.0065;
  double y = lat - 0.006;
  double z = sqrt(x * x + y * y) - 0.00002 * sin(y * _x_PI_);
  double theta = atan2(y, x) - 0.000003 * cos(x * _x_PI_);
  gcj_lng = z * cos(theta);
  gcj_lat = z * sin(theta);  

  PyObject* list = PyTuple_New(2);
  PyTuple_SET_ITEM(list, 0, PyFloat_FromDouble(gcj_lat));
  PyTuple_SET_ITEM(list, 1, PyFloat_FromDouble(gcj_lng));
  return list;
}

static PyObject *PyGCJ2BD09(PyObject *self, PyObject *args) {
  (void)self;
  double lat, lng;
  if (!PyArg_ParseTuple(args, "dd", &lat, &lng)) {
    Py_RETURN_NONE;
  }
  double bd_lat = 0, bd_lng = 0;

  double z = sqrt(lng * lng + lat * lat) + 0.00002 * sin(lat * _x_PI_);
  double theta = atan2(lat, lng) + 0.000003 * cos(lng * _x_PI_);
  bd_lng = z * cos(theta) + 0.0065;
  bd_lat = z * sin(theta) + 0.006;
  
  PyObject* list = PyTuple_New(2);
  PyTuple_SET_ITEM(list, 0, PyFloat_FromDouble(bd_lat));
  PyTuple_SET_ITEM(list, 1, PyFloat_FromDouble(bd_lng));
  return list;
}

// kemey.RU
static bool GetPointsFromArgsTuple3(PyObject *list, std::vector<DPoint> &points) {
  Py_ssize_t size = PyList_Size(list);
  for (int i = 0; i < size; ++i) {
    PyObject *ll = PyList_GetItem(list, i);
    double lat = PyFloat_AsDouble(PyTuple_GetItem(ll, 0));
    double lng = PyFloat_AsDouble(PyTuple_GetItem(ll, 1));
    int cnt = std::max<int>(1, PyInt_AsLong(PyTuple_GetItem(ll, 2)));
    points.push_back(DPoint(lat, lng, cnt));
  }
  return true;
}

static bool GetPointsFromArgsTuple2(PyObject *list, std::vector<DPoint> &points) {
  Py_ssize_t size = PyList_Size(list);
  for (int i = 0; i < size; ++i) {
    PyObject *ll = PyList_GetItem(list, i);
    double lat = PyFloat_AsDouble(PyTuple_GetItem(ll, 0));
    double lng = PyFloat_AsDouble(PyTuple_GetItem(ll, 1));
    points.push_back(DPoint(lat, lng));
  }
  return true;
}

static bool GetPointsFromArgsList3(PyObject *list, std::vector<DPoint> &points) {
  Py_ssize_t size = PyList_Size(list);
  for (int i = 0; i < size; ++i) {
    PyObject *ll = PyList_GetItem(list, i);
    double lat = PyFloat_AsDouble(PyList_GetItem(ll, 0));
    double lng = PyFloat_AsDouble(PyList_GetItem(ll, 1));
    int cnt = std::max<int>(1, PyInt_AsLong(PyList_GetItem(ll, 2)));
    // for(int c = 0; c < cnt; ++ c)
    points.push_back(DPoint(lat, lng, cnt));
  }
  return true;
}

static bool GetPointsFromArgsList2(PyObject *list, std::vector<DPoint> &points) {
  Py_ssize_t size = PyList_Size(list);
  for (int i = 0; i < size; ++i) {
    PyObject *ll = PyList_GetItem(list, i);
    double lat = PyFloat_AsDouble(PyList_GetItem(ll, 0));
    double lng = PyFloat_AsDouble(PyList_GetItem(ll, 1));
    points.push_back(DPoint(lat, lng));
  }
  return true;
}

static bool GetPointsFromArgs(PyObject *args, std::vector<DPoint> &points) {
  PyObject *list;
  if (!PyArg_ParseTuple(args, "O", &list)) {
    return false;
  }
  if (!PyList_Check(list)) {
    return false;
  }
  Py_ssize_t size = PyList_Size(list);
  if (size == 0) return false;

  PyObject *ll = PyList_GetItem(list, 0);
  if (PyTuple_Check(ll)) {
    if (3 <= PyTuple_Size(ll)) 
      return GetPointsFromArgsTuple3(list, points);
    else
      return GetPointsFromArgsTuple2(list, points);
  }
  else {
    if (3 <= PyList_Size(ll)) 
      return GetPointsFromArgsList3(list, points);
    else 
      return GetPointsFromArgsList2(list, points);
  }
}

static PyObject *PyIsMobile(PyObject *self, PyObject *args) {
  (void)self;
  std::vector<DPoint> points;
  if (!GetPointsFromArgs(args, points)) 
    return NULL;
  MobileDetectOption opt;
  return PyInt_FromLong(is_mobile(points, opt));
}

static PyObject *PyPrune(PyObject *self, PyObject *args) {
  (void)self;
  std::vector<DPoint> points;
  if (!GetPointsFromArgs(args, points)) 
    Py_RETURN_NONE;
  PruneOption opt;
  prune(points, opt);
  PyObject *list = PyList_New(points.size());
  int i = 0;
  for(auto pt: points) {
    PyObject *e = PyList_New(3);
    PyList_SET_ITEM(e, 0, PyFloat_FromDouble(pt.lat));
    PyList_SET_ITEM(e, 1, PyFloat_FromDouble(pt.lng));
    PyList_SET_ITEM(e, 2, PyInt_FromLong(pt.samples));
    PyList_SET_ITEM(list, i, e);
    i ++;
  }
  return list;
}

static PyObject *PyDBScan(PyObject *self, PyObject *args) {
  (void)self;
  PyObject *list;
  int dbscan_min_pts = 1;
  double dbscan_min_radius = 45;
  double merge_factor = 1.5;

  if (!PyArg_ParseTuple(args, "Oidd",
        &list, &dbscan_min_pts,
        &dbscan_min_radius,
        &merge_factor)) {
    return NULL;
  }
  // printf("%d, %g\n", dbscan_min_pts, dbscan_min_radius);
  if (!PyList_Check(list)) {
    Py_RETURN_NONE;
  }
  Py_ssize_t size = PyList_Size(list);
  if (size == 0) return NULL;

  std::vector<DPoint> points;
  points.reserve(size);

  PyObject *ll = PyList_GetItem(list, 0);
  if (PyTuple_Check(ll)) {
    if (3 <= PyTuple_Size(ll)) 
      GetPointsFromArgsTuple3(list, points);
    else
      GetPointsFromArgsTuple2(list, points);
  }
  else {
    if (3 <= PyList_Size(ll)) 
      GetPointsFromArgsList3(list, points);
    else 
      GetPointsFromArgsList2(list, points);
  }

  dbscan(points, dbscan_min_radius, dbscan_min_pts);
  if (merge_factor > 1) {
    merge_cluster_0(points, dbscan_min_radius);
    merge_cluster_1(points, dbscan_min_radius*merge_factor);
  }

  PyObject *cls = PyList_New(points.size());
  for(int i = 0, l = points.size(); i < l; ++ i) {
    PyList_SET_ITEM(cls, i, PyInt_FromLong(points[i].cluster));
  }
  // printf("%d - %d - %d\n", size, points.size(), PyList_Size(cls));
  return cls;
}

static PyObject *PyLatLng2XY(PyObject *self, PyObject *args) {
  (void)self;
  double lat, lng;
  int level;
  if (!PyArg_ParseTuple(args, "ddi", &lat, &lng, &level)) {
    Py_RETURN_NONE;
  }
  MapPoint pt(lat, lng);
  PyObject *r = PyList_New(2);
  PyList_SET_ITEM(r, 0, PyInt_FromLong(pt.x));
  PyList_SET_ITEM(r, 1, PyInt_FromLong(pt.y));
  return r;
}

static PyObject *PyBoundBox(PyObject *self, PyObject *args) {
  (void)(self);

  PyObject *pts;
  if (!PyArg_ParseTuple(args, "O", &pts)) {
    Py_RETURN_NONE;
  }
  int len = PyList_Size(pts);

  std::vector<S2Point> points;
  points.reserve(len);

  for(int i = 0; i < len; ++ i) {
    PyObject *t = PyList_GetItem(pts, i);
    if (PyTuple_Check(t)) {
      double lat = PyFloat_AsDouble(PyTuple_GetItem(t, 0));
      double lng = PyFloat_AsDouble(PyTuple_GetItem(t, 1));
      points.push_back(S2LatLng::FromDegrees(lat, lng).ToPoint());
    } else {
      double lat = PyFloat_AsDouble(PyList_GetItem(t, 0));
      double lng = PyFloat_AsDouble(PyList_GetItem(t, 1));
      points.push_back(S2LatLng::FromDegrees(lat, lng).ToPoint());
    } 
  }

  S2Loop loop;
  loop.Init(points);
  loop.Normalize();

  S2LatLngRect bound = loop.GetRectBound();

  if (bound.is_full()) {
    loop.Invert();
    bound = loop.GetRectBound();
  }

  PyObject *res = PyList_New(4);
  PyList_SET_ITEM(res, 0, PyFloat_FromDouble(bound.lat_lo().degrees()));
  PyList_SET_ITEM(res, 1, PyFloat_FromDouble(bound.lng_lo().degrees()));
  PyList_SET_ITEM(res, 2, PyFloat_FromDouble(bound.lat_hi().degrees()));
  PyList_SET_ITEM(res, 3, PyFloat_FromDouble(bound.lng_hi().degrees()));

  return res;
}

static PyObject *PyEncodeLoop(PyObject *self, PyObject *args) {
  (void)(self);

  PyObject *pts;
  if (!PyArg_ParseTuple(args, "O", &pts)) {
    Py_RETURN_NONE;
  }
  int len = PyList_Size(pts);

  std::vector<S2Point> points;
  points.reserve(len);

  for(int i = 0; i < len; ++ i) {
    PyObject *t = PyList_GetItem(pts, i);
    if (PyTuple_Check(t)) {
      double lat = PyFloat_AsDouble(PyTuple_GetItem(t, 0));
      double lng = PyFloat_AsDouble(PyTuple_GetItem(t, 1));
      points.push_back(S2LatLng::FromDegrees(lat, lng).ToPoint());
    } else {
      double lat = PyFloat_AsDouble(PyList_GetItem(t, 0));
      double lng = PyFloat_AsDouble(PyList_GetItem(t, 1));
      points.push_back(S2LatLng::FromDegrees(lat, lng).ToPoint());
    } 
  }

  S2Loop loop;
  loop.Init(points);
  loop.Normalize();

  S2LatLngRect bound = loop.GetRectBound();
  if (bound.is_full()) {
    loop.Invert();
    bound = loop.GetRectBound();
  }

  Encoder coder;
  loop.Encode(&coder);

  PyObject *res = PyList_New(5);
  PyList_SET_ITEM(res, 0, PyFloat_FromDouble(bound.lat_lo().degrees()));
  PyList_SET_ITEM(res, 1, PyFloat_FromDouble(bound.lng_lo().degrees()));
  PyList_SET_ITEM(res, 2, PyFloat_FromDouble(bound.lat_hi().degrees()));
  PyList_SET_ITEM(res, 3, PyFloat_FromDouble(bound.lng_hi().degrees()));
  PyList_SET_ITEM(res, 4, PyString_FromStringAndSize(coder.base(), coder.length()));

  return res;
}

static PyObject* PySimplify(PyObject* self, PyObject* args) {
  (void)self;
  typedef boost::geometry::model::point<double, 2, boost::geometry::cs::spherical_equatorial<boost::geometry::degree> > latlng;

  PyObject *pts;
  double distance = 100;

  if (!PyArg_ParseTuple(args, "Od", &pts, &distance)) {
    Py_RETURN_NONE;
  }

  boost::geometry::model::linestring<latlng> line;
  int len = PyList_Size(pts);
  for(int i = 0; i < len; ++ i) {
    PyObject *t = PyList_GetItem(pts, i);
    if (PyTuple_Check(t)) {
      double lat = PyFloat_AsDouble(PyTuple_GetItem(t, 0));
      double lng = PyFloat_AsDouble(PyTuple_GetItem(t, 1));
      line.push_back(latlng(lat, lng));
    } else {
      double lat = PyFloat_AsDouble(PyList_GetItem(t, 0));
      double lng = PyFloat_AsDouble(PyList_GetItem(t, 1));
      line.push_back(latlng(lat, lng));
    } 
  }

  // Simplify it, using distance of 0.5 units
  boost::geometry::model::linestring<latlng> simplified;
  boost::geometry::simplify(line, simplified, distance/6378137.0);

  PyObject *res = PyList_New(simplified.size());
  for(size_t i = 0; i < simplified.size(); ++ i) {
    PyObject *t = PyList_New(2);
    PyList_SET_ITEM(t, 0, PyFloat_FromDouble(simplified[i].get<0>()));
    PyList_SET_ITEM(t, 1, PyFloat_FromDouble(simplified[i].get<1>()));
    PyList_SET_ITEM(res, i, t); 
  }
  return res;
}

static PyMethodDef kMethods[] = {
  {"Id2PosFaceLevel", PyId2PosFaceLevel, METH_VARARGS, "Id2PosFaceLevel(id) -> (pos, face, level)"},
  {"ClosestLevel", PyGetClosestLevel, METH_VARARGS, "GetClosestLevel(distance) -> level"},
  {"Distance",     PyDistanceOfLL,    METH_VARARGS, "Distance(lat1, lng1, lat2, lng2)"},
  {"DistanceOfId", PyDistanceOfId,    METH_VARARGS, "Distance(cell_id1, cell_id2)"},
  {"AdvanceId",    PyAdvanceId,       METH_VARARGS, "AdvanceId(cell_id, step)"},
  {"LatLng2Id",    PyLatLng2Id,       METH_VARARGS, "LatLng2Id(lat, lng, level) -> id"},
  {"Id2LatLng",    PyId2LatLng,       METH_VARARGS, "Id2LatLng(id) -> (lat, lng, level)"},
  {"LatLng2UTM",   PyLatLng2UTM,      METH_VARARGS, "LatLng2UTM(lat, lng) -> x,y"},
  {"UTM2LatLng",   PyLatLng2UTM,      METH_VARARGS, "UTM2LatLng(x, y) -> lat,lng"},
  {"LatLng2XY",    PyLatLng2XY,       METH_VARARGS, "LatLng2XY(lat, lng) -> i,j"},
  {"Id2QuadKey",   PyId2QuadKey,      METH_VARARGS, "Id2QuadKey(id) -> (face, ..)"},
  {"Indexes",      PyIndexCells,      METH_VARARGS, "Indexes(lat, lng, max_level, min_level) -> (cell_id, ..)"},
  {"Neighbors",    PyNeighborCells,   METH_VARARGS, "Neighbors(lat, lng, radius, min_level, max_level) -> (cell_id, ..)"},
  {"Adjacent",     PyAdjacentCells,   METH_VARARGS, "Adjacent(id) -> (cell_id, ..)\nreturn 8 adjacent cells"},
  {"Next",         PyNextId,          METH_VARARGS, "Next(id) -> next cell id in hillbert curve space"},
  {"Prev",         PyPrevId,          METH_VARARGS, "Prev(id) -> next cell id in hillbert curve space"},
  {"Parent",       PyParentId,        METH_VARARGS, "Parent(id) -> parent id at level"},
  {"BoundBox",     PyBoundBox,        METH_VARARGS, "BBox((lat1, lng1), ...) -> lat_lo, lng_lo, lat_hi, lng_hi"},
  {"IsPointInBox", PyIsPointInBox,    METH_VARARGS|METH_KEYWORDS, "IsPointInBox(box_upleft_lat,box_upleft_lng,box_downright_lat,box_downright_lng,lat,lng) -> 0/1\nIsPointInBox(39.98596,116.47112,39.9852,116.47277,39.9856,116.47192)"},
  {"EncodeArea",   PyEncodeLoop,      METH_VARARGS, "EncodeLoop((lat1, lng1), ...)"},
  {"Simplify",     PySimplify,        METH_VARARGS, "PySimplify((lat1, lng1), ...)"},
  {"GPS2GCJ",      PyGPS2GCJ,         METH_VARARGS, ""},
  {"GCJ2BD09",     PyGCJ2BD09,        METH_VARARGS, "GCJ2BD09(gcj_lat, gcj_lng) -> (bd_lat, bd_lng)"},
  {"BD092GCJ",     PyBD092GCJ,        METH_VARARGS, "BD092GCJ(bd_lat, bd_lng) -> (gcj_lat, gcj_lng)"},
  {"IsMobile",     PyIsMobile,        METH_VARARGS|METH_KEYWORDS, ""},
  {"Prune",        PyPrune,           METH_VARARGS|METH_KEYWORDS, ""},
  {"dbscan",       PyDBScan,          METH_VARARGS|METH_KEYWORDS, ""},
  {NULL, NULL, 0, NULL}
};


static PyObject* _exception = NULL;
PyMODINIT_FUNC initqgis(void) {
  PyObject* _module = Py_InitModule((char*)"qgis", kMethods);
  if (_module == 0) {
    return;
  }
  
  PyObject *mversion = PyString_FromString("2.0.0");  
  PyObject *mauthor  = PyString_FromString("kemey@163.com"); 
  
  PyModule_AddObject(_module, "Version", mversion);
  PyModule_AddObject(_module, "Author", mauthor);   

  // add custom exception
  _exception = PyErr_NewException((char*)"qgis.Error", 0, 0);

  if (_exception == 0) {
    Py_DECREF(_module);
    return;
  }

  if (PyModule_AddObject(_module, (char*)"Error", _exception) != 0) {
    Py_DECREF(_module);
    return;
  }

}
