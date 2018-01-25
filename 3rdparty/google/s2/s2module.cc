#include <Python.h>
#include "s2helper.h"
#include "s2latlng.h"
#include "s2module.h"

static PyObject* PyGetClosestLevel(PyObject* self, PyObject* args) {
    double meters;

    if (!PyArg_ParseTuple(args, "d", &meters, &meters)) {
        return NULL;
    }

    int result = GetClosestLevel(meters);

    return PyInt_FromLong(result);
}

static PyObject* PyDistanceOfLL(PyObject* self, PyObject* args) {
    double lat1, lng1, lat2, lng2;

    if (!PyArg_ParseTuple(args, "dddd", &lat1, &lng1, &lat2, &lng2)) {
        return NULL;
    }

    double result = DistanceBetweenLocations(lat1, lng1, lat2, lng2);

    return PyFloat_FromDouble(result);
}

static PyObject* PyDistanceOfId(PyObject* self, PyObject* args) {
    uint64_t c1 = 0, c2 = 0;

    if (!PyArg_ParseTuple(args, "ll", &c1, &c2)) {
        return NULL;
    }

    double result = RadiansToEarthMeters(
            S2CellId(c1).ToLatLng().GetDistance(
                S2CellId(c2).ToLatLng()).radians());

    return PyFloat_FromDouble(result);
}


static PyObject* PyIndexCells(PyObject* self, PyObject* args) {
    double lat, lng;
    int min_level, max_level;

    if (!PyArg_ParseTuple(args, "ddii", &lat, &lng, &min_level, &max_level)) {
        return NULL;
    }

    std::vector<S2CellId> cells;
    IndexCells(lat, lng, min_level, max_level, cells);

    PyObject* list = PyList_New(cells.size());
    for (int i = 0; i < cells.size(); i++) {
        PyList_SET_ITEM(list, i, PyInt_FromLong(cells[i].id()));
    }
    return list;
}

static PyObject* PyNeighborCells(PyObject* self, PyObject* args) {
    double lat, lng, radius;
    int min_level, max_level;

    if (!PyArg_ParseTuple(args, "dddii", &lat, &lng, &radius, &min_level, &max_level)) {
        return NULL;
    }

    std::vector<S2CellId> cells;
    NeighborCells(lat, lng, radius, min_level, max_level, cells);

    PyObject* list = PyList_New(cells.size());
    for (int i = 0; i < cells.size(); i++) {
        PyList_SET_ITEM(list, i, PyInt_FromLong(cells[i].id()));
    }
    return list;
}


static PyObject* PyAdjacentCells(PyObject* self, PyObject* args) {
    uint64_t id;
    if (!PyArg_ParseTuple(args, "l", &id)) {
        return NULL;
    }

    S2CellId cid(id);
    std::vector<S2CellId> cells;
    cid.AppendAllNeighbors(cid.level(), &cells);

    PyObject* list = PyList_New(cells.size());
    for (int i = 0; i < cells.size(); i++) {
        PyList_SET_ITEM(list, i, PyInt_FromLong(cells[i].id()));
    }
    return list;
}

static PyObject* PyLatLng2CellId(PyObject* self, PyObject* args) {
    double lat, lng;
    int level = 25;
    if (!PyArg_ParseTuple(args, "ddi", &lat, &lng, &level)) {
        return NULL;
    }
    S2CellId id = S2CellId::FromLatLng(S2LatLng::FromDegrees(lat, lng));
    return PyInt_FromLong(id.parent(level).id());
}

static PyObject* PyCellId2LatLng(PyObject* self, PyObject* args) {
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
    

static PyObject* PyCellId2QuadKey(PyObject* self, PyObject* args) {
    uint64_t id;
    if (!PyArg_ParseTuple(args, "l", &id)) {
        return NULL;
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
    uint64_t id;
    if (!PyArg_ParseTuple(args, "l", &id)) {
        return NULL;
    }
    S2CellId cid(id);
    return PyInt_FromLong(cid.next().id());
}


static PyObject* PyPrevId(PyObject* self, PyObject* args) {
    uint64_t id;
    if (!PyArg_ParseTuple(args, "l", &id)) {
        return NULL;
    }
    S2CellId cid(id);
    return PyInt_FromLong(cid.prev().id());
}

static PyObject* PyParentId(PyObject* self, PyObject* args) {
    uint64_t id;
    int level;
    if (!PyArg_ParseTuple(args, "li", &id, &level)) {
        return NULL;
    }
    S2CellId cid(id);
    return PyInt_FromLong(cid.parent(level).id());
}

static PyObject* PyCellId2PosFaceLevel(PyObject* self, PyObject* args) {
    uint64_t id;
    if (!PyArg_ParseTuple(args, "l", &id)) {
        return NULL;
    }
    S2CellId cid(id);
    PyObject* list = PyList_New(3);
    PyList_SET_ITEM(list, 0, PyInt_FromLong(cid.pos()));
    PyList_SET_ITEM(list, 1, PyInt_FromLong(cid.face()));
    PyList_SET_ITEM(list, 2, PyInt_FromLong(cid.level()));
    return list;
}

static PyMethodDef kS2Methods[] = {
    {"ClosestLevel", PyGetClosestLevel, METH_VARARGS, ""},
    {"Distance", PyDistanceOfLL, METH_VARARGS, "Distance(lat1, lng1, lat2, lng2)"},
    {"DistanceOfId", PyDistanceOfId, METH_VARARGS, "Distance(cell_id1, cell_id2)"},
    {"LatLng2Id", PyLatLng2CellId, METH_VARARGS, "LatLng2Id(lat, lng, level) -> id"},
    {"Id2LatLng", PyCellId2LatLng, METH_VARARGS, "Id2LatLng(id) -> (lat, lng, level)"},
    {"Id2QuadKey", PyCellId2QuadKey, METH_VARARGS, "Id2QuadKey(id) -> (face, ..)"},
    {"Id2PosFaceLevel", PyCellId2PosFaceLevel, METH_VARARGS, "Id2PosFaceLevel(id) -> (pos, face, level)"},
    {"Indexes", PyIndexCells, METH_VARARGS, "Indexes(lat, lng, max_level, min_level) -> (cell_id, ..)"},
    {"Neighbors", PyNeighborCells, METH_VARARGS, "Neighbors(lat, lng, radius, min_level, max_level) -> (cell_id, ..)"},
    {"Adjacent", PyAdjacentCells, METH_VARARGS, "Ajacent(id) -> (cell_id, ..)\nreturn 8 adjacent cells"},
    {"Next", PyNextId, METH_VARARGS, "Next(id) -> next cell id in hillbert curve space"},
    {"Prev", PyPrevId, METH_VARARGS, "Prev(id) -> next cell id in hillbert curve space"},
    {"Parent", PyParentId, METH_VARARGS, "Parent(id) -> parent id at level"},
    {NULL, NULL, 0, NULL},
};

PyObject* s2_exception = 0;

PyMODINIT_FUNC inits2(void) {
    // Py_InitModule("s2", kS2Methods);
    PyObject* _module = Py_InitModule3((char*)"s2", kS2Methods, 0);
    if (_module == 0) {
        return;
    }

    // add custom exception
    s2_exception = PyErr_NewException((char*)"s2.S2Error", 0, 0);

    if (s2_exception == 0) {
        Py_DECREF(_module);
        return;
    }

    if (PyModule_AddObject(_module, (char*)"S2Error", s2_exception) != 0) {
        Py_DECREF(_module);
        return;
    }

    if (PyType_Ready(&PyS2Loop_Type) < 0) {
        Py_DECREF(_module);
        return;
    }

    // add custom types to the different modules
    Py_INCREF(&PyS2Loop_Type);

    if (PyModule_AddObject(_module, (char*)"S2Loop", (PyObject*)&PyS2Loop_Type) != 0) {
        Py_DECREF(_module);
        return;
    }
}
