// Copyright (c) Arni Mar Jonsson.
// See LICENSE for details.

#include "s2module.h"
#include "s2loop.h"
#include "s2latlng.h"


static void PyS2_set_error(const char *err)
{
	PyErr_SetString(s2_exception, err);
}

static void PyS2Loop_dealloc(PyS2Loop* self)
{
    ((PyObject*)self)->ob_type->tp_free((PyObject*)self);
    delete self->bound_;
    self->bound_ = NULL;
}

static PyObject* PyS2Loop_new(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
	PyS2Loop* self = (PyS2Loop*)type->tp_alloc(type, 0);
    if (self) {
        self->bound_ = new S2Loop();
    }
    return (PyObject*)self;
}

static bool Contains(S2Loop &bound, double lat, double lng) {
    if (lng < 72.004 || lng > 137.8347)  
        return false;  
    if (lat < 20.231 || lat > 55.8271)  
        return false;  
    S2Point p = S2LatLng::FromDegrees(lat, lng).Normalized().ToPoint();
    return bound.Contains(p);
}

static PyObject* PyS2Loop_Contains(PyS2Loop* self, PyObject* args, PyObject* kwds)
{
    double lat, lng;
	if (!PyArg_ParseTuple(args, "dd", &lat, &lng))
		return NULL;

    bool y = Contains(*(self->bound_), lat, lng);
    return PyInt_FromLong(y?1:0);
}

static PyObject* PyS2Loop_ContainsId(PyS2Loop* self, PyObject* args, PyObject* kwds)
{
    uint64_t id = 0;
	if (!PyArg_ParseTuple(args, "l", &id))
		return NULL;

    S2CellId cid(id);
    S2LatLng ll = cid.ToLatLng();
    double lat = ll.lat().degrees();
    double lng = ll.lng().degrees();

    bool y = Contains(*(self->bound_), lat, lng);
    return PyInt_FromLong(y?1:0);
}


static PyObject* PyS2Loop_Init(PyS2Loop* self, PyObject* args, PyObject* kwds) {
    PyObject *list;
    if (!PyArg_ParseTuple(args, "O", &list)) {
        return PyInt_FromLong(0);
    }
    if (!PyList_Check(list)) {
        return PyInt_FromLong(0);
    }
    std::vector<S2Point> points;
    Py_ssize_t size = PyList_Size(list);
    for(int i = 0; i < size; ++ i) {
        PyObject *ll = PyList_GetItem(list, i);
        if (PyList_Check(ll)) {
            double lat = PyFloat_AsDouble(PyList_GetItem(ll, 0));
            double lng = PyFloat_AsDouble(PyList_GetItem(ll, 1));
            points.push_back(S2LatLng::FromDegrees(lat, lng).Normalized().ToPoint());
        }
        else if (PyTuple_Check(ll)) {
            double lat = PyFloat_AsDouble(PyTuple_GetItem(ll, 0));
            double lng = PyFloat_AsDouble(PyTuple_GetItem(ll, 1));
            points.push_back(S2LatLng::FromDegrees(lat, lng).Normalized().ToPoint());
        }
    }
    if (points.size() < 3) {
        return PyInt_FromLong(0);
    }
    self->bound_->Init(points);
    if (self->bound_->is_hole()) {
        return PyInt_FromLong(0);
    }
    return PyInt_FromLong(1);
}

static PyMethodDef PyS2Loop_methods[] = {
	{(char*)"Contains",    (PyCFunction)PyS2Loop_Contains,      METH_VARARGS, (char*)"" },
	{(char*)"ContainsId",  (PyCFunction)PyS2Loop_ContainsId,    METH_VARARGS, (char*)"" },
	{(char*)"Init",        (PyCFunction)PyS2Loop_Init,          METH_VARARGS, (char*)"Init([(lat, lng),...])" },
	{NULL}
};

static int PyS2Loop_init(PyS2Loop* self, PyObject* args, PyObject* kwds)
{
	return 0;
}

PyDoc_STRVAR(PyS2Loop_doc, "");

PyTypeObject PyS2Loop_Type = {
	PyObject_HEAD_INIT(NULL)
	0,
	(char*)"radar.S2Loop",      /*tp_name*/
	sizeof(PyS2Loop),           /*tp_basicsize*/
	0,                             /*tp_itemsize*/
	(destructor)PyS2Loop_dealloc, /*tp_dealloc*/
	0,                             /*tp_print*/
	0,                             /*tp_getattr*/
	0,                             /*tp_setattr*/
	0,                             /*tp_compare*/
	0,                             /*tp_repr*/
	0,                             /*tp_as_number*/
	0,                             /*tp_as_sequence*/
	0,                             /*tp_as_mapping*/
	0,                             /*tp_hash */
	0,                             /*tp_call*/
	0,                             /*tp_str*/
	0,                             /*tp_getattro*/
	0,                             /*tp_setattro*/
	0,                             /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,            /*tp_flags*/
	(char*)PyS2Loop_doc,        /*tp_doc */
	0,                             /*tp_traverse */
	0,                             /*tp_clear */
	0,                             /*tp_richcompare */
	0,                             /*tp_weaklistoffset */
	0,                             /*tp_iter */
	0,                             /*tp_iternext */
	PyS2Loop_methods,           /*tp_methods */
	0,                             /*tp_members */
	0,                             /*tp_getset */
	0,                             /*tp_base */
	0,                             /*tp_dict */
	0,                             /*tp_descr_get */
	0,                             /*tp_descr_set */
	0,                             /*tp_dictoffset */
	(initproc)PyS2Loop_init,    /*tp_init */
	0,                             /*tp_alloc */
	PyS2Loop_new,               /*tp_new */
};
