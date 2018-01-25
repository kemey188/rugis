#ifndef __H_QMAP_S2MODULE_H__
#define __H_QMAP_S2MODULE_H__

extern "C" {

#include <Python.h>

#include "structmember.h"

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
}

#include <vector>

class S2Loop;
typedef struct {
    PyObject_HEAD
    S2Loop *bound_;
} PyS2Loop;

// custom types
extern PyTypeObject PyS2Loop_Type;

#define PyS2Loop_Check(op) PyObject_TypeCheck(op, &PyS2Loop_Type)

extern PyObject* s2_exception;

#endif
