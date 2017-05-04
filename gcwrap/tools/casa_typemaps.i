%{
#include <memory>
#include <string>
#include <vector>
#include <complex>
#include <stdcasa/record.h>
#include <tools/swigconvert_python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#if PY_MAJOR_VERSION >= 3
# define PYTEXT_CHECK PyUnicode_Check
# define PYTEXT_ASDATA PyUnicode_AsUTF8String
# define PYBYTES_CHECK PyBytes_Check
# define PYBYTES_ASDATA PyBytes_AsString
# define PYSTRING_FROM_STRING PyUnicode_FromString
#else
# define PYTEXT_CHECK PyUnicode_Check
# define PYTEXT_ASDATA PyUnicode_AsUTF8String
# define PYBYTES_CHECK PyString_Check
# define PYBYTES_ASDATA PyString_AsString
# define PYSTRING_FROM_STRING PyString_FromString
#endif

using casac::record;
using casac::variant;
using namespace casac;
%}

// As far as I can tell, the "deleter" construct adds a method-scope
// unique_ptr to delete the object after the method exits.

// Integer scalars

%typemap(in) int {
    $1 = PyInt_AsLong($input);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be an integer");
        return NULL;
    }
}

%typemap(in) long {
    $1 = PyInt_AsLong($input);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be an integer");
        return NULL;
    }
}

%typemap(in) long long {
    $1 = PyInt_AsLong($input);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be an integer");
        return NULL;
    }
}

%typemap(argout) int& OUTARGINT {
    PyObject *o = PyLong_FromLong(*$1);

    if (!$result || $result == Py_None) {
        $result = o;
    } else {
        // What the hell is this code path?
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result, 0, o2);
        }

        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3, 0, o);
        o2 = $result;
        $result = PySequence_Concat(o2, o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}


// Floating scalars

%typemap(in) float {
    $1 = PyFloat_AsDouble($input);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be a number");
        return NULL;
    }
}

%typemap(in) double {
    $1 = PyFloat_AsDouble($input);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be a number");
        return NULL;
    }
}

%typemap(in) complex {
    Py_complex c = PyComplex_AsCComplex($input);
    if (PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be a number");
        return NULL;
    }

    $1 = std::complex<double>(c.real, c.imag);
}

%typemap(out) complex {
    $result = PyComplex_FromDouble($1.real(), $1.imag());
}

%typemap(argout) double& OUTARGDBL {
    PyObject *o = PyFloat_FromDouble(*$1);

    if (!$result || $result == Py_None) {
        $result = o;
    } else {
        // Again: what??
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result, 0, o2);
        }

        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3, 0, o);
        o2 = $result;
        $result = PySequence_Concat(o2, o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}


// std::vector<> of numeric scalars

%typemap(in) std::vector<bool>& (std::unique_ptr<std::vector<bool> > deleter) {
    if (!$1) {
        deleter.reset (new std::vector<bool>(0));
        $1 = deleter.get();
    } else {
        $1->resize(0);
    }

    std::vector<int> shape;

    if (casac::pyarray_check($input)) {
        casac::numpy2vector((PyArrayObject*) $input, *$1, shape);
    } else {
        if (PYTEXT_CHECK($input)) {
            $1->push_back(0);
            PyErr_SetString(PyExc_TypeError, "argument $1_name must be a string");
            return NULL;
        } else if (PyBool_Check($input)) {
            $1->push_back(bool(PyInt_AsLong($input)));
        } else if (PyInt_Check($input)) {
            $1->push_back(bool(PyInt_AsLong($input)));
        } else if (PyLong_Check($input)) {
            $1->push_back(bool(PyLong_AsLong($input)));
        } else if (PyFloat_Check($input)) {
            $1->push_back(bool(PyInt_AsLong(PyNumber_Long($input))));
        } else {
            shape.push_back(PyList_Size($input));
            casac::pylist2vector($input,  *$1, shape);
        }
    }
}

%typemap(in) std::vector<int>& (std::unique_ptr<std::vector<int> > deleter) {
    if (!$1) {
        deleter.reset (new std::vector<int>(0));
        $1 = deleter.get();
    } else {
        $1->resize(0);
    }

    std::vector<int> shape;

    if (casac::pyarray_check($input)) {
        casac::numpy2vector((PyArrayObject*) $input, *$1, shape);
    } else {
        if (PYTEXT_CHECK($input)) {
            $1->push_back(-1);
            PyErr_SetString(PyExc_TypeError, "argument $1_name must not be a string");
            return NULL;
        } else if (PyInt_Check($input)) {
            $1->push_back(int(PyInt_AsLong($input)));
        } else if (PyLong_Check($input)) {
            $1->push_back(PyLong_AsLong($input));
        } else if (PyFloat_Check($input)) {
            $1->push_back(PyInt_AsLong(PyNumber_Long($input)));
        } else {
            shape.push_back(PyList_Size($input));
            casac::pylist2vector($input,  *$1, shape);
        }
    }
}

%typemap(in) std::vector<long>& (std::unique_ptr<std::vector<long> > deleter) {
    if (!$1) {
        deleter.reset (new std::vector<long>(0));
        $1 = deleter.get():
    } else {
        $1->resize(0);
    }

    std::vector<int> shape;

    if (casac::pyarray_check($input)) {
        casac::numpy2vector((PyArrayObject*) $input, *$1, shape);
    } else {
        if (PYTEXT_CHECK($input)) {
            $1->push_back(-1);
            PyErr_SetString(PyExc_TypeError, "argument $1_name must not be a string");
            return NULL;
        } else if (PyInt_Check($input)) {
            $1->push_back(int(PyInt_AsLong($input)));
        } else if (PyLong_Check($input)) {
            $1->push_back(PyLong_AsLong($input));
        } else if (PyFloat_Check($input)) {
            $1->push_back(PyInt_AsLong(PyNumber_Long($input)));
        } else {
            shape.push_back(PyList_Size($input));
            casac::pylist2vector($input,  *$1, shape);
        }
    }
}

%typemap(in) std::vector<long long>& (std::unique_ptr<std::vector<long long> > deleter) {
    if (!$1) {
        deleter.reset (new std::vector<long long>(0));
        $1 = deleter.get();
    } else {
        $1->resize(0);
    }

    std::vector<int> shape;

    if (casac::pyarray_check($input)) {
        casac::numpy2vector((PyArrayObject*) $input, *$1, shape);
    } else {
        if (PYTEXT_CHECK($input)) {
            $1->push_back(-1);
            PyErr_SetString(PyExc_TypeError, "argument $1_name must not be a string");
            return NULL;
        } else if (PyInt_Check($input)) {
            $1->push_back(int(PyInt_AsLong($input)));
        } else if (PyLong_Check($input)) {
            $1->push_back(PyLong_AsLong($input));
        } else if (PyFloat_Check($input)) {
            $1->push_back(PyInt_AsLong(PyNumber_Long($input)));
        } else {
            shape.push_back(PyList_Size($input));
            casac::pylist2vector($input,  *$1, shape);
        }
    }
}

%typemap(in) std::vector<double>& (std::unique_ptr<std::vector<double> > deleter) {
    if (!$1) {
        deleter.reset (new std::vector<double>(0));
        $1 = deleter.get();
    } else {
        $1->resize(0);
    }

    std::vector<int> shape;

    if (casac::pyarray_check($input)) {
        casac::numpy2vector((PyArrayObject*) $input, *$1, shape);
    } else {
        if (PYTEXT_CHECK($input)) {
            $1->push_back(-1);
        } else if (PyInt_Check($input)) {
            $1->push_back(double(PyInt_AsLong($input)));
        } else if (PyLong_Check($input)) {
            $1->push_back(PyLong_AsDouble($input));
        } else if (PyFloat_Check($input)) {
            $1->push_back(PyFloat_AsDouble($input));
        } else {
            shape.push_back(PyList_Size($input));
            casac::pylist2vector($input,  *$1, shape);
        }
    }
}

%typemap(out) std::vector<bool> {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<bool>& {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<int> {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<int>& {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<long long> {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<long long>& {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<double> {
    $result = casac::map_vector($1);
}

%typemap(out) std::vector<double>& {
    $result = casac::map_vector($1);
}

%typemap(argout) std::vector<int>& OUTARGVEC {
    PyObject *o = casac::map_vector(*$1);

    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}

%typemap(argout) std::vector<double>& OUTARGVEC {
    PyObject *o= casac::map_vector(*$1);
    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}


// CASA integer vector types

%typemap(in) BoolVec (std::unique_ptr<BoolVec> deleter) {
    deleter.reset (new casac::BoolAry);
    $1 = deleter.get();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) BoolVec& (std::unique_ptr<BoolVec> deleter) {
    deleter.reset (new casac::BoolAry);
    $1 = deleter.get();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) IntVec (std::unique_ptr<IntVec> deleter) {
    deleter.reset (new casac::IntAry);
    $1 = deleter.get();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) IntVec& (std::unique_ptr<IntVec> deleter) {
    deleter.reset (new casac::IntAry);
    $1 = deleter.get();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) DoubleVec (std::unique_ptr<DoubleVec> deleter) {
    deleter.reset (new casac::DoubleAry);
    $1 = deleter.reset();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) DoubleVec& (std::unique_ptr<DoubleVec> deleter) {
    if (!$1) {
        deleter.reset (new casac::DoubleAry);
        $1 = deleter.reset();
    }
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) ComplexVec (std::unique_ptr<ComplexVec> deleter) {
    deleter.reset (new casac::ComplexAry);
    $1 = deleter.get();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(in) ComplexVec& (std::unique_ptr<ComplexVec> deleter) {
    deleter.reset (new casac::ComplexAry);
    $1 = deleter.get();
    if (pyarray_check($input)) {
        numpy2vector((PyArrayObject*) $input, $1->value, $1->shape);
    } else {
        shape.push_back(PyList_Size($input));
        pylist2vector($input,  $1->value, $1->shape);
    }
}

%typemap(out) BoolVec {
    $result = ::map_array($1.value, $1.shape);
}

%typemap(out) IntVec {
    $result = ::map_array($1.value, $1.shape);
}

%typemap(out) ComplexVec {
    $result = ::map_array($1.value, $1.shape);
}

%typemap(out) DoubleVec {
    $result = ::map_array($1.value, $1.shape);
}


// Strings

%typemap(typecheck) string {
    $1 = PYTEXT_CHECK($input) || PYBYTES_CHECK($input);
}

%typemap(typecheck) string& {
    $1 = PYTEXT_CHECK($input) || PYBYTES_CHECK($input);
}

%typemap(typecheck) const string& {
    $1 = PYTEXT_CHECK($input) || PYBYTES_CHECK($input);
}

%typemap(in) string {
    if (PYTEXT_CHECK($input)) {
        $1 = string(PYTEXT_ASDATA($input));
    } else if (PYBYTES_CHECK($input)) {
        $1 = string(PYBYTES_ASDATA($input));
    } else {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be a string");
        return NULL;
    }
}

%typemap(in) string& (std::unique_ptr<string> deleter) {
    if (PYTEXT_CHECK($input)) {
        if (!$1) {
            deleter.reset (new string(PYTEXT_ASDATA($input)));
            $1 = deleter.get();
        } else {
            *$1 = string(PYTEXT_ASDATA($input));
        }
    } else if (PYBYTES_CHECK($input)) {
        if (!$1) {
            deleter.reset (new string(PYBYTES_ASDATA($input)));
            $1 = deleter.get();
        } else {
            *$1 = string(PYBYTES_ASDATA($input));
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be a string");
        return NULL;
    }
}

%typemap(in) const string& (std::unique_ptr<string> deleter) {
    if (PYTEXT_CHECK($input)) {
        if (!$1) {
            deleter.reset (new string(PYTEXT_ASDATA($input)));
            $1 = deleter.get();
        } else {
            *$1 = string(PYTEXT_ASDATA($input));
        }
    } else if (PYBYTES_CHECK($input)) {
        if (!$1) {
            deleter.reset (new string(PYBYTES_ASDATA($input)));
            $1 = deleter.get();
        } else {
            *$1 = string(PYBYTES_ASDATA($input));
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "argument $1_name must be a string");
        return NULL;
    }
}

%typemap(freearg) const string& columnname {
    if ($1) {
        delete $1;
    }
}

%typemap(out) string {
    $result = PYSTRING_FROM_STRING($1.c_str());
}

%typemap(out) string* {
    if ($1)
        $result = PYSTRING_FROM_STRING($1->c_str());
    else
        $result = Py_None;
    delete $1;
}

%typemap(argout) string& OUTARGSTR {
    PyObject *o = PYSTRING_FROM_STRING($1->c_str());
    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}

%typemap(argout) std::string& OUTARGSTR {
    PyObject *o = PYSTRING_FROM_STRING($1->c_str());
    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}


// String vectors

%typemap(in) std::vector<std::string>& (std::unique_ptr<std::vector<std::string> > deleter) {
    if (PyList_Check($input)) {
        Py_ssize_t size = PyList_Size($input);
        if (!$1) {
            deleter.reset (new std::vector<std::string>(size));
            $1 = deleter.get();
        }
        for (Py_ssize_t i = 0; i < size; i++) {
            PyObject *o = PyList_GetItem($input,i);
            if (PYTEXT_CHECK(o))
                if (i < (Py_ssize_t)($1->size()))
                    (*$1)[i] = PYTEXT_ASDATA(PyList_GetItem($input,i));
                else
                    $1->push_back(PYTEXT_ASDATA(PyList_GetItem($input,i)));
            else {
                PyErr_SetString(PyExc_TypeError, "list $1_name must contain strings");
                return NULL;
            }
        }
    } else {
        if (PYTEXT_CHECK($input)) {
            if (!$1) {
                deleter.reset (new std::vector<std::string>(1));
                $1 = deleter.get();
            }
            if (!$1->size())
                $1->push_back(PYTEXT_ASDATA($input));
            else
                (*$1)[0] = PYTEXT_ASDATA($input);
        } else {
            PyErr_SetString(PyExc_TypeError, "$1_name is not a list");
            return NULL;
        }
    }
}

%typemap(in) StringVec {
    if (PyList_Check($input)) {
        Py_ssize_t size = PyList_Size($input);
        for (Py_ssize_t i = 0; i < size; i++) {
            PyObject *o = PyList_GetItem($input,i);
            if (PYTEXT_CHECK(o))
                $1.value.push_back(PYTEXT_ASDATA(PyList_GetItem($input,i)));
            else {
                PyErr_SetString(PyExc_TypeError, "list $1_name must contain strings");
                return NULL;
            }
        }
    } else {
        if (PYTEXT_CHECK($input)) {
            $1.value.push_back(PYTEXT_ASDATA($input));
        } else {
            PyErr_SetString(PyExc_TypeError, "$1_name is not a list");
            return NULL;
        }
    }
}

%typemap(out) std::vector<std::string> {
    $result = PyList_New($1.size());

    for(std::vector<std::string>::size_type i=0;i<$1.size();i++)
        PyList_SetItem($result, i, PYSTRING_FROM_STRING($1[i].c_str()));
}

%typemap(out) StringVec {
    $result = PyList_New($1.size());

    for(StringVec::size_type i=0;i<$1.size();i++)
        PyList_SetItem($result, i, PYSTRING_FROM_STRING($1[i].c_str()));
}

%typemap(argout) std::vector<std::string>& OUTARGVEC {
    PyObject *o = PyList_New($1.size());

    for(std::vector<std::string>::size_type i=0;i<$1.size();i++)
        PyList_SetItem($result, i, PYSTRING_FROM_STRING($1[i].c_str()));

    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}

// CASA "variant" type

%typemap(in) variant {
    if ($1) {
        (* $1) = variant(pyobj2variant($input, true));
    } else {
        PyErr_SetString (PyExc_RuntimeError, "BugCheck: Argument not initialized???");
        return nullptr;
    }
}

%typemap(in) variant& (std::unique_ptr<variant> deleter) {
    deleter.reset (new variant (pyobj2variant ($input, true)));
    $1 = deleter.get ();
}

%typemap(in) variant* (std::unique_ptr<variant> deleter) {
    deleter.reset (new variant (pyobj2variant ($input, true)));
    $1 = deleter.get ();
}

%typemap(out) variant {
    $result = variant2pyobj($1);
}

%typemap(out) variant& {
    $result = variant2pyobj($1);
}

%typemap(out) variant* {
    if ($1) {
        $result = variant2pyobj(*$1);
    } else {
        variant temp_v;
        $result = variant2pyobj(temp_v);
    }
    delete $1;
}


// CASA "Quantity" type

%typemap(in) Quantity {
    if (PyDict_Check($input)) {
        PyObject *theUnits = PyDict_GetItemString($input, "unit");
        PyObject *theVal = PyDict_GetItemString($input, "value");
        if ( theUnits && theVal) {
            std::vector<int> shape;
            std::vector<double> myVals;
            if (casac::pyarray_check(theVal)) {
                casac::numpy2vector((PyArrayObject*)theVal, myVals, shape);
            } else {
                if (PYTEXT_CHECK(theVal)) {
                    myVals.push_back(-1);
                } else if (PyInt_Check(theVal)) {
                    myVals.push_back(double(PyInt_AsLong(theVal)));
                } else if (PyLong_Check(theVal)) {
                    myVals.push_back(PyLong_AsDouble(theVal));
                } else if (PyFloat_Check(theVal)) {
                    myVals.push_back(PyFloat_AsDouble(theVal));
                } else {
                    shape.push_back(PyList_Size(theVal));
                    casac::pylist2vector(theVal,  myVals, shape);
                }
            }
            $1 = Quantity(myVals, PYTEXT_ASDATA(theUnits));
        }
    } else if (PYTEXT_CHECK($input)) {
        std::string inpstring(PYTEXT_ASDATA($input));
        double val;
        std::string units;
        istringstream iss(inpstring);
        iss >> val >> units;
        myVals.push_back(val);
        $1 = Quantity(myVals,units.c_str());
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name is not a dictionary Dictionary");
        return NULL;
    }
}

%typemap(in) Quantity * (std::unique_ptr<Quantity> deleter) {
    if (PyDict_Check($input)) {
        PyObject *theUnits = PyDict_GetItemString($input, "unit");
        PyObject *theVal = PyDict_GetItemString($input, "value");
        if ( theUnits && theVal) {
            std::vector<int> shape;
            std::vector<double> myVals;
            if (casac::pyarray_check(theVal)) {
                casac::numpy2vector((PyArrayObject*)theVal, myVals, shape);
            } else {
                if (PYTEXT_CHECK(theVal)) {
                    myVals.push_back(-1);
                } else if (PyInt_Check(theVal)) {
                    myVals.push_back(double(PyInt_AsLong(theVal)));
                } else if (PyLong_Check(theVal)) {
                    myVals.push_back(PyLong_AsDouble(theVal));
                } else if (PyFloat_Check(theVal)) {
                    myVals.push_back(PyFloat_AsDouble(theVal));
                } else {
                    shape.push_back(PyList_Size(theVal));
                    casac::pylist2vector(theVal,  myVals, shape);
                }
            }
            $1 = new Quantity(myVals,PYTEXT_ASDATA(theUnits));
        }
    } else if (PYTEXT_CHECK($input)) {
        std::string inpstring(PYTEXT_ASDATA($input));
        double val;
        std::string units;
        istringstream iss(inpstring);
        iss >> val >> units;
        myVals.push_back(val);
        deleter.reset (new Quantity(myVals,units.c_str()));
        $1 = deleter.get();
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name is not a dictionary");
        return NULL;
    }
}

%typemap(in) Quantity& (std::unique_ptr<Quantity> deleter) {
    if (PyDict_Check($input)) {
        PyObject *theUnits = PyDict_GetItemString($input, "unit");
        PyObject *theVal = PyDict_GetItemString($input, "value");
        if ( theUnits && theVal) {
            std::vector<int> shape;
            std::vector<double> myVals;
            if (casac::pyarray_check(theVal)) {
                casac::numpy2vector((PyArrayObject*)theVal, myVals, shape);
            } else {
                if (PYTEXT_CHECK(theVal)) {
                    myVals.push_back(-1);
                } else if (PyInt_Check(theVal)) {
                    myVals.push_back(double(PyInt_AsLong(theVal)));
                } else if (PyLong_Check(theVal)) {
                    myVals.push_back(PyLong_AsDouble(theVal));
                } else if (PyFloat_Check(theVal)) {
                    myVals.push_back(PyFloat_AsDouble(theVal));
                } else {
                    shape.push_back(PyList_Size(theVal));
                    casac::pylist2vector(theVal,  myVals, shape);
                }
            }
            deleter.reset (new Quantity(myVals, PYTEXT_ASDATA(theUnits)));
            $1 = deleter.get();
        }
    } else if (PYTEXT_CHECK($input)) {
        std::vector<double> myVals;
        std::string inpstring(PYTEXT_ASDATA($input));
        double val;
        std::string units;
        istringstream iss(inpstring);
        iss >> val >> units;
        myVals.push_back(val);
        deleter.reset (new Quantity(myVals,units.c_str()));
        $1 = deleter.get();
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name is not a dictionary");
        return NULL;
    }
}

%typemap(out) Quantity {
    $result = PyDict_New();
    PyDict_SetItem($result, PYSTRING_FROM_STRING("unit"), PYSTRING_FROM_STRING($1.units.c_str()));
    PyObject *v = casac::map_vector($1.value);
    PyDict_SetItem($result, PYSTRING_FROM_STRING("value"), v);
    Py_DECREF(v);
}

%typemap(out) Quantity& {
    $result = PyDict_New();
    PyDict_SetItem($result, PYSTRING_FROM_STRING("unit"), PYSTRING_FROM_STRING($1.units.c_str()));
    PyObject *v = casac::map_vector($1.value);
    PyDict_SetItem($result, PYSTRING_FROM_STRING("value"), v);
    Py_DECREF(v);
}

%typemap(out) Quantity* {
    $result = PyDict_New();
    PyDict_SetItem($result, PYSTRING_FROM_STRING("unit"), PYSTRING_FROM_STRING($1->units.c_str()));
    PyObject *v = casac::map_vector($1->value);
    PyDict_SetItem($result, PYSTRING_FROM_STRING("value"), v);
    Py_DECREF(v);
    delete $1;
}

%typemap(argout) Quantity& OUTARGQUANTITY {
    PyObject *o = PyDict_New();
    PyDict_SetItem(o, PYSTRING_FROM_STRING("unit"), PYSTRING_FROM_STRING($1->units.c_str()));
    PyObject *v = casac::map_vector($1->value);
    PyDict_SetItem(o, PYSTRING_FROM_STRING("value"), v);
    Py_DECREF(v);
    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}


// CASA "record" type

%typemap(typecheck) record& {
    if ($input)
        $1 = PyDict_Check($input);
    else
        $1 = 1;
}

%typemap(in) record {
    if (PyDict_Check($input)) {
        $1 = pyobj2variant($input, true).asRecord();
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name is not a dictionary");
        return NULL;
    }
}

%typemap(in) record * (std::unique_ptr<record> deleter) {
    if (PyDict_Check($input)) {
        deleter.reset (new record(pyobj2variant($input, true).asRecord()));
        $1 = deleter.get();
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name is not a dictionary");
        return NULL;
    }
}

%typemap(in) record& (std::unique_ptr<record> deleter) {
    if (PyDict_Check($input)) {
        deleter.reset (new record(pyobj2variant($input, true).asRecord()));
        $1 = deleter.get();
    } else {
        PyErr_SetString(PyExc_TypeError, "$1_name is not a dictionary");
        return NULL;
    }
}

%typemap(out) record {
    $result = PyDict_New();
    for(record::const_iterator iter = $1.begin(); iter != $1.end(); ++iter) {
        const std::string &key = (*iter).first;
        const casac::variant &val = (*iter).second;
        PyObject *v = casac::variant2pyobj(val);
        PyDict_SetItem($result, PYSTRING_FROM_STRING(key.c_str()), v);
        Py_DECREF(v);
    }
}

%typemap(out) record& {
    $result = PyDict_New();
    for(record::const_iterator iter = $1->begin(); iter != $1->end(); ++iter) {
        const std::string &key = (*iter).first;
        const casac::variant &val = (*iter).second;
        PyObject *v = casac::variant2pyobj(val);
        PyDict_SetItem($result, PYSTRING_FROM_STRING(key.c_str()), v);
        Py_DECREF(v);
    }
}

%typemap(out) record* {
    $result = PyDict_New();
    if ($1) {
        for(record::const_iterator iter = $1->begin(); iter != $1->end(); ++iter) {
            const std::string &key = (*iter).first;
            const casac::variant &val = (*iter).second;
            PyObject *v = casac::variant2pyobj(val);
            PyDict_SetItem($result, PYSTRING_FROM_STRING(key.c_str()), v);
            Py_DECREF(v);
        }
        delete $1;
    }
}

%typemap(argout) record& OUTARGREC {
    PyObject *o = PyDict_New();
    for(record::const_iterator iter = $1->begin(); iter != $1->end(); ++iter) {
        const std::string &key = (*iter).first;
        const casac::variant &val = (*iter).second;
        PyObject *v = casac::variant2pyobj(val);
        PyDict_SetItem(o, PYSTRING_FROM_STRING(key.c_str()), v);
        Py_DECREF(v);
    }

    if ((!$result) || ($result == Py_None)) {
        $result = o;
    } else {
        PyObject *o2 = $result;
        if (!PyTuple_Check($result)) {
            $result = PyTuple_New(1);
            PyTuple_SetItem($result,0,o2);
        }
        PyObject *o3 = PyTuple_New(1);
        PyTuple_SetItem(o3,0,o);
        o2 = $result;
        $result = PySequence_Concat(o2,o3);
        Py_DECREF(o2);
        Py_DECREF(o3);
    }
}

// CASA "RecordVec" type

%typemap(out) RecordVec {
    $result = PyList_New($1.size());
    for(RecordVec::size_type i=0;i<$1.size();i++) {
        const record &val = $1[i];
        PyObject *r = record2pydict(val);
        PyList_SetItem($result, i, r);
    }
}
