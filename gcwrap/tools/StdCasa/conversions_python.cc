#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <tools/casaswig_types.h>
#include <tools/swigconvert_python.h>
#include <Python.h>
#include <string.h>
#if USING_NUMPY_ARRAYS
#include <numpy/arrayobject.h>
#endif
#include <patchlevel.h>

#include <stdcasa/StdCasa/string_conversions>
#include <iostream>

#if PY_MAJOR_VERSION >= 3
# define PYINT_CHECK(x) 0
# define PY_INTEGER_FROM_LONG(x) PyLong_FromLong(x)
# define PYTEXT_CHECK PyUnicode_Check
# define PYTEXT_ASDATA PyUnicode_AsUTF8
#else
# define PYINT_CHECK(x) PyInt_Check(x)
# define PY_INTEGER_FROM_LONG(x) PyInt_FromLong(x)
# define PYTEXT_CHECK PyString_Check
# define PYTEXT_ASDATA PyString_AsString
#endif

STRINGTOCOMPLEX_DEFINITION(casac::complex,stringtoccomplex)

#if ( (PY_MAJOR_VERSION <= 1) || (PY_MINOR_VERSION <= 4) )
    #define MYPYSIZE int
#else
    #define MYPYSIZE Py_ssize_t
#endif

#define NODOCOMPLEX(x) (x)
#define DOCOMPLEX(x) casac::complex(x,0)
#define DOCOMPLEXCOND(cond) casac::complex((cond ? 1 : 0),0)

#if USING_NUMPY_ARRAYS

static PyArray_Descr *get_string_description(unsigned int size) {
    PyArray_Descr *string_descr = (PyArray_Descr*) malloc(sizeof(PyArray_Descr));
    memcpy(string_descr,PyArray_DescrFromType(NPY_STRING),sizeof(PyArray_Descr));
    string_descr->elsize = size;
    return string_descr;
}

#define npycomplextostring(TYPE,FORMAT)				\
inline std::string complextostring( const TYPE &c ) {		\
    char buff[512];						\
    sprintf( buff, FORMAT, c.real, c.imag );			\
    return std::string( (const char *) buff );			\
}

npycomplextostring(npy_cfloat,"(%f,%f)")
npycomplextostring(npy_cdouble,"(%f,%f)")
npycomplextostring(npy_clongdouble,"(%Lf,%Lf)")

#define CPXREALPART(cpx) cpx.real
#define CPXNONZERO(cpx) (cpx.real != 0 && cpx.imag != 0) ? true : false
#define CPXTOCCOMPLEX(cpx) casac::complex(cpx.real,cpx.imag)
#define CPXTOCOMPLEX(cpx) std::complex<double>(cpx.real,cpx.imag)

#define NUMPY2VECTOR_PLACEIT(TYPE,DATATYPE,BOOLVALUES,CPXFUNC,CPXFUNC2,STRFUNC) \
    if (PyArray_CHKFLAGS(obj,NPY_ARRAY_F_CONTIGUOUS)) {			\
        DATATYPE *data = (DATATYPE*) PyArray_DATA(obj);			\
	for ( unsigned int D=0; D < size; ++D ) {			\
	    vec[D] = (TYPE) STRFUNC CPXFUNC CPXFUNC2 (data[D] BOOLVALUES); \
	}								\
    } else {								\
        DATATYPE *data = (DATATYPE*) PyArray_DATA(obj);			\
	/*** convert c-style to fortran order ***/			\
	npy_intp counts[NPY_MAXDIMS];					\
	npy_intp strides[NPY_MAXDIMS];					\
	for ( int i=0; i < ndims; ++i ) {				\
	    counts[i] = 0;						\
	    strides[i] = PyArray_STRIDE(obj,i) / itemsize;		\
	}								\
	for ( unsigned int D=0; D < size; ++D ) {			\
	    int offset = 0;						\
	    for ( int i=0; i < ndims; ++i ) {				\
		offset += counts[i] * strides[i];			\
	    }								\
	    vec[D] = (TYPE) STRFUNC CPXFUNC CPXFUNC2 (data[offset] BOOLVALUES); \
	    for ( int j=0; j < ndims; ++j ) {				\
		if ( ++counts[j] < dims[j] ) break;			\
		counts[j] = 0;						\
	    }								\
	}								\
    }

#define NUMPY2VECTOR_PLACESTR(TYPE,CONVERT)				\
    {									\
    char *buf = new char[itemsize+1];					\
    buf[itemsize] = '\0';						\
    char *data = (char*) PyArray_DATA(obj);				\
    if (PyArray_CHKFLAGS(obj,NPY_ARRAY_F_CONTIGUOUS)) {			\
	for ( unsigned int D=0; D < size; ++D ) {			\
	    strncpy(buf,&data[D*itemsize],itemsize);			\
	    vec[D] =  CONVERT(buf);					\
	}								\
    } else {								\
	/*** convert c-style to fortran order ***/			\
	npy_intp counts[NPY_MAXDIMS];					\
	npy_intp strides[NPY_MAXDIMS];					\
	for ( int i=0; i < ndims; ++i ) {				\
	    counts[i] = 0;						\
	    strides[i] = PyArray_STRIDE(obj,i);				\
	}								\
	for ( unsigned int D=0; D < size; ++D ) {			\
	    int offset = 0;						\
	    for ( int i=0; i < ndims; ++i ) {				\
		offset += counts[i] * strides[i];			\
	    }								\
	    strncpy(buf,&data[offset],itemsize);			\
	    vec[D] = CONVERT(buf);					\
	    for ( int j=0; j < ndims; ++j ) {				\
		if ( ++counts[j] < dims[j] ) break;			\
		counts[j] = 0;						\
	    }								\
	}								\
    }									\
    delete [] buf;							\
    }

#define COPY_PYNSTRING(TYPE,NPYTYPE,START, END, ITR)			\
    char *buf = new char[itemsize+1];					\
    buf[itemsize] = '\0';						\
    std::vector<std::string>::iterator iter = ITR;			\
    for ( const char *ptr = START; ptr < END; ptr += itemsize ) {	\
	strncpy( buf, ptr, itemsize );					\
	*iter++ = (const char*) buf;					\
    }									\
    delete [] buf;

#define COPY_BUILTIN(TYPE,NPYTYPE, START, END, ITR)				\
    if (sizeof(TYPE) == sizeof (NPYTYPE))					\
	copy( START, END, ITR );						\
    else {									\
	NPYTYPE *from = (NPYTYPE*) data;					\
	for ( std::vector<TYPE>::iterator to = ITR; from < END; ++from, ++to ) {\
	    *to = (TYPE) *from;							\
	}									\
    }

#define COPY_COMPLEX(TYPE,NPYTYPE, START, END, ITR)				\
    if (sizeof(TYPE) == sizeof (NPYTYPE))					\
 	copy( (TYPE*) START, (TYPE*) END, ITR );				\
    else {									\
	NPYTYPE *from = data;							\
	for ( std::vector<TYPE>::iterator to = ITR; from < END; ++from, ++to ) {\
	    *to = TYPE((*from).real,(*from).imag);                              \
	}									\
    }

#define COPY_CCOMPLEX(TYPE,NPYTYPE, START, END, ITR)				\
    if (sizeof(TYPE) == sizeof (NPYTYPE))					\
 	copy( (TYPE*) START, (TYPE*) END, ITR );				\
    else {									\
	NPYTYPE *from = data;							\
	for ( std::vector<TYPE>::iterator to = ITR; from < END; ++from, ++to ) {\
	    (*to).re = (*from).real;						\
	    (*to).im = (*from).imag;						\
	}									\
    }

#define NUMPY2VECTOR(TYPE,NPYTYPE,NPY_TYPE,DIRECT_COPY_COND,BOOLCVT,COPY,STRCVT,CPXCVT,CPXIMAG,BTOSTR,ITOSTR,DTOSTR,CTOSTR,ASSIGN,INCR) \
void casac::numpy2vector( PyArrayObject *obj, std::vector<TYPE > &vec, std::vector<int> &shape ) { \
    npy_intp *dims = (npy_intp *)PyArray_DIMS(obj);					\
    int ndims = (int)PyArray_NDIM(obj);							\
    npy_intp itemsize = (npy_intp)PyArray_ITEMSIZE(obj);				\
    npy_uintp size = 1;									\
    for (int d=0; d < ndims; ++d) size *= dims[d];					\
    shape.resize(ndims);								\
    vec.resize(size);									\
											\
    copy( dims, dims+ndims, shape.begin() );						\
    if (DIRECT_COPY_COND && PyArray_CHKFLAGS(obj,NPY_ARRAY_F_CONTIGUOUS)) {		\
        NPY_TYPE *data = (NPY_TYPE*) PyArray_DATA(obj);					\
	COPY( TYPE, NPY_TYPE, data, data+size, vec.begin() );				\
	return;										\
    }											\
											\
    switch ( PyArray_TYPE(obj) ) {							\
    case NPY_BOOL:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_bool,BOOLCVT,,CPXIMAG,BTOSTR)			\
	break;										\
    case NPY_BYTE:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_byte,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_UBYTE:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_ubyte,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_SHORT:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_short,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_USHORT:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_ushort,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_INT:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_int,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_UINT:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_uint,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_LONG:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_long,,,CPXIMAG,ITOSTR)				\
	break;										\
    case NPY_ULONG:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_ulong,,,CPXIMAG, ITOSTR)				\
	break;										\
    case NPY_LONGLONG:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_longlong,,,CPXIMAG,ITOSTR)			\
	break;										\
    case NPY_ULONGLONG:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_ulonglong,,,CPXIMAG,ITOSTR)			\
	break;										\
    case NPY_FLOAT:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_float,,,CPXIMAG,DTOSTR)				\
	break;										\
    case NPY_DOUBLE:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_double,,,CPXIMAG,DTOSTR)				\
	break;										\
    case NPY_LONGDOUBLE:								\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_longdouble,,,CPXIMAG,DTOSTR)			\
	break;										\
    case NPY_CFLOAT:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_cfloat,,CPXCVT,,CTOSTR)				\
	break;										\
    case NPY_CDOUBLE:									\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_cdouble,,CPXCVT,,CTOSTR)				\
	break;										\
    case NPY_CLONGDOUBLE:								\
	NUMPY2VECTOR_PLACEIT(TYPE,npy_clongdouble,,CPXCVT,,CTOSTR) 			\
	break;										\
    case NPY_STRING:									\
	NUMPY2VECTOR_PLACESTR(TYPE,STRCVT)						\
	break;										\
    case NPY_OBJECT:									\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_OBJECT\n" );		\
	break;										\
    case NPY_CHAR:									\
	vec.resize(0);									\
	shape.resize(0);								\
 	fprintf( stderr, "cannot handle numpy arrays of: NPY_CHAR\n" ); 		\
	break;										\
    case NPY_UNICODE:									\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_UNICODE\n" );		\
	break;										\
    case NPY_VOID:									\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_VOID\n" ); 		\
	break;										\
    case NPY_NTYPES:									\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_NTYPES\n" );		\
	break;										\
    case NPY_NOTYPE:									\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_NOTYPE\n" );		\
	break;										\
    case NPY_USERDEF:									\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_USERDEF\n" );		\
	break;										\
    default:										\
	vec.resize(0);									\
	shape.resize(0);								\
	fprintf( stderr, "cannot handle numpy arrays of: %d\n", PyArray_TYPE(obj) );	\
	break;										\
    }											\
}

NUMPY2VECTOR(int,int,npy_int32,(PyArray_TYPE((PyArrayObject*)obj) == NPY_INT),? 1 : 0,COPY_BUILTIN,stringtoint,CPXREALPART,NODOCOMPLEX,,,,,*to = *from, 1)
NUMPY2VECTOR(double,double,npy_double,(PyArray_TYPE((PyArrayObject*)obj) == NPY_DOUBLE),? 1 : 0,COPY_BUILTIN,stringtodouble,CPXREALPART,NODOCOMPLEX,,,,,*to = *from, 1)
NUMPY2VECTOR(bool,bool,npy_bool,(PyArray_TYPE((PyArrayObject*)obj) == NPY_BOOL),? 1 : 0,COPY_BUILTIN,stringtobool,CPXNONZERO,NODOCOMPLEX,,,,,*to = (bool) *from,1)

// here we're depending upon correspondence between npy_cdouble and std::complex<double> structs --+
//                                                                                                 |
//                                vvvvvvvvvvvvvvvvvvvv---------------------------------------------+
//NUMPY2VECTOR(std::complex<double>,std::complex<double>,(PyArray_TYPE((PyArrayObject*)obj) == NPY_CDOUBLE),? 1 : 0,copy,stringtocomplex,CPXTOCOMPLEX,,,,)
NUMPY2VECTOR(std::complex<double>,std::complex<double>,npy_cdouble,(PyArray_TYPE((PyArrayObject*)obj) == NPY_CDOUBLE),? 1 : 0,COPY_COMPLEX,stringtocomplex,CPXTOCOMPLEX,,,,,,(*to).real() = (*from).real; (*to).imag() = (*from).imag,1)
NUMPY2VECTOR(casac::complex,casac::complex,npy_cdouble,(PyArray_TYPE((PyArrayObject*)obj) == NPY_CDOUBLE),? 1 : 0,COPY_CCOMPLEX,stringtoccomplex,CPXTOCCOMPLEX,DOCOMPLEX,,,,,(*to).re = (*from).real; (*to).im = (*from).imag,1)
NUMPY2VECTOR(std::string,char,const char,(false),,COPY_PYNSTRING,,,NODOCOMPLEX,booltostring,inttostring,doubletostring,complextostring,*to = *from,itemsize)

static unsigned int initialized_numpy_ = 0;
inline void initialize_numpy( ) {
    if ( initialized_numpy_ == 0 ) {
	++initialized_numpy_;
	import_array1();
    }
}
int casac::pyarray_check(PyObject *obj) {
    if ( initialized_numpy_ == 0 ) initialize_numpy( );
    return PyArray_Check(obj);
}
#endif


#define PYLIST2VECTOR_PLACEIT(VECTOR,INDEX,BOOLCVT,INTCVT,DOUBLECVT,COMPLEXCVT,STRINGCVT) \
{										\
    if (PyBool_Check(ele)) {							\
										\
	VECTOR[INDEX] = BOOLCVT( ele == Py_True );				\
										\
    } else if (PYINT_CHECK(ele)) {						\
										\
	long l = PyInt_AsLong(ele);						\
	/*** need range check ***/						\
	VECTOR[INDEX] = INTCVT( l );						\
										\
    } else if ( PyLong_Check(ele)) {						\
	long l_result = PyLong_AsLong(ele);					\
	PyObject *err = PyErr_Occurred();					\
	if ( ! err ) {								\
	    /*** need range check ***/						\
	    VECTOR[INDEX] = INTCVT(l_result);					\
	} else if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError)) {	\
	    PyErr_Clear();							\
	    long long ll_result = PyLong_AsLongLong(ele);			\
	    err = PyErr_Occurred();						\
	    if (err) {								\
		if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError))	\
		    throw casac::variant::overflow("tovariant(PyObject*,std::vector<bool>*,int,int)");\
		else								\
		    throw casac::variant::error("tovariant(PyObject*,std::vector<bool>*,int,int): unknown error");\
	    }									\
	    /*** need range check ***/						\
	    VECTOR[INDEX] = INTCVT(ll_result);					\
	}									\
    } else if ( PyFloat_Check(ele)) {						\
	VECTOR[INDEX] = DOUBLECVT(PyFloat_AsDouble(ele));			\
    } else if ( PyComplex_Check(ele)) {						\
	Py_complex c = PyComplex_AsCComplex(ele);				\
	VECTOR[INDEX] = COMPLEXCVT(c); 						\
    } else if (PYTEXT_CHECK(ele)) {						\
	VECTOR[INDEX] = STRINGCVT(PYTEXT_ASDATA(ele));			\
    } else { \
	    if (PyNumber_Check(ele)){ \
		    if(!strncmp(ele->ob_type->tp_name, "numpy.int", 9)){ \
		       VECTOR[INDEX] = INTCVT(PyLong_AsLong(PyNumber_Long(ele)));  		\
		    }else if(!strncmp(ele->ob_type->tp_name, "numpy.float", 11)){ \
		       VECTOR[INDEX] = DOUBLECVT(PyFloat_AsDouble(PyNumber_Float(ele)));  		\
		    } \
            }										\
    }										\
}										\

#define BOOLCVT_INT(cond) (cond ? 1 : 0)
#define BOOLCVT_BOOL(cond) (cond ? true : false)
#define INTCVT_BOOL(cond) (cond == 0 ? false : true)
#define CPXCVT_BOOL(cond) ((cond.real == 0 && cond.imag == 0) ? false : true)
#define CPXCVT_CPX(val) (std::complex<double>(val.real,val.imag))
#define CPXCVT_CCPX(val) (casac::complex(val.real,val.imag))
#define CVT_PASS(val) (val)
#define CPXCVT_REAL(val) (val.real)

#define PYLIST2VECTOR(TYPE,BOOLCVT,INTCVT,DOUBLECVT,COMPLEXCVT,STRINGCVT)		\
int casac::pylist2vector( PyObject *array, std::vector<TYPE> &vec, std::vector<int> &shape, int stride, int offset ) { \
											\
    if ( PyList_Check(array) || PyTuple_Check(array) ) {				\
	int number_elements = -1;							\
	bool list_elements = false;							\
	bool singleton_elements = false;						\
	bool array_is_tuple = PyTuple_Check(array) ? true : false;			\
	bool initialized_shape = false;							\
	int array_size = array_is_tuple ? PyTuple_Size(array) : PyList_Size(array);	\
											\
	int plane_offset = 1;								\
	for (int i=0; i < array_size; ++i) {						\
	    PyObject *ele = array_is_tuple ? PyTuple_GetItem(array,i) : PyList_GetItem(array,i);\
	    if ( PyList_Check(ele) || PyTuple_Check(ele) ) {				\
		bool is_tuple = PyTuple_Check(ele) ? true : false;			\
		int element_size = is_tuple ? PyTuple_Size(ele) : PyList_Size(ele);	\
											\
	      	if ( singleton_elements )						\
		    return 0;								\
		list_elements = true;							\
		if ( number_elements < 0 )						\
		    number_elements = element_size;					\
		if ( element_size != number_elements )					\
		    return 0;								\
											\
		std::vector<int> element_shape;						\
		int result = casac::pylist2vector( ele, vec, element_shape, stride*array_size, i*stride+offset ); \
											\
		if ( result == 0 ) return 0;						\
											\
		if ( element_shape.size() > 1 && plane_offset == 1 ) {			\
		    plane_offset = 1;							\
		    for ( unsigned int x=0; x < element_shape.size(); ++x )		\
		      plane_offset *= element_shape[x];					\
		}									\
											\
		/*** handle shape of the sub-array ***/					\
		if ( ! initialized_shape ) {						\
		    shape.resize( element_shape.size() + 1 );				\
		    for ( unsigned int i = 0; i < element_shape.size(); ++i )		\
			shape[i+1] = element_shape[i];					\
		    shape[0] = array_size;						\
		} else {								\
		    if ( shape.size() != element_shape.size() + 1 )			\
			return 0;							\
		    if ( element_shape.size() == 1 ) {					\
			if ( element_shape[0] != shape[1] )				\
			    return 0;							\
		    } else {								\
			for ( unsigned int i = 0; i < element_shape.size(); ++i )	\
			    if ( element_shape[i] != shape[i] )				\
				return 0;						\
		    }									\
		}									\
											\
	    } else {									\
											\
		if ( list_elements )							\
		    return 0;								\
		if ( singleton_elements == false ) {					\
		    singleton_elements = true;						\
		    shape.resize(1);							\
		    shape[0] =  array_size;						\
											\
		    int bulk_size = (array_size-1) * stride + 1;			\
											\
		    if (vec.size() < (unsigned) (offset + bulk_size))			\
			vec.resize(offset + bulk_size);					\
		}									\
											\
		PYLIST2VECTOR_PLACEIT(vec,offset+(i*stride),BOOLCVT,INTCVT,DOUBLECVT,COMPLEXCVT,STRINGCVT); \
	    }										\
	}										\
											\
	return 1;									\
    }											\
											\
    return 0;										\
}

PYLIST2VECTOR(int,BOOLCVT_INT,CVT_PASS,(int) CVT_PASS, (int) CPXCVT_REAL, stringtoint)
PYLIST2VECTOR(double,(double) BOOLCVT_INT,(double) CVT_PASS,CVT_PASS,CPXCVT_REAL, stringtodouble)
PYLIST2VECTOR(bool,BOOLCVT_BOOL,INTCVT_BOOL,INTCVT_BOOL,CPXCVT_BOOL, stringtobool)
PYLIST2VECTOR(std::complex<double>,(double) BOOLCVT_INT,(double) CVT_PASS,CVT_PASS,CPXCVT_CPX,stringtocomplex)
PYLIST2VECTOR(casac::complex,DOCOMPLEXCOND,DOCOMPLEX,DOCOMPLEX,CPXCVT_CCPX,stringtoccomplex)
PYLIST2VECTOR(std::string,booltostring BOOLCVT_BOOL,inttostring,doubletostring,complextostring CPXCVT_CPX,CVT_PASS)

namespace casac {
int convert_idl_complex_from_python_complex(PyObject *obj,void *s) {
    if ( PyComplex_Check(obj) ) {
	casac::complex *to = (casac::complex*) s;
	to->re = PyComplex_RealAsDouble(obj);
	to->im = PyComplex_ImagAsDouble(obj);
	return 1;
    }

    PyErr_SetString( PyExc_TypeError, "not a complex" );
    return 0;
}
PyObject *convert_idl_complex_to_python_complex(const casac::complex &from) {
    return PyComplex_FromDoubles(from.re,from.im);
}

#define RECORD2PYDICT													\
    for ( record::const_iterator iter = rec.begin(); iter != rec.end(); ++iter ) {					\
	const std::string &key = (*iter).first;										\
	const variant &val = (*iter).second;										\
	PyObject *v = variant2pyobj( val );										\
	PyDict_SetItem(result, PyUnicode_FromString(key.c_str()), v);							\
	Py_DECREF(v);													\
    }															\
															\
    return result;

PyObject *record2pydict(const record &rec) {
    PyObject *result = PyDict_New();
    RECORD2PYDICT
}



#define PLACEIT(VARIANT,INDEX)							\
{										\
    if (PyBool_Check(ele)) {							\
										\
	VARIANT.place( ele == Py_True ? true : false, INDEX );			\
										\
    } else if (PYINT_CHECK(ele)) {						\
										\
	long l = PyInt_AsLong(ele);						\
	/*** need range check ***/						\
	VARIANT.place((int)l,INDEX);						\
										\
    } else if ( PyLong_Check(ele)) {						\
	long l_result = PyLong_AsLong(ele);					\
	PyObject *err = PyErr_Occurred();					\
	if ( ! err ) {								\
	    /*** need range check ***/						\
	    VARIANT.place((int)l_result, INDEX);				\
	} else if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError)) {	\
	    PyErr_Clear();							\
	    long long ll_result = PyLong_AsLongLong(ele);			\
	    err = PyErr_Occurred();						\
	    if (err) {								\
		if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError))	\
		    throw casac::variant::overflow("tovariant(PyObject*,std::vector<bool>*,int,int)");\
		else								\
		    throw casac::variant::error("tovariant(PyObject*,std::vector<bool>*,int,int): unknown error");\
	    }									\
	    /*** need range check ***/						\
	    VARIANT.place((int)ll_result,INDEX);				\
	}									\
    } else if ( PyFloat_Check(ele)) {						\
	VARIANT.place(PyFloat_AsDouble(ele),INDEX);				\
    } else if ( PyComplex_Check(ele)) {						\
      static int M = 0; ++M; \
	Py_complex c = PyComplex_AsCComplex(ele);				\
	VARIANT.place(std::complex<double>(c.real, c.imag),INDEX);		\
    } else if ( PyNumber_Check(ele)) {						\
	if(!strncmp(ele->ob_type->tp_name, "numpy.int", 9)){ \
	   VARIANT.place((int)PyLong_AsLong(PyNumber_Long(ele)),INDEX);  		\
	}else if(!strncmp(ele->ob_type->tp_name, "numpy.float", 11)){ \
	   VARIANT.place(double(PyFloat_AsDouble(PyNumber_Float(ele))),INDEX);  		\
	} \
    } else if (PYTEXT_CHECK(ele)) {						\
	VARIANT.place(std::string(PYTEXT_ASDATA(ele)),INDEX);		\
    }										\
}										\


#define FWD_DECL_map_array_pylist( TYPE )							\
static PyObject *map_array_pylist( const std::vector<TYPE> &vec, const std::vector<int> &shape, \
				   int stride=1, int offset=0 );				\
static PyObject *map_vector_pylist( const std::vector<TYPE> &vec );

FWD_DECL_map_array_pylist(int)
FWD_DECL_map_array_pylist(unsigned int)
FWD_DECL_map_array_pylist(bool)
FWD_DECL_map_array_pylist(double)
FWD_DECL_map_array_pylist(long long)
FWD_DECL_map_array_pylist(std::complex<double>)
FWD_DECL_map_array_pylist(casac::complex)
FWD_DECL_map_array_pylist(std::string)

#if USING_NUMPY_ARRAYS
//static int unmap_array_numpy( PyObject *array, std::vector<int> &shape, casac::variant &vnt ) {

#define MAP_ARRAY_NUMPY(TYPE,NPYTYPE,NUMPY_TYPE,ASSIGN)						\
PyObject *map_vector_numpy(const std::vector<TYPE> &vec) {				\
    initialize_numpy( );									\
    PyArray_Descr *type = PyArray_DescrFromType(NUMPY_TYPE);					\
    npy_intp dim[1];										\
    dim[0] = vec.size();									\
    PyObject *ary =  PyArray_NewFromDescr( &PyArray_Type, type, 1, (npy_intp*) &dim,		\
					   NULL /*npy_intp *strides*/,				\
					   NULL /*char *data*/,					\
					   1 /* zero=>c-order, nonzero=>fortran */,		\
					   NULL /*PyObject *obj*/ );				\
    if (sizeof(TYPE) == sizeof (NPYTYPE))							\
	copy( vec.begin(), vec.end(), (TYPE*) PyArray_DATA((PyArrayObject*)ary) );		\
    else {											\
	std::vector<TYPE>::const_iterator from = vec.begin();					\
	for ( NPYTYPE *to = (NPYTYPE*) PyArray_DATA((PyArrayObject*)ary); from != vec.end();    \
              ++from, ++to ) {	                                                                \
	    ASSIGN;										\
	}											\
    }												\
    return ary;											\
}												\
												\
PyObject *map_array_numpy( const std::vector<TYPE> &vec, const std::vector<int> &shape ) { \
    initialize_numpy( );									\
    PyArray_Descr *type = PyArray_DescrFromType(NUMPY_TYPE);					\
    npy_intp dim[NPY_MAXDIMS];                 							\
    for (size_t i = 0; i < shape.size(); i++) {                                                 \
        dim[i] = (npy_intp)shape[i];                                                            \
    }                                                                                           \
    PyObject *ary =  PyArray_NewFromDescr( &PyArray_Type, type, shape.size(), dim,		\
					   NULL /*npy_intp *strides*/,				\
					   NULL /*char *data*/,					\
					   1 /* zero=>c-order, nonzero=>fortran */,		\
					   NULL /*PyObject *obj*/ );				\
    if (sizeof(TYPE) == sizeof (NPYTYPE))							\
	copy( vec.begin(), vec.end(), (TYPE*) PyArray_DATA((PyArrayObject*)ary) );		\
    else {											\
	std::vector<TYPE>::const_iterator from = vec.begin();					\
	for ( NPYTYPE *to = (NPYTYPE*) PyArray_DATA((PyArrayObject*)ary); from != vec.end();    \
              ++from, ++to ) {	                                                                \
	    ASSIGN;										\
	}											\
    }												\
    return ary;											\
}

MAP_ARRAY_NUMPY(int, npy_int32, NPY_INT,*to = (npy_int32) *from)
MAP_ARRAY_NUMPY(unsigned int, npy_uint, NPY_UINT,*to = (npy_uint) *from)
MAP_ARRAY_NUMPY(long long, npy_int64, NPY_INT64,*to = (npy_int64) *from)
MAP_ARRAY_NUMPY(double, npy_double,NPY_DOUBLE,*to = (npy_double) *from)
MAP_ARRAY_NUMPY(std::complex<double>, npy_cdouble, NPY_CDOUBLE,(*to).real = (*from).real(); (*to).imag = (*from).imag())
MAP_ARRAY_NUMPY(casac::complex, npy_cdouble, NPY_CDOUBLE,(*to).real = (*from).re; (*to).imag = (*from).im)
MAP_ARRAY_NUMPY(bool, npy_bool, NPY_BOOL,*to = (npy_bool) *from)

PyObject *map_vector_numpy( const std::vector<std::string> &vec ) {
    initialize_numpy( );
    unsigned int size = 0;
    for ( std::vector<std::string>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter ) {
	unsigned int slen = (*iter).size( );
	if ( slen > size ) size = slen;
    }
    npy_intp dim[1];
    dim[0] = vec.size();
    PyObject *ary = PyArray_NewFromDescr( &PyArray_Type, get_string_description(size+1), 1, dim,
					  NULL, NULL, 1, NULL );
    char *data = (char*) PyArray_DATA((PyArrayObject*)ary);
    memset(data, '\0', (size+1) * vec.size());
    for ( std::vector<std::string>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter ) {
	strcpy(data,(*iter).c_str( ));
	data += size+1;
    }
    return ary;
}

PyObject *map_array_numpy( const std::vector<std::string> &vec, const std::vector<int> &shape ) {
    initialize_numpy( );
    unsigned int size = 0;
    for ( std::vector<std::string>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter ) {
	unsigned int slen = (*iter).size( );
	if ( slen > size ) size = slen;
    }
    npy_intp *dim = new npy_intp[shape.size()];
    copy( shape.begin(), shape.end(), dim );
    PyObject *ary = PyArray_NewFromDescr( &PyArray_Type, get_string_description(size+1), shape.size(), dim,
					  NULL, NULL, 1, NULL );
    char *data = (char*) PyArray_DATA((PyArrayObject*)ary);
    memset(data, '\0', (size+1) * vec.size());
    for ( std::vector<std::string>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter ) {
	strcpy(data,(*iter).c_str( ));
	data += size+1;
    }
    delete [] dim;
    return ary;
}

#define HANDLE_NUMPY_ARRAY_CASE(NPYTYPE,CVTTYPE)				\
    case NPYTYPE:								\
        casac::numpy2vector((PyArrayObject*)obj,result.as ## CVTTYPE ## Vec(),  \
                            result.arrayshape()); 	                        \
	break;

#define HANDLE_NUMPY_ARRAY							\
else if ( casac::pyarray_check(obj) ) {						\
    switch ( PyArray_TYPE((PyArrayObject*)obj) ) {				\
    HANDLE_NUMPY_ARRAY_CASE(NPY_BOOL,Bool)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_BYTE,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_UBYTE,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_SHORT,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_USHORT,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_INT,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_UINT,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_LONG,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_ULONG,Int)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_LONGLONG,Double)				\
    HANDLE_NUMPY_ARRAY_CASE(NPY_ULONGLONG,Double)				\
    HANDLE_NUMPY_ARRAY_CASE(NPY_FLOAT,Double)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_DOUBLE,Double)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_LONGDOUBLE,Double)				\
    HANDLE_NUMPY_ARRAY_CASE(NPY_CFLOAT,Complex)					\
    HANDLE_NUMPY_ARRAY_CASE(NPY_CDOUBLE,Complex)				\
    HANDLE_NUMPY_ARRAY_CASE(NPY_CLONGDOUBLE,Complex)				\
    HANDLE_NUMPY_ARRAY_CASE(NPY_STRING,String)					\
    case NPY_OBJECT:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_OBJECT\n" );	\
	break;									\
    case NPY_CHAR:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_CHAR\n" );		\
	break;									\
    case NPY_UNICODE:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_UNICODE\n" );	\
	break;									\
    case NPY_VOID:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_VOID\n" );		\
	break;									\
    case NPY_NTYPES:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_NTYPES\n" );	\
	break;									\
    case NPY_NOTYPE:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_NOTYPE\n" );	\
	break;									\
    case NPY_USERDEF:								\
	fprintf( stderr, "cannot handle numpy arrays of: NPY_USERDEF\n" );	\
	break;									\
    default:									\
	fprintf( stderr, "cannot handle numpy arrays of: %d\n", PyArray_TYPE((PyArrayObject*)obj) ); \
	break;									\
    }										\
}

#define NOT_NUMPY_ARRAY(OBJ) && ! casac::pyarray_check(OBJ)

#else

#define FWD_map_array_pylist( TYPE )								\
inline PyObject *map_array( const std::vector<TYPE> &vec, const std::vector<int> &shape ) {	\
    return map_array_pylist(vec, shape);							\
}												\
												\
inline PyObject *map_vector( const std::vector<TYPE> &vec ) {					\
    return map_vector_pylist( vec );								\
}

FWD_map_array_pylist(int)
FWD_map_array_pylist(unsigned int)
FWD_map_array_pylist(bool)
FWD_map_array_pylist(double)
FWD_map_array_pylist(std::complex<double>)
FWD_map_array_pylist(casac::complex)
FWD_map_array_pylist(std::string)

#define HANDLE_NUMPY_ARRAY
#define NOT_NUMPY_ARRAY(OBJ)

#endif

//
// returns non-zero upon success
//
static int unmap_array_pylist( PyObject *array, std::vector<int> &shape, casac::variant &vnt, int stride = 1, int offset = 0 ) {

    if ( PyList_Check(array) || PyTuple_Check(array) ) {
	int number_elements = -1;
	bool list_elements = false;
	bool singleton_elements = false;
	bool array_is_tuple = PyTuple_Check(array) ? true : false;
	bool initialized_shape = false;
	int array_size = array_is_tuple ? PyTuple_Size(array) : PyList_Size(array);

	int plane_offset = 1;
	for (int i=0; i < array_size; ++i) {
	    PyObject *ele = array_is_tuple ? PyTuple_GetItem(array,i) : PyList_GetItem(array,i);
	    if ( PyList_Check(ele) || PyTuple_Check(ele) ) {
		bool is_tuple = PyTuple_Check(ele) ? true : false;
		int element_size = is_tuple ? PyTuple_Size(ele) : PyList_Size(ele);

	      	if ( singleton_elements ){
		    return 0;
		}
		list_elements = true;
		if ( number_elements < 0 )
		    number_elements = element_size;
                // to allow irregular shaped python list 
                // (e.g. [[1,2,3],[4]]
		//if ( element_size != number_elements )
		//    return 0;

		std::vector<int> element_shape;
		int result = unmap_array_pylist( ele, element_shape, vnt, stride*array_size, i*stride+offset );

		if ( result == 0 ) return 0;

		if ( element_shape.size() > 1 && plane_offset == 1 ) {
		    plane_offset = 1;
		    for ( unsigned int x=0; x < element_shape.size(); ++x )
		      plane_offset *= element_shape[x];
		}

		//
		// handle shape of the sub-array
		//
		if ( ! initialized_shape ) {
		    shape.resize( element_shape.size() + 1 );
		    for ( unsigned int i = 0; i < element_shape.size(); ++i )
			shape[i+1] = element_shape[i];
		    shape[0] = array_size;
		} else {
		    if ( shape.size() != element_shape.size() + 1 )
			return 0;
		    if ( element_shape.size() == 1 ) {
			if ( element_shape[0] != shape[1] )
			    return 0;
		    } else {
			for ( unsigned int i = 0; i < element_shape.size(); ++i )
			    if ( element_shape[i] != shape[i] )
				return 0;
		    }
		}

	    } else {

		if ( list_elements )
		    return 0;
		if ( singleton_elements == false ) {
		    singleton_elements = true;
		    shape.resize(1);
		    shape[0] =  array_size;

		    int bulk_size = (array_size-1) * stride + 1;

		    if (vnt.size() < (offset + bulk_size)) {
			vnt.resize(offset + bulk_size);
		    }
		}

		PLACEIT(vnt,offset+(i*stride));
	    }
	}

	return 1;
    }

    return 0;
}

#define PYOBJ2VARIANT(SINGLETON,CREATION,DO_THROW)					\
    if ( PyBool_Check(obj) )								\
	SINGLETON(obj == Py_True ? true : false );					\
											\
    else if ( PYINT_CHECK(obj) )							\
        /*** need range check ***/							\
	SINGLETON((int) PyInt_AsLong(obj) );						\
											\
    else if ( PyLong_Check(obj) ) {							\
	long l_result = PyLong_AsLong(obj);						\
	PyObject *err = PyErr_Occurred();						\
	if ( ! err )									\
            /*** need range check ***/							\
	    SINGLETON((int)l_result );							\
	else if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError)) {		\
	    PyErr_Clear();								\
	    long long ll_result = PyLong_AsLongLong(obj);				\
	    err = PyErr_Occurred();							\
	    if (err) {									\
		if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError))		\
		    throw variant::overflow("convert_single(long long)");		\
		else									\
		    throw variant::error("convert_single(long long): unknown error");	\
	    }										\
            /*** need range check ***/							\
	    SINGLETON((int)ll_result );							\
	}										\
    }											\
											\
    else if ( PyFloat_Check(obj) )							\
	SINGLETON(PyFloat_AsDouble(obj) );						\
											\
    else if ( PYTEXT_CHECK(obj) NOT_NUMPY_ARRAY(obj) )				\
	SINGLETON(std::string(PYTEXT_ASDATA(obj)));					\
											\
    else if ( PyComplex_Check(obj) ) {							\
											\
	Py_complex c = PyComplex_AsCComplex(obj);					\
	SINGLETON(std::complex<double>(c.real,c.imag) );				\
											\
    }											\
    else if (PyNumber_Check(obj)) {					\
	if(!strncmp(obj->ob_type->tp_name, "numpy.int", 9)){ \
	   SINGLETON((int)PyLong_AsLong(PyNumber_Long(obj)));  		\
	}else if(!strncmp(obj->ob_type->tp_name, "numpy.float", 11)){ \
	   SINGLETON(double(PyFloat_AsDouble(PyNumber_Float(obj))));  		\
	} \
    } \
											\
    CREATION										\
											\
    if (PyList_Check(obj) || PyTuple_Check(obj)) {					\
											\
	bool is_list = PyList_Check(obj) ? true : false;				\
	bool done = false;								\
											\
	if ( ( is_list ? PyList_Size(obj) > 0 : PyTuple_Size(obj) > 0 ) &&              \
		( is_list ? PyList_Check(PyList_GetItem(obj,0)) || PyTuple_Check(PyList_GetItem(obj,0)) :       \
			    PyList_Check(PyTuple_GetItem(obj,0)) || PyTuple_Check(PyTuple_GetItem(obj,0))) ) {	\
	    std::vector<int> oshape;							\
	    int outcome = unmap_array_pylist( obj, oshape, result );			\
	    if ( outcome != 0 ) {							\
		result.arrayshape() = oshape;						\
		done = true;								\
	    } else {									\
		if ( DO_THROW )								\
		    throw std::string("all list elements do not confrom to array shape"); \
		result.asBool();							\
	    }										\
	}										\
											\
	if ( ! done ) {									\
											\
	    for ( int i=0; ( is_list ? i < PyList_Size(obj) : i < PyTuple_Size(obj) ); ++ i) {		\
		PyObject *ele = is_list ? PyList_GetItem(obj,i) : PyTuple_GetItem(obj, i);              \
											\
		if (PyBool_Check(ele)) {						\
											\
		    result.push( ele == Py_True ? true : false );			\
											\
		} else if (PYINT_CHECK(ele)) {						\
											\
		    long l = PyInt_AsLong(ele);						\
		    /*** need range check ***/						\
		    result.push((int)l);						\
											\
		} else if ( PyLong_Check(ele)) {					\
		    long l_result = PyLong_AsLong(ele);					\
		    PyObject *err = PyErr_Occurred();					\
		    if ( ! err ) {							\
			/*** need range check ***/					\
			result.push((int)l_result);					\
		    } else if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError)) {	\
			PyErr_Clear();							\
			long long ll_result = PyLong_AsLongLong(ele);			\
			err = PyErr_Occurred();						\
			if (err) {							\
			    if (PyErr_GivenExceptionMatches(err,PyExc_OverflowError))	\
				throw variant::overflow("tovariant(PyObject*,std::vector<bool>*,int,int)");\
			    else							\
				throw variant::error("tovariant(PyObject*,std::vector<bool>*,int,int): unknown error");\
			}								\
			/*** need range check ***/					\
			result.push((int)ll_result);					\
		    }									\
		} else if ( PyFloat_Check(ele)) {					\
		    result.push(PyFloat_AsDouble(ele));					\
		} else if ( PyComplex_Check(ele)) {					\
		    Py_complex c = PyComplex_AsCComplex(ele);				\
		    result.push(std::complex<double>(c.real, c.imag));			\
		} else if (PYTEXT_CHECK(ele)) {					\
		    result.push(std::string(PYTEXT_ASDATA(ele)));			\
		} else if (PyNumber_Check(ele)) {					\
		    if(!strncmp(ele->ob_type->tp_name, "numpy.int", 9)){ \
		       result.push((int)PyLong_AsLong(PyNumber_Long(ele)));  		\
		    }else if(!strncmp(ele->ob_type->tp_name, "numpy.float", 11)){ \
		       result.push(double(PyFloat_AsDouble(PyNumber_Float(ele))));  		\
		    } \
		} else if (PyList_Check(ele) || PyTuple_Check(ele)) {			\
		    pyobj2variant(ele,result);						\
		}									\
	    }										\
	}										\
    } 											\
											\
    HANDLE_NUMPY_ARRAY									\
											\
    else if (PyDict_Check(obj)) {							\
	PyObject *key, *val;								\
        MYPYSIZE pos = 0;                                                               \
	record &rec = result.asRecord( );						\
	while ( PyDict_Next(obj, &pos, &key, &val) ) {					\
	    const char *str = 0;							\
	    PyObject *strobj = 0;							\
	    if (PYTEXT_CHECK(key)) {							\
		str = PYTEXT_ASDATA(key);						\
	    } else {									\
		strobj = PyObject_Str(key);						\
		str = PYTEXT_ASDATA(strobj);					\
	    }										\
            if(PyBool_Check(val) || PYINT_CHECK(val) || PyLong_Check(val) || \
               PyFloat_Check(val) || PYTEXT_CHECK(val) || PyComplex_Check(val) ||     \
               PyList_Check(val) || PyTuple_Check(val) || PyDict_Check(val) ||          \
               (casac::pyarray_check(val) &&                                            \
                (PyArray_TYPE((PyArrayObject*)val) == NPY_BOOL ||                       \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_BYTE ||                       \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_UBYTE ||                      \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_SHORT ||                      \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_USHORT ||                     \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_INT ||                        \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_UINT ||                       \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_LONG ||                       \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_ULONG ||                      \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_LONGLONG ||                   \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_ULONGLONG||                   \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_FLOAT ||                      \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_DOUBLE ||                     \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_LONGDOUBLE ||                 \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_CFLOAT ||                     \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_CDOUBLE ||                    \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_CLONGDOUBLE ||                \
                 PyArray_TYPE((PyArrayObject*)val) == NPY_STRING))) {                   \
	      rec.insert(str,pyobj2variant(val,DO_THROW));				\
            }                                                                           \
                                                                                        \
	    Py_XDECREF(strobj);								\
	}										\
    }

void pyobj2variant(PyObject *obj, variant &result) {
    PYOBJ2VARIANT(result.push,,false)
}

#define ARRAY2PYOBJ(TYPE,ELEMENT_CTOR,RHS,INCREF,STAGE)						\
static PyObject *map_array_vector_pylist( const std::vector<TYPE> &vec, unsigned int stride=1,	\
			     unsigned int offset=0, unsigned int len= 0 ) {			\
												\
    if ( len == 0 )										\
	len = vec.size();									\
												\
    PyObject *result = PyTuple_New(len);							\
												\
    for ( unsigned int i=0; i < len && (offset+(i*stride)) < vec.size(); ++i ) {		\
	TYPE val = vec[offset+(i*stride)];							\
	PyObject *ele = ELEMENT_CTOR;								\
	PyTuple_SetItem(result, i, ele);							\
    }												\
    return result;										\
}												\
												\
static PyObject *map_vector_pylist( const std::vector<TYPE> &vec ) {				\
    int cnt = 0;										\
    PyObject *list = PyList_New(vec.size());							\
    for ( std::vector<TYPE>::const_iterator iter = vec.begin(); iter != vec.end(); ++iter ) {	\
	STAGE											\
	PyObject *vec_val = RHS;								\
	INCREF											\
	PyList_SetItem(list, cnt++, vec_val);							\
    }												\
    return list;										\
}												\
												\
static PyObject *map_array_pylist( const std::vector<TYPE> &vec, const std::vector<int> &shape, int stride, int offset ) { \
    PyObject *result = NULL;									\
    if ( shape.size() == 0 )									\
	return PyList_New(0);									\
    else if ( shape.size() == 1 )								\
	return map_array_vector_pylist( vec, 1, 0, shape[0] );					\
    else if ( shape.size() == 2 ) {								\
	result = PyList_New(shape[0]);								\
	for ( int i = 0; i < shape[0]; ++i ) {							\
	    PyList_SetItem(result,i, map_array_vector_pylist( vec, stride*shape[0],		\
						 i*stride+offset, shape[1]));			\
	}											\
    } else {											\
	result = PyList_New(shape[0]);								\
	std::vector<int> plane_shape(shape.size()-1);						\
	unsigned int plane_size = 1;								\
	for (unsigned int i=0; i < shape.size()-1; ++i) {					\
	    plane_size *= shape[i+1];								\
	    plane_shape[i] = shape[i+1];							\
	}											\
	for (int i=0; i < shape[0]; ++i) {							\
	    PyList_SetItem(result,i,map_array_pylist( vec, plane_shape,				\
					       stride*shape[0], i*stride+offset ));		\
	}											\
    }												\
												\
    return result;										\
}

ARRAY2PYOBJ(int,PY_INTEGER_FROM_LONG(val),PY_INTEGER_FROM_LONG(*iter),,)
ARRAY2PYOBJ(unsigned int,PY_INTEGER_FROM_LONG(val),PY_INTEGER_FROM_LONG(*iter),,)
ARRAY2PYOBJ(bool,(val == 0 ? Py_False : Py_True); Py_INCREF(ele),(*iter == false ? Py_False : Py_True),Py_INCREF(vec_val);,)
ARRAY2PYOBJ(double,PyFloat_FromDouble(val),PyFloat_FromDouble(*iter),,)
ARRAY2PYOBJ(std::complex<double> ,PyComplex_FromDoubles(val.real(),val.imag()),PyComplex_FromDoubles(cpx.real(),cpx.imag()),,std::complex<double> cpx = *iter;)
ARRAY2PYOBJ(casac::complex ,PyComplex_FromDoubles(val.re,val.im),PyComplex_FromDoubles(cpx.re,cpx.im),,casac::complex cpx = *iter;)
ARRAY2PYOBJ(std::string,PyUnicode_FromString(val.c_str()),PyUnicode_FromString((*iter).c_str()),,)

#define HANDLEVEC2(TYPE,FETCH)												\
{															\
	const std::vector<TYPE> &vec = val.FETCH();									\
	const std::vector<int> &shape = val.arrayshape( );								\
															\
	unsigned int shape_size = 1;											\
	if ( shape.size() > 1 ) {											\
	    for ( unsigned int i = 0; i < shape.size(); ++i ) {								\
		if ( shape[i] <= 0 ) {											\
		    shape_size = 0;											\
		    break;												\
		}													\
		shape_size *= shape[i];											\
	    }														\
	}														\
															\
	if ( shape.size() > (unsigned) 1 && shape_size == vec.size() ) {						\
	    result = map_array( vec, shape );										\
	} else {													\
	    result = map_vector( vec );											\
	}														\
    break;														\
}

#define VARIANT2PYOBJ													\
	switch (val.type()) {												\
	    case variant::RECORD:											\
		result =  record2pydict(val.getRecord());								\
		break;													\
	    case variant::BOOL:												\
		{													\
		PyObject *bool_val = val.toBool() == false ? Py_False : Py_True;					\
		Py_INCREF(bool_val);											\
		result =  bool_val;											\
		break;													\
		}													\
	    case variant::INT:												\
		result = PY_INTEGER_FROM_LONG(val.toInt());									\
		break;													\
        case variant::UINT:                                             \
        result = PY_INTEGER_FROM_LONG(val.touInt());                          \
        break;                                                          \
	    case variant::LONG:												\
		result = PY_INTEGER_FROM_LONG(val.toLong());									\
		break;													\
	    case variant::DOUBLE:											\
		result = PyFloat_FromDouble(val.toDouble());								\
		break;													\
	    case variant::COMPLEX:											\
		{													\
		const std::complex<double> &c = val.getComplex( );							\
		result = PyComplex_FromDoubles(c.real(),c.imag());							\
		break;													\
		}													\
	    case variant::STRING:											\
		result = PyUnicode_FromString(val.toString().c_str());							\
		break;													\
															\
	    case variant::BOOLVEC:											\
		HANDLEVEC2(bool,getBoolVec)										\
															\
	    case variant::INTVEC:											\
		HANDLEVEC2(int,getIntVec)										\
															\
        case variant::UINTVEC:                                          \
        HANDLEVEC2(unsigned int,getuIntVec)                             \
                                                                        \
	    case variant::LONGVEC:											\
		HANDLEVEC2(long long,getLongVec)										\
															\
	    case variant::DOUBLEVEC:											\
		HANDLEVEC2(double,getDoubleVec)										\
															\
	    case variant::COMPLEXVEC:											\
		HANDLEVEC2(std::complex<double>,getComplexVec) 								\
															\
	    case variant::STRINGVEC:											\
		HANDLEVEC2(std::string,getStringVec)									\
															\
	    default:													\
		fprintf( stderr, "encountered unknown variant type in pyobj2variant()!\n" );				\
	}

PyObject *variant2pyobj(const variant &val) {
    PyObject *result = 0;
    VARIANT2PYOBJ
    return result;
}

variant pyobj2variant(PyObject *obj, bool throw_error) {
    variant result;

    PYOBJ2VARIANT(return variant,,throw_error)

    return result;
}

#if USING_NUMPY_ARRAYS
#define ARYSTRUCT_CONVERT(ARY)								\
static int python_ ## ARY ## Ary_check_pytuple(PyObject *);				\
											\
static int python_ ## ARY ## Ary_check_pylist(PyObject *obj) {				\
    for ( int i=0; i < PyList_Size(obj); ++i ) {					\
	PyObject *t = PyList_GetItem(obj,i);						\
	if (PyList_Check(obj)) {							\
	    if (! python_ ## ARY ## Ary_check_pylist(t))				\
		return 0;								\
	} else if (PyTuple_Check(obj)) {						\
	    if (! python_ ## ARY ## Ary_check_pytuple(t))				\
		return 0;								\
	} else if ( ! PyBool_Check(t)) {						\
		return 0;								\
	}										\
    }											\
    return 1;										\
}											\
											\
static int python_ ## ARY ## Ary_check_pytuple(PyObject *obj) {				\
    for ( int i=0; i < PyTuple_Size(obj); ++i ) {					\
	PyObject *t = PyTuple_GetItem(obj,i);						\
	if (PyList_Check(obj)) {							\
	    if (! python_ ## ARY ## Ary_check_pylist(t))				\
		return 0;								\
	} else if (PyTuple_Check(obj)) {						\
	    if (! python_ ## ARY ## Ary_check_pytuple(t))				\
		return 0;								\
	} else if ( ! PyBool_Check(t)) {						\
		return 0;								\
	}										\
    }											\
    return 1;										\
}											\
											\
int python_ ## ARY ## Ary_check(PyObject *obj) {					\
    if (casac::pyarray_check(obj)) return 1;							\
    if (PyList_Check(obj)) return python_ ## ARY ## Ary_check_pylist(obj);		\
    if (PyTuple_Check(obj)) return python_ ## ARY ## Ary_check_pytuple(obj);		\
    return 0;										\
}
											\
#else
#define ARYSTRUCT_CONVERT(ARY)								\
											\
int python_ ## ARY ## Ary_check(PyObject *obj) {					\
	return PyList_Check(obj) || PyTuple_Check(obj);					\
}
}
#endif
ARYSTRUCT_CONVERT(Bool)
ARYSTRUCT_CONVERT(Int)
ARYSTRUCT_CONVERT(Double)
//ARYSTRUCT_CONVERT(Complex)
ARYSTRUCT_CONVERT(String)

int is_intvec_compatible_numpy_array( PyObject *obj ) {
    if ( casac::pyarray_check(obj) &&
	 ( PyArray_TYPE((PyArrayObject*)obj) == NPY_BOOL ||
	   PyArray_TYPE((PyArrayObject*)obj) == NPY_BYTE ||
	   PyArray_TYPE((PyArrayObject*)obj) == NPY_UBYTE ||
	   PyArray_TYPE((PyArrayObject*)obj) == NPY_SHORT ||
	   PyArray_TYPE((PyArrayObject*)obj) == NPY_INT ) ) {
	return 1;
    } else {
	return 0;
    }
}

int convert_intvec_from_compatible_numpy_array( PyObject *obj, void *s ) {
    if ( is_intvec_compatible_numpy_array(obj) ) {
	std::vector<int> *to = (std::vector<int>*) s;
	std::vector<int> shape;
	casac::numpy2vector((PyArrayObject*)obj,*to, shape);
	return 1;
    }
    return 0;
}

}
