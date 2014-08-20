from c_src cimport *

cdef class PyIOwningList_pol:
    cdef IOwningList[IPolynomial_uint64*]* thisptr

cdef class PyIPolynomial_uint64:
    cdef IPolynomial[long long unsigned int]* thisptr

cdef class PyMonomialFactory_uint64:
    cdef MonomialFactory[long long unsigned int]* thisptr

cdef class PyPolynomialFactory_uint64:
    cdef PolynomialFactory[long long unsigned int]* thisptr

cdef class PyMatrixFactory_uint64:
    cdef IMatrixFactory[long long unsigned int]* thisptr

cdef class PyMatrixFactory_Fn_uint64(PyMatrixFactory_uint64):
    pass

cdef class PyBorderBasisTools_uint64:
    cdef BorderBasisTools[long long unsigned int]* thisptr
    cpdef getStatistics(self)
    cpdef calculateBasis(self,generators)
    cdef PyIOwningList_pol _toNativePolList(self,pythonList)
    cdef _fromNativePolList(self,PyIOwningList_pol nativeList,ring)
    cdef _fromNativePol(self,PyIPolynomial_uint64 nativePol,ring)

