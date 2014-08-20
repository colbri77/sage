from c_src cimport *

cdef class PyIOwningList_pol:
    cdef IOwningList[IPolynomial_uint64*]* thisptr

cdef class PyIPolynomial_uint64:
    cdef IPolynomial[uint64_t]* thisptr

cdef class PyMonomialFactory_uint64:
    cdef MonomialFactory[uint64_t]* thisptr

cdef class PyPolynomialFactory_uint64:
    cdef PolynomialFactory[uint64_t]* thisptr

cdef class PyMatrixFactory_uint64:
    cdef IMatrixFactory[uint64_t]* thisptr

cdef class PyMatrixFactory_Fn_uint64(PyMatrixFactory_uint64):
    pass

cdef class PyBorderBasisTools_uint64:
    cdef BorderBasisTools[uint64_t]* thisptr
    cdef PyMatrixFactory_uint64 matrixFactory,
    cdef PyPolynomialFactory_uint64 polFactory,
    cdef PyMonomialFactory_uint64 monFactory,
    cdef indet
    cdef optimizations
    cpdef getStatistics(self)
    cpdef calculateBasis(self,generators)
    cdef PyIOwningList_pol _toNativePolList(self,pythonList)
    cdef _fromNativePolList(self,PyIOwningList_pol nativeList,ring)
    cdef _fromNativePol(self,PyIPolynomial_uint64 nativePol,ring)

