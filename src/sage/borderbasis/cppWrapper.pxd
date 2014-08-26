from c_src cimport *

cdef class PyIOwningList_pol:
    cdef IOwningList[IPolynomial_uint64*]* thisptr

cdef class PyIPolynomial_uint64:
    cdef IPolynomial[uint64_t]* thisptr

cdef class PyMonomialFactory:
    cdef MonomialFactory* thisptr

cdef class PyPolynomialFactory_uint64:
    cdef PolynomialFactory[uint64_t]* thisptr

cdef class PyField_uint64:
    cdef IField[uint64_t]* thisptr

cdef class PyFieldFn(PyField_uint64):
    pass

cdef class PyMatrixFactory_uint64:
    cdef IMatrixFactory[uint64_t]* thisptr

cdef class PyMatrixFactory_Fn_uint64(PyMatrixFactory_uint64):
    pass

cdef class PyBorderBasisTools_uint64:
    cdef BorderBasisTools[uint64_t]* thisptr
    cdef PyField_uint64 field
    cdef PyMatrixFactory_uint64 matrixFactory
    cdef PyPolynomialFactory_uint64 polFactory
    cdef PyMonomialFactory monFactory
    cdef indet
    cdef optimizations
    cpdef get_statistics(self)
    cpdef calculate_basis(self,generators)
    cdef PyIOwningList_pol _to_native_pol_list(self,pythonList)
    cdef _from_native_pol_list(self,PyIOwningList_pol nativeList,ring,variables)
    cdef _from_native_pol(self,PyIPolynomial_uint64 nativePol,ring,variables)
    cdef _get_dict(self,polynomial,variables)

