from c_src cimport *
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

cdef class PyIOwningList_pol:
    def __init__(self):
        pass
    def __cinit__(self):
        pass

cdef class PyIPolynomial_uint64:
    def __init__(self):
        pass
    def __cinit__(self):
        pass

cdef class PyMonomialFactory_uint64:
    def __init__(self):
        pass
    def __cinit__(self):
        self.thisptr = new MonomialFactory[long long unsigned int](MONOMIALTYPE_DEGLEX)
    def __dealloc__(self):
        del self.thisptr

cdef class PyPolynomialFactory_uint64:
    def __init__(self):
        pass
    def __cinit__(self):
        self.thisptr = new PolynomialFactory[long long unsigned int](POLTYPE_VECTOR)
    def __dealloc__(self):
        del self.thisptr

cdef class PyMatrixFactory_uint64:
    def __init__(self):
        pass
    def __cinit__(self):
        pass

cdef class PyMatrixFactory_Fn_uint64(PyMatrixFactory_uint64):
    def __init__(self):
        pass
    def __cinit__(self,minPolynomial):
        self.thisptr = <IMatrixFactory[long long unsigned int]*>(new MatrixFactory_Fn(minPolynomial))
    def __dealloc__(self):
        del self.thisptr

cdef class PyBorderBasisTools_uint64:
    def __init__(self):
        pass

    def __cinit__(self,PyMatrixFactory_uint64 matrixFactory,PyPolynomialFactory_uint64 polFactory,PyMonomialFactory_uint64 monFactory,indeterminates,optimizations):
        self.matrixFactory = matrixFactory
        self.polFactory = polFactory
        self.monFactory = monFactory
        self.indet = indeterminates
        if(optimizations == "none"):
            self.optimizations = NONE
        elif(optimizations == "enhanced"):
            self.optimizations = ENHANCED
        elif(optimizations == "optimistic"):
            self.optimizations = OPTIMISTIC
        elif(optimizations == "experimental"):
            self.optimizations = EXPERIMENTAL
        else:
            raise ValueError("optimization value \""+optimizations+"\" unknown")

        self.thisptr = new BorderBasisTools[long long unsigned int](matrixFactory.thisptr,polFactory.thisptr,monFactory.thisptr,indeterminates,self.optimizations)

    def __dealloc__(self):
        del self.thisptr

    cpdef getStatistics(self):
        cdef Statistics* stats = new Statistics()
        self.thisptr.getStatistics(stats)
        return {'maxMatrix': {'rows': stats.maxMatrix.rows, 'columns': stats.maxMatrix.columns}}

    cpdef calculateBasis(self,generators):
        nativeInWrapper = self._toNativeList(generators)
        cdef IOwningList[IPolynomial_uint64*]* nativeIn = (<PyIOwningList_pol>(nativeInWrapper)).thisptr
        cdef IOwningList[IPolynomial_uint64*]* nativeOut = <IOwningList[IPolynomial_uint64*]*>(new OwningVector[IPolynomial_uint64*]())
        cdef IPolynomial[long long unsigned int]* nativeOrderIdeal = (<PyPolynomialFactory_uint64>(self.polFactory)).thisptr.create(self.indet)

        self.thisptr.calculateBasis(<void*>(nativeIn),<void*>(nativeOut),<void*>(nativeOrderIdeal))  

        nativeOutWrapper = PyIOwningList_pol()
        nativeOrderIdealWrapper = PyIPolynomial_uint64()
        nativeOutWrapper.thisptr = nativeOut
        nativeOrderIdealWrapper.thisptr = nativeOrderIdeal

        polynomials = self._fromNativePolList(nativeOutWrapper,generators.ring())
        orderIdeal = self._fromNativePol(nativeOrderIdealWrapper,generators.ring())

        del nativeIn
        del nativeOut
        del nativeOrderIdeal

        return (polynomials, orderIdeal)
 
    cdef PyIOwningList_pol _toNativePolList(self,pythonList):
        result = <IOwningList[IPolynomial_uint64*]*>(new OwningVector[IPolynomial_uint64*]())
        for generator in pythonList:
            values = pythonList.dict()
            pol = (<PyPolynomialFactory_uint64>(self.polFactory)).thisptr.create(self.indet)
            for exponents in values:
                monomial = (<PyMonomialFactory_uint64>(self.monFactory)).thisptr.create(self.indet)
                coef = values[exponents]
                for expPos in range(0,self.indet):
                    exp = exponents[expPos]
                    monomial.set(expPos,exp)
                term = new Term[long long unsigned int](coef,monomial)
                pol.push(term)
            result.push_back(<IPolynomial_uint64*>(pol))
        wrappedResult = PyIOwningList_pol()
        wrappedResult.thisptr = result
        return wrappedResult

    cdef _fromNativePolList(self,PyIOwningList_pol nativeList,ring):
        polList = []
        nativeListSize = nativeList.thisptr.size()
        for polPos in range(0,nativeListSize):
            pol = nativeList.thisptr.at(polPos)
            polWrapper = PyIPolynomial_uint64()
            polWrapper.thisptr = <IPolynomial[long long unsigned int]*>(pol)
            sagePol = self._fromNativePol(polWrapper,ring)
            polList.append(sagePol)
        result = PolynomialSequence(ring,polList)
        return result

    cdef _fromNativePol(self,PyIPolynomial_uint64 nativePol,ring):
        sagePol = 0
        variables = ring.variable_names_recursive()
        termSize = nativePol.thisptr.size()
        for termPos in range(0,termSize):
            term = nativePol.thisptr.at(termPos)
            coef = term.getCoef()
            monomial = term.getMonomial()
            monomialSize = monomial.getIndet()
            sageMon = 1
            for expPos in range(0,monomialSize):
                exp = monomial.at(expPos)
                sageMon = sageMon * variables[expPos]
            sagePol = sagePol + coef * sageMon
        return sagePol

