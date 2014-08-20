from c_src cimport *
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.calculus.var import var
from sage.symbolic.expression_conversions import polynomial

cdef class PyIOwningList_pol:
    def __cinit__(self):
        pass

cdef class PyIPolynomial_uint64:
    def __cinit__(self):
        pass

cdef class PyMonomialFactory_uint64:
    def __cinit__(self):
        self.thisptr = new MonomialFactory[uint64_t](MONOMIALTYPE_DEGLEX)
    def __dealloc__(self):
        del self.thisptr

cdef class PyPolynomialFactory_uint64:
    def __cinit__(self):
        self.thisptr = new PolynomialFactory[uint64_t](POLTYPE_VECTOR)
    def __dealloc__(self):
        del self.thisptr

cdef class PyMatrixFactory_uint64:
    def __cinit__(self):
        pass

cdef class PyMatrixFactory_Fn_uint64(PyMatrixFactory_uint64):
    def __cinit__(self,minPolynomial):
        self.thisptr = <IMatrixFactory[uint64_t]*>(new MatrixFactory_Fn(minPolynomial))
    def __dealloc__(self):
        del self.thisptr

cdef class PyBorderBasisTools_uint64:
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

        self.thisptr = new BorderBasisTools[uint64_t](matrixFactory.thisptr,polFactory.thisptr,monFactory.thisptr,indeterminates,self.optimizations)

    def __dealloc__(self):
        del self.thisptr

    cpdef getStatistics(self):
        cdef Statistics* stats = new Statistics()
        self.thisptr.getStatistics(stats)
        return {'maxMatrix': {'rows': stats.maxMatrix.rows, 'columns': stats.maxMatrix.columns}}

    cpdef calculateBasis(self,generators):
        nativeInWrapper = self._toNativePolList(generators)
        cdef IOwningList[IPolynomial_uint64*]* nativeIn = (<PyIOwningList_pol>(nativeInWrapper)).thisptr
        cdef IOwningList[IPolynomial_uint64*]* nativeOut = <IOwningList[IPolynomial_uint64*]*>(new OwningVector[IPolynomial_uint64*]())
        cdef IPolynomial[uint64_t]* nativeOrderIdeal = (<PyPolynomialFactory_uint64>(self.polFactory)).thisptr.create(self.indet)

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
            values = generator.dict()
            pol = (<PyPolynomialFactory_uint64>(self.polFactory)).thisptr.create(self.indet)
            for exponents in values:
                monomial = (<PyMonomialFactory_uint64>(self.monFactory)).thisptr.create(self.indet)
                coef = values[exponents]
                for expPos in range(0,self.indet):
                    exp = exponents[expPos]
                    monomial.set(expPos,exp)
                term = new Term[uint64_t](coef,monomial)
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
            polWrapper.thisptr = <IPolynomial[uint64_t]*>(pol)
            sagePol = self._fromNativePol(polWrapper,ring)
            polList.append(sagePol)
        #result = PolynomialSequence(polList,ring)
        # There are like 100 different types of polynomials in sage, and all seem to be incompatible to each other -.-
        # currently I have no idea how to get from this polynomial to one that the PolynomialSequence accepts...
        result = polList
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
                sageMon = sageMon * (var(variables[expPos])**exp)
            sagePol = sagePol + coef * sageMon
        sagePol = polynomial(sagePol,ring)
        return sagePol

