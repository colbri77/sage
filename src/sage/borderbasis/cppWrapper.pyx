r"""
Wrapping classes for the border basis algorithm in C++, located under src/

Internal use only!

AUTHORS:

- Christian Olbrich (2014): initial version
"""
from c_src cimport *
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.calculus.var import var
from sage.symbolic.expression_conversions import polynomial

cdef class PyIOwningList_pol:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self):
        pass

cdef class PyIPolynomial_uint64:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self):
        pass

cdef class PyMonomialFactory_uint64:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self):
        self.thisptr = new MonomialFactory[uint64_t](MONOMIALTYPE_DEGLEX)
    def __dealloc__(self):
        del self.thisptr

cdef class PyPolynomialFactory_uint64:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self):
        self.thisptr = new PolynomialFactory[uint64_t](POLTYPE_VECTOR)
    def __dealloc__(self):
        del self.thisptr

cdef class PyMatrixFactory_uint64:
    r"""
    Parameter class, used to send C++ classes to python methods
    """    
    def __cinit__(self):
        pass

cdef class PyMatrixFactory_Fn_uint64(PyMatrixFactory_uint64):
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self,minPolynomial):
        self.thisptr = <IMatrixFactory[uint64_t]*>(new MatrixFactory_Fn(minPolynomial))
    def __dealloc__(self):
        del self.thisptr

cdef class PyBorderBasisTools_uint64:
    r"""
    Wrapping class for the C++ class ``BorderBasisTools``

    INPUT::

        - ``matrixFactory`` -- the matrix factory to use during calculation
        - ``polFactory`` -- the polynomial factory to use for generating polynomials
        - ``monFactory`` -- the monomial factory to use for generating monomials
        - ``optimizations`` -- a string describing which algorithm to use
    
    EXAMPLES::

        sage: from sage.borderbasis.cppWrapper import *
        sage: from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

        sage: R.<x,y> = PolynomialRing(GF(2),2)
        sage: F = PolynomialSequence([x*y,y**2+x],R)

        sage: matrixFactory = PyMatrixFactory_Fn_uint64(2)
        sage: polynomialFactory = PyPolynomialFactory_uint64()
        sage: monFactory = PyMonomialFactory_uint64()

        sage: PyBorderBasisTools_uint64(matrixFactiry,polynomialFactory,monFactory,F.nvariables(),'enhanced')
        <sage.borderbasis.cppWrapper.PyBorderBasisTools_uint64>

    """
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

    cpdef get_statistics(self):
        r"""
        Returns the statistics collected during the last run of ``calculate_basis``

        OUTPUT::

            A map with the following structure:
            {'maxMatrix': {'columns': <int>, 'rows': <int>}}
            ``maxMatrix`` describes the biggest matrix the algorithm had to handle during calculation
        
        EXAMPLES::
            sage: from sage.borderbasis.cppWrapper import *
            sage: sr = mq.SR(2,1,1,4,gf2=True,polybori=False)
            sage: F,s = sr.polynomial_system()

            sage: matrixFactory = PyMatrixFactory_Fn_uint64(2)
            sage: polynomialFactory = PyPolynomialFactory_uint64()
            sage: monFactory = PyMonomialFactory_uint64()

            sage: bbt = PyBorderBasisTools_uint64(matrixFactiry,polynomialFactory,monFactory,F.nvariables(),'enhanced')
            sage: Basis,orderIdeal = bbt.calculate_basis(F)

            sage: bbt.get_statistics()
            {'maxMatrix': {'columns': 6576L, 'rows': 21203L}}
        """
        cdef Statistics* stats = new Statistics()
        self.thisptr.getStatistics(stats)
        return {'maxMatrix': {'rows': stats.maxMatrix.rows, 'columns': stats.maxMatrix.columns}}

    cpdef calculate_basis(self,generators):
        r"""
        Calculates a border basis from the given generator polynomials

        INPUT::

            - ``generators`` -- A PolynomialSequence of generator polynomials

        OUTPUT::

            A 2-touple consisting of objects in this order:

            1. A list of polynomials, the calculated border basis
            2. A polynomial, the according order ideal
        """
        nativeInWrapper = self._to_native_pol_list(generators)
        cdef IOwningList[IPolynomial_uint64*]* nativeIn = (<PyIOwningList_pol>(nativeInWrapper)).thisptr
        cdef IOwningList[IPolynomial_uint64*]* nativeOut = <IOwningList[IPolynomial_uint64*]*>(new OwningVector[IPolynomial_uint64*]())
        cdef IPolynomial[uint64_t]* nativeOrderIdeal = (<PyPolynomialFactory_uint64>(self.polFactory)).thisptr.create(self.indet)

        self.thisptr.calculateBasis(<void*>(nativeIn),<void*>(nativeOut),<void*>(nativeOrderIdeal))  

        nativeOutWrapper = PyIOwningList_pol()
        nativeOrderIdealWrapper = PyIPolynomial_uint64()
        nativeOutWrapper.thisptr = nativeOut
        nativeOrderIdealWrapper.thisptr = nativeOrderIdeal

        polynomials = self._from_native_pol_list(nativeOutWrapper,generators.ring())
        orderIdeal = self._from_native_pol(nativeOrderIdealWrapper,generators.ring())

        del nativeIn
        del nativeOut
        del nativeOrderIdeal

        return (polynomials, orderIdeal)
 
    cdef PyIOwningList_pol _to_native_pol_list(self,pythonList):
        r"""
        Converts a python list of sage polynomials to a PyIOwningList_pol, which contains a pointer used for native calls

        No example, since intermixed with C++ objects and not applicable to pure python sage

        INPUT::

            - ``pythonList`` -- A PolynomialSequence containing sage polynomials to parse

        OUTPUT::

            A PyIOwningList_pol, equivalent to the provided ``pythonList``
        """
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

    cdef _from_native_pol_list(self,PyIOwningList_pol nativeList,ring):
        r"""
        Converts a C++ list of C++ polynomials to a python list of sage polynomials

        No example, since handling C++ objects (not applicable to sage) 

        INPUT::

            - ``nativeList`` -- a PyIOwningList_pol, containing a C++ list of C++ polynomials to parse
            - ``ring`` -- the polynomial ring associated with the polynomials

        OUTPUT::

            A python list of sage polynomials equivalent to ``nativeList``
        """
        polList = []
        nativeListSize = nativeList.thisptr.size()
        for polPos in range(0,nativeListSize):
            pol = nativeList.thisptr.at(polPos)
            polWrapper = PyIPolynomial_uint64()
            polWrapper.thisptr = <IPolynomial[uint64_t]*>(pol)
            sagePol = self._from_native_pol(polWrapper,ring)
            polList.append(sagePol)
        #result = PolynomialSequence(polList,ring)
        # There are like 100 different types of polynomials in sage, and all seem to be incompatible to each other -.-
        # currently I have no idea how to get from this polynomial to one that the PolynomialSequence accepts...
        result = polList
        return result

    cdef _from_native_pol(self,PyIPolynomial_uint64 nativePol,ring):
        r"""
        Converts a C++ polynomial to a sage polynomial

        No examples, since handling function of C++ objects (not applicable for sage)

        INPUT::

            - ``nativePol`` -- a PyIPolynomial_uint64, containing the C++ polynomial to parse
            - ``ring`` -- the polynomial ring associated with the polynomial

        OUTPUT::

            A sage polynomial equivalent to the C++ polynomial
        """
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

