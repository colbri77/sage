r"""
Wrapping classes for the border basis algorithm in C++, located under src/

Internal use only!

AUTHORS:

- Christian Olbrich (2014): initial version
"""
cimport cython
from libc.stdlib cimport malloc, free
from c_src cimport *
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
from sage.calculus.var import var
from sage.symbolic.expression_conversions import polynomial
from sage.rings.polynomial.pbori import BooleanPolynomial

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

cdef class PyMonomialFactory:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self,use_positions,indet,gf2,order):
        if(gf2):
            if use_positions and order=="deglex":
                self.thisptr = <IMonomialFactory*>(new MonomialFactoryDegLexGF2(indet))
            elif (not use_positions) and order=="deglex":
                self.thisptr = <IMonomialFactory*>(new MonomialFactoryNoOrderPosGF2(indet))
            elif use_positions and order=="degrevlex":
                self.thisptr = <IMonomialFactory*>(new MonomialFactoryDegRevLexGF2(indet))
        else:
            if use_positions and order=="deglex":
                self.thisptr = <IMonomialFactory*>(new MonomialFactoryDegLex(indet))
            elif (not use_positions) and order=="deglex":
                self.thisptr = <IMonomialFactory*>(new MonomialFactoryNoOrderPos(indet))
            elif use_positions and order=="degrevlex":
                self.thisptr = <IMonomialFactory*>(new MonomialFactoryDegRevLex(indet))
    def __dealloc__(self):
        del self.thisptr

cdef class PyPolynomialFactory_uint64:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self,gf2):
        if(gf2):
            self.thisptr = new PolynomialFactory[uint64_t](POLTYPE_VECTOR_GF2)
        else:
            self.thisptr = new PolynomialFactory[uint64_t](POLTYPE_VECTOR)
    def __dealloc__(self):
        del self.thisptr

cdef class PyField_uint64:
    r"""
    Parameter class, used to send C++ classes to python methods
    """    
    def __cinit__(self):
        pass

cdef class PyFieldFn(PyField_uint64):
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self,isNull,minPolynomial):
        if isNull:
            self.thisptr = <IField[uint64_t]*>NULL
        else:
            self.thisptr = <IField[uint64_t]*>(new FieldFn(minPolynomial))
    def __dealloc__(self):
        if self.thisptr != <IField[uint64_t]*>NULL:
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
    def __cinit__(self,isNull,minPolynomial): 
        if isNull:
            self.thisptr = <IMatrixFactory[uint64_t]*>NULL
        else:
            self.thisptr = <IMatrixFactory[uint64_t]*>(new MatrixFactory_Fn(minPolynomial))
    def __dealloc__(self):
        if self.thisptr != <IMatrixFactory[uint64_t]*>NULL:
            del self.thisptr

cdef class PyBBConfig:
    r"""
    Parameter class, used to send C++ classes to python methods
    """
    def __cinit__(self,indet,OptLevel opt,PyField_uint64 field,PyMatrixFactory_uint64 mFac,PyPolynomialFactory_uint64 polFac,PyMonomialFactory monFac,use_pol_ex,use_variable_exclusion,variable_exclusions,use_gf2_reductions,min_mutants_limit):
        self.exclusions = <bool*>malloc(len(variable_exclusions)*cython.sizeof(bool))
        for i in xrange(len(variable_exclusions)):
            self.exclusions[i] = variable_exclusions[i]
        self.thisptr = <BBConfig*>(new BBConfig(indet,opt,field.thisptr,mFac.thisptr,polFac.thisptr,monFac.thisptr,use_pol_ex,use_variable_exclusion,self.exclusions,use_gf2_reductions,min_mutants_limit))
    def __dealloc__(self):
        del self.thisptr
        free(self.exclusions)

cdef class PyBorderBasisTools_uint64:
    r"""
    Wrapping class for the C++ class ``BorderBasisTools``

    INPUT::

        - ``field`` -- the field to use during calculation (if matrix is used with a null ptr)
        - ``matrixFactory`` -- the matrix factory to use during calculation (if field is used with a null ptr)
        - ``polFactory`` -- the polynomial factory to use for generating polynomials
        - ``monFactory`` -- the monomial factory to use for generating monomials
        - ``optimizations`` -- a string describing which algorithm to use
        - ``use_pol_exclusion`` -- whether to exclude multiples of minimal, detected polynomials
        - ``use_variable_exclusion`` -- whether to exclude variables from the calculation on the fly
        - ``variable_exclusions`` -- if variable exclusion is used, which variables to keep
        - ``use_gf2_reductions`` -- whether the algorithm should reduce polynomials by the field pols (x²+x)
        - ``min_mutants_limit`` -- the minimal amount of mutants to be used for extension in mutant algorithms (in fraction of columns)   
 
    EXAMPLES::

        sage: from sage.borderbasis.cppWrapper import *
        sage: from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

        sage: R.<x,y> = PolynomialRing(GF(2),2)
        sage: F = PolynomialSequence([x*y,y**2+x],R)

        sage: field = PyFieldFn(false,2)
        sage: matrixFactory = PyMatrixFactory_Fn_uint64(True,2)
        sage: polynomialFactory = PyPolynomialFactory_uint64(False)
        sage: monFactory = PyMonomialFactory(True,2,False,"deglex")

        sage: PyBorderBasisTools_uint64(field,polynomialFactory,monFactory,F.nvariables(),'enhanced')
        <sage.borderbasis.cppWrapper.PyBorderBasisTools_uint64>

    """
    def __cinit__(self,PyField_uint64 field,PyMatrixFactory_uint64 matrixFactory,PyPolynomialFactory_uint64 polFactory,PyMonomialFactory monFactory,indeterminates,optimizations,use_pol_exclusion=False,use_variable_exclusion=False,variable_exclusions=[False,],use_gf2_reductions=True,min_mutants_limit=0):
        self.field = field
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
        elif(optimizations == "mutant"):
            self.optimizations = MUTANT
        elif(optimizations == "improved_mutant"):
            self.optimizations = IMPROVED_MUTANT
        elif(optimizations == "improved_mutant_linear"):
            self.optimizations = IMPROVED_MUTANT_LINEAR
        elif(optimizations == "improved_mutant_optimistic"):
            self.optimizations = IMPROVED_MUTANT_OPTIMISTIC
        else:
            raise ValueError("optimization value \""+optimizations+"\" unknown")

        self.cfg = PyBBConfig(indeterminates,self.optimizations,field,matrixFactory,polFactory,monFactory,use_pol_exclusion,use_variable_exclusion,variable_exclusions,use_gf2_reductions,min_mutants_limit)
        
        self.thisptr = new BorderBasisTools[uint64_t](self.cfg.thisptr)

    def __dealloc__(self):
        del self.thisptr

    cpdef get_statistics(self):
        r"""
        Returns the statistics collected during the last run of ``calculate_basis``

        OUTPUT::

            A map with the following structure:
            {'maxComparisons': <int>,'maxMatrix': {'rows': <int>, 'columns': <int>}}
            ``maxComparisons`` describes the biggest amount of comparisons between terms the algorithm had to handle during one reduction step - only filled if not matrices are used
            ``maxMatrix`` describes the biggest matrix the algorithm had to calculate during one reduction step - only filled if no field is used
        
        EXAMPLES::
            sage: from sage.borderbasis.cppWrapper import *
            sage: sr = mq.SR(2,1,1,4,gf2=True,polybori=False)
            sage: F,s = sr.polynomial_system()

            sage: field = PyFieldFn(True,2)
            sage: matrixFactory = PyMatrixFactory_Fn_uint64(False,2)
            sage: polynomialFactory = PyPolynomialFactory_uint64(False)
            sage: monFactory = PyMonomialFactory(True,F.nvariables(),False,"degrevlex")

            sage: bbt = PyBorderBasisTools_uint64(field,matrixFactory,polynomialFactory,monFactory,F.nvariables(),'enhanced')
            sage: Basis,orderIdeal = bbt.calculate_basis(F)

            sage: bbt.get_statistics()
            {'maxComparisons': 17623L}
        """
        cdef Statistics* stats = new Statistics()
        self.thisptr.getStatistics(stats)
        return {'maxComparisons': stats.max_comparisons_in_reduction,'maxMatrix': {"rows": stats.maxMatrix.rows, "columns": stats.maxMatrix.columns}}

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

        polynomials = self._from_native_pol_list(nativeOutWrapper,generators.ring(),self._get_variables(generators))
        orderIdeal = self._from_native_pol(nativeOrderIdealWrapper,generators.ring(),self._get_variables(generators))

        del nativeIn
        del nativeOut
        del nativeOrderIdeal

        return (polynomials, orderIdeal)

    cdef _get_variables(self,pythonList):
        existing = pythonList.variables()
        ordered = pythonList.ring().gens()
        result = []
        for i in ordered:
            for e in existing:
                if ("%s" % (i)) == ("%s" % (e)):
                    result.append(i)
                    break
        return result
 
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
            values = self._get_dict(generator,self._get_variables(pythonList))
            pol = (<PyPolynomialFactory_uint64>(self.polFactory)).thisptr.create(self.indet)
            for exponents in values:
                monomial = (<PyMonomialFactory>(self.monFactory)).thisptr.create()
                coef = values[exponents]
                try:
                    coef = (int)(coef.int_repr())
                except:
                    pass
                for expPos in range(0,self.indet):
                    exp = exponents[expPos]
                    if(exp>0):
                        monomial = <IMonomial*>(monomial.set(expPos,exp))
                term = new Term[uint64_t](coef,monomial)
                pol.push(term)
            result.push_back(<IPolynomial_uint64*>(pol))
        wrappedResult = PyIOwningList_pol()
        wrappedResult.thisptr = result
        return wrappedResult

    cdef _from_native_pol_list(self,PyIOwningList_pol nativeList,ring,variables):
        r"""
        Converts a C++ list of C++ polynomials to a python list of sage polynomials

        No example, since handling C++ objects (not applicable to sage) 

        INPUT::

            - ``nativeList`` -- a PyIOwningList_pol, containing a C++ list of C++ polynomials to parse
            - ``ring`` -- the polynomial ring associated with the polynomials
            - ``variables`` -- a list of variables to assign use for the conversion (ordered)

        OUTPUT::

            A python list of sage polynomials equivalent to ``nativeList``
        """
        polList = []
        nativeListSize = nativeList.thisptr.size()
        for polPos in range(0,nativeListSize):
            pol = nativeList.thisptr.at(polPos)
            polWrapper = PyIPolynomial_uint64()
            polWrapper.thisptr = <IPolynomial[uint64_t]*>(pol)
            sagePol = self._from_native_pol(polWrapper,ring,variables)
            polList.append(sagePol)
        result = polList
        return result

    cdef _from_native_pol(self,PyIPolynomial_uint64 nativePol,ring,variables):
        r"""
        Converts a C++ polynomial to a sage polynomial

        No examples, since handling function of C++ objects (not applicable for sage)

        INPUT::

            - ``nativePol`` -- a PyIPolynomial_uint64, containing the C++ polynomial to parse
            - ``ring`` -- the polynomial ring associated with the polynomial
            - ``variables`` -- a list of variable names to assign to the variable positions

        OUTPUT::

            A sage polynomial equivalent to the C++ polynomial
        """
        sagePol = var(variables[0])-var(variables[0])
        termSize = nativePol.thisptr.size()
        for termPos in range(0,termSize):
            term = nativePol.thisptr.at(termPos)
            coef = term.getCoef()
            try:
                coef = ring.base_ring().fetch_int(coef)
            except:
                pass
            monomial = term.getMonomial()
            monomialSize = monomial.getIndet()
            sageMon = 1
            for expPos in range(0,monomialSize):
                exp = monomial.at(expPos)
                sageMon = sageMon * (var(variables[expPos])**exp)
            sagePol = sagePol + coef * sageMon
        sagePol = polynomial(sagePol,ring)
        return sagePol

    cdef _get_dict(self,polynomial,variables):
        r"""
        Builds a dictionary from the given polynomial.

        INPUT::

            - ``polynomial`` -- a polynomial of undefined class
            - ``variables`` -- a list of variables used for the conversion

        OUTPUT::

            A dictionary describing the polynomial
        """
        if(type(polynomial)==BooleanPolynomial):
            result = {}
            mapping = {}
            curPos = 0
            for var in variables:
                mapping[var] = curPos
                curPos = curPos + 1
            for term in polynomial.terms():
                key = [0]*(len(variables))
                for var in term.variables():
                    key[mapping[var]] = 1 #only option in BooleanPolynomial
                key = tuple(key)
                result[key] = 1
            return result
        else:
            result = {}
            mapping = {}
            curPos = 0
            for var in variables:
                mapping[var] = curPos
                curPos = curPos + 1
            termCtr = 0
            coefficients = polynomial.coefficients()
            for term in polynomial.monomials():
                key = [0]*(len(variables))
                exps = term.exponents()[0]
                pos = -1
                for var in term.variables():
                    pos = pos + 1
                    while exps[pos]==0 and pos<len(exps):
                        pos = pos + 1
                    key[mapping[var]] = exps[pos]
                key = tuple(key)
                coef = 1
                try:
                    coef = (int)(coefficients[termCtr].int_repr())
                except:
                    try:
                        coef = (int)(coefficients[termCtr])
                    except:
                        pass
                result[key] = coef
                termCtr = termCtr + 1
            return result



