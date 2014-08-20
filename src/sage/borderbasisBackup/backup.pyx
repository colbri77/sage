from sage.ext.sage_object import SageObject
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

class Generator(SageObject):
    def __init__(self, optimization="enhanced"):
        self.optimization = optimization
    
    def calcBasis(self,generators,modPolynomial):
        indet = generators.nvariables()
        matrixFactory = <IMatrixFactory[long long unsigned int]*>(new MatrixFactory_Fn(modPolynomial))
        polFactory = new PolynomialFactory[long long unsigned int](POLTYPE_VECTOR)
        monFactory = new MonomialFactory[long long unsigned int](MONOMIALTYPE_DEGLEX)
        
        cdef optimization = NONE
        if(self.optimization == "none"):
            optimization = NONE
        elif(self.optimization == "enhanced"):
            optimization = ENHANCED
        elif(self.optimization == "optimistic"):
            optimization = OPTIMISTIC
        elif(self.optimization == "experimental"):
            optimization = EXPERIMENTAL
        else:
            raise ValueError("optimization value \""+self.optimization+"\" unknown")

        bbt = new BorderBasisTools[long long unsigned int](matrixFactory,polFactory,monFactory,indet,optimization)

        nativeIn = self._toNativePolList(generators)
        nativeOut = (new OwningVector[IPolynomial*]())
        nativeOrderIdeal = polFactory.create(indet)
        
        bbt.calculateBasis(nativeIn,nativeOut,nativeOrderIdeal)

        polynomials = self._fromNativePolList(nativeOut,generators)
        orderIdeal = self._fromNativePol(nativeOrderIdeal,generators)

        del bbt
        del nativeIn
        del nativeOut
        del nativeOrderIdeal

        return (polynomials, orderIdeal)

    cdef _Object _toNativePolList(self,generators) except *:
        cdef IOwningList[IPolynomial*]* result = <IOwningList[IPolynomial*]*>(new OwningVector[IPolynomial*]())
        indet = generators.nvariables()
        monFactory = new MonomialFactory[long long unsigned int](MONOMIALTYPE_DEGLEX)
        polFactory = new PolynomialFactory[long long unsigned int](POLTYPE_VECTOR)
        for generator in generators:
            values = generator.dict()
            pol = polFactory.create(indet)
            for exponents in values:
                monomial = monFactory.create(indet)
                coef = values[exponents]
                for expPos in range(0,indet):
                    exp = exponents[expPos]
                    monomial.set(index=expPos,value=exp)
                term = new Term[long long unsigned int](coef,monomial)
                pol.push(term)
            result.push_back(pol)
        return result

    def _fromNativePolList(self,nativeList,generators):
        polList = []
        ring = generators.ring()

        nativeListSize = nativeList.size()
        for polPos in range(0,nativeListSize):
            pol = nativeList.at(polPos)
            sagePol = self._fromNativePol(pol,generators)    
            polList.append(sagePol)
        result = PolynomialSequence(ring,polList)
        return result

    def _fromNativePol(self,nativePol,generators):
        variables = generators.ring().variable_names_recursive()
        sagePol = 0
        termSize = nativePol.size()
        for termPos in range(0,termSize):
            term = nativePol.at(termPos)
            coef = term.getCoef()
            monomial = term.getMonomial()
            monomialSize = monomial.getIndet()
            sageMon = 1
            for expPos in range(0,monomialSize):
                exp = monomial.at(expPos)
                sageMon = sageMon * variables[expPos]
            sagePol = sagePol + coef * sageMon
        return sagePol




