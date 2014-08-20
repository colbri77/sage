cimport cppWrapper
from sage.structure.sage_object import SageObject

class Generator(SageObject):
    def __init__(self, optimization="enhanced"):
        self.optimization = optimization

    def __cinit__(self):
        pass
    
    def calcBasis(self,generators,modPolynomial):
        matrixFactory = cppWrapper.PyMatrixFactory_uint64(modPolynomial)
        polynomialFactory = cppWrapper.PyPolynomialFactory_uint64()
        monFactory = cppWrapper.PyMonomialFactory_uint64()
        bbt = cppWrapper.PyBorderBasisTools_uint64(matrixFactory,polynomialFactory,monFactory,generators.nvariables(),self.optimization)

        return bbt.calculateBasis(generators)
