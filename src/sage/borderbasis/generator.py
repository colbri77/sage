from sage.borderbasis.cppWrapper import *
from sage.structure.sage_object import SageObject

class BBGenerator(SageObject):
    def __init__(self, optimization="enhanced"):
        self.optimization = optimization

    def __cinit__(self):
        pass
    
    def calcBasis(self,generators,modPolynomial):
        matrixFactory = PyMatrixFactory_Fn_uint64(modPolynomial)
        polynomialFactory = PyPolynomialFactory_uint64()
        monFactory = PyMonomialFactory_uint64()
        bbt = PyBorderBasisTools_uint64(matrixFactory,polynomialFactory,monFactory,generators.nvariables(),self.optimization)

        return bbt.calculateBasis(generators)
