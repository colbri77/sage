r"""
Implementation of the border basis algorithm in 4 different optimization levels

AUTHORS::

- Christian Olbrich (2014): initial version

EXAMPLES::

Small scale example:

    sage: from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
    sage: from sage.borderbasis.generator import BBGenerator

    sage: R.<x,y> = PolynomialRing(GF(2),2)
    sage: F = PolynomialSequence([x*y,y**2+x],R)
    sage: gen = BBGenerator(optimization='none')
    sage: basis,orderIdeal,statistics = gen.calc_basis(generators=F,modPolynomial=2)

    sage: basis
    [y^2 + x, x*y, x^2]
    sage: orderIdeal
    x + y + 1
    sage: statistics
    {'maxMatrix': {'columns': 8L, 'rows': 7L}}

From mini AES:

    sage: from sage.borderbasis.generator import BBGenerator
    sage: sr = mq.SR(2,1,1,4,gf2=True,polybori=False)
    sage: F,s = sr.polynomial_system()
    sage: gen = BBGenerator(optimization='enhanced')
    sage: basis,orderIdeal,statistics = gen.calcBasis(F,2)

    sage: basis
    [k003,
     k002,
     k001 + 1,
     k000,
     s003 + 1,
     s002,
     s001 + 1,
     s000 + 1,
     w103 + 1,
     w102 + 1,
     w101 + 1,
     w100,
     x103,
     x102 + 1,
     x101 + 1,
     x100,
     k103 + 1,
     k102 + 1,
     k101,
     k100,
     s103,
     s102 + 1,
     s101 + 1,
     s100 + 1,
     w203 + 1,
     w202,
     w201,
     w200 + 1,
     x203,
     x202 + 1,
     x201,
     x200,
     k203,
     k202 + 1,
     k201 + 1,
     k200]
"""

#*****************************************************************************
#       Copyright (C) 2014 Christian Olbrich <inforigine@web.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.borderbasis.cppWrapper import *
from sage.structure.sage_object import SageObject

class BBGenerator(SageObject):
    r"""
    The ``BBGenerator`` class is used to calculate a border basis

    EXAMPLES::

        sage: from sage.borderbasis.generator import BBGenerator

        sage: BBGenerator()
        BBGenerator(optimization='enhanced')
        sage: BBGenerator('none')
        BBGenerator(optimization='none')
        sage: BBGenerator(optimization='enhanced')
        BBGenerator(optimization='enhanced')
        sage: BBGenerator(optimization='optimistic')
        BBGenerator(optimization='optimistic')
        sage: BBGenerator(optimization='experimental')
        BBGenerator(optimization='experimental')

    NOTE::

        The default optimization level is 'enhanced'. For the levels 'optimistic' and 'experimental', termination can no longer be proven.
    """
    def __init__(self, optimization="enhanced"):
        r"""
        Generates a ``BBGenerator`` and initializes it with the chosen optimization level

        INPUT::

            - ``optimization`` -- (default 'enhanced') which algorithm variant should be used
        
        EXAMPLES::

            sage: from sage.borderbasis.generator import BBGenerator

            sage: BBGenerator()
            BBGenerator(optimization='enhanced')
        """
        self.optimization = optimization

    def __cinit__(self):
        pass
    
    def calc_basis(self,generators,modPolynomial):
        r"""
        Calculates the border basis of the generator polynomials

        INPUT::

            - ``generators`` -- a ``PolynomialSequence``, containing the generator polynomials
            - ``modPolynomial`` -- the GF-polynomial used in the calculation as integer

        OUTPUT::

            A 3-touple consisting of objects in this order:

            1. The calculated borderbasis as list of polynomials
            2. The according order ideal as polynomial
            3. A map containing statistics and benchmarks, currently only the biggest matrix calculated

        EXAMPLES:

            sage: from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
            sage: from sage.borderbasis.generator import BBGenerator

            sage: R.<x,y> = PolynomialRing(GF(2),2)
            sage: F = PolynomialSequence([x*y,y**2+x],R)
            sage: gen = BBGenerator(optimization='none')
            sage: basis,orderIdeal,statistics = gen.calc_basis(generators=F,modPolynomial=2)

            sage: basis
            [y^2 + x, x*y, x^2]
            sage: orderIdeal
            x + y + 1
            sage: statistics
            {'maxMatrix': {'columns': 8L, 'rows': 7L}}

            
            sage: gen = BBGenerator(optimization='optimistic')
            sage: basis,orderIdeal,statistics = gen.calcBasis(F,2)
            # doesn't seem to terminate: Termination not proven with levels above 'enhanced'

            sage: F = PolynomialSequence([x**2,y**2+1],2)
            sage: basis,orderIdeal,statistics = gen.calcBasis(F,2)
            sage: basis
            [y^2 + 1, x^2, x*y^2 + x, x^2*y]  # 'optimistic' does terminate in this case, though

        WARNING::

            Currently, it is only possible to calculate border bases of polynomials in the galois field.
        """
        matrixFactory = PyMatrixFactory_Fn_uint64(modPolynomial)
        polynomialFactory = PyPolynomialFactory_uint64()
        monFactory = PyMonomialFactory_uint64()
        bbt = PyBorderBasisTools_uint64(matrixFactory,polynomialFactory,monFactory,generators.nvariables(),self.optimization)

        basis,orderIdeal = bbt.calculate_basis(generators)
        statistics = bbt.get_statistics()
        return (basis,orderIdeal,statistics)

    def _latex_(self):
        r"""
        Return latex representation of the BBGenerator. This just displays the class name with the current optimization.

        EXAMPLES::

            sage: from sage.borderbasis.generator import BBGenerator
            sage: gen = BBGenerator(optimization='none')
            sage: latex(gen)
            BBGenerator(optimization='none')
        """
        return "BBGenerator(optimization='"+self.optimization+"')"

    def _repr_(self):
        r"""
        Return the string representation of the BBGenerator. Currently equivalent to _latex_(self).

        EXAMPLES::

            sage: from sage.borderbasis.generator import BBGenerator
            sage: BBGenerator(optimization='experimental')
            BBGenerator(optimization='experimental')
        """
        return self._latex_()

