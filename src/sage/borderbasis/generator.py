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
    sage: gen = BBGenerator(optimization='none',use_positions=True)
    sage: basis,orderIdeal,statistics = gen.calc_basis(generators=F,modPolynomial=2)

    sage: basis
    [y^2 + x, x*y, x^2]
    sage: orderIdeal
    x + y + 1
    sage: statistics
    {'maxComparisons': }

From mini AES:

    sage: from sage.borderbasis.generator import BBGenerator
    sage: sr = mq.SR(2,1,1,4,gf2=True,polybori=False)
    sage: F,s = sr.polynomial_system()
    sage: gen = BBGenerator(optimization='enhanced')
    sage: basis,orderIdeal,statistics = gen.calcBasis(F,2)

    sage: basis
    [k003,
     k002 + 1,
     k000,
     s003 + 1,
     k001 + s002,
     k001 + s001,
     k001 + s000 + 1,
     w103 + 1,
     w102 + 1,
     k001 + w101 + 1,
     w100,
     x103,
     x102 + 1,
     x101 + 1,
     k001 + x100,
     k103,
     k001 + k102,
     k101 + 1,
     k100,
     s103 + 1,
     k001 + s102,
     s101 + 1,
     k001 + s100 + 1,
     w203,
     w202 + 1,
     k001 + w201 + 1,
     k001 + w200 + 1,
     x203 + 1,
     k001 + x202 + 1,
     x201,
     k001 + x200,
     k001 + k203,
     k202,
     k001 + k201,
     k200,
     k001*k003,
     k001*k002 + k001,
     k001^2 + k001,
     k000*k001,
     k001*s003 + k001,
     k001*s002 + k001,
     k001*s001 + k001,
     k001*s000,
     k001*w103 + k001,
     k001*w102 + k001,
     k001*w101,
     k001*w100,
     k001*x103,
     k001*x102 + k001,
     k001*x101 + k001,
     k001*x100 + k001,
     k001*k103,
     k001*k102 + k001,
     k001*k101 + k001,
     k001*k100,
     k001*s103 + k001,
     k001*s102 + k001,
     k001*s101 + k001,
     k001*s100,
     k001*w203,
     k001*w202 + k001,
     k001*w201,
     k001*w200,
     k001*x203 + k001,
     k001*x202,
     k001*x201,
     k001*x200 + k001,
     k001*k203 + k001,
     k001*k202,
     k001*k201 + k001,
     k001*k200]
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
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence as PS
from sage.rings.polynomial.pbori import BooleanPolynomial
from sage.rings.polynomial.pbori import BooleanMonomial
from sage.rings.polynomial.pbori import BooleanMonomialMonoid
from sage.rings.polynomial.pbori import BooleanPolynomialRing

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
        BBGenerator(optimization='enhanced',use_positions=False)
        sage: BBGenerator(optimization='optimistic')
        BBGenerator(optimization='optimistic',use_positions=True)
        sage: BBGenerator(optimization='experimental')
        BBGenerator(optimization='experimental')

    NOTE::

        The default optimization level is 'enhanced'. For the levels 'optimistic' and 'experimental', termination can no longer be proven.
    """
    def __init__(self, optimization="optimistic", use_positions=True, use_matrix=True, use_autoreduction=True, use_pol_exclusion=False, order="degrevlex"):
        r"""
        Generates a ``BBGenerator`` and initializes it with the chosen optimization level

        INPUT::

            - ``optimization`` -- (default 'enhanced') which algorithm variant should be used
            - ``use_positions`` -- (default True) whether to use calculated DegLex-Positions in the algorithm        

        EXAMPLES::

            sage: from sage.borderbasis.generator import BBGenerator

            sage: BBGenerator()
            BBGenerator(optimization='enhanced')
        """
        self.order = order
        self.optimization = optimization
        self.use_positions = use_positions
        self.use_matrix = use_matrix
        self.use_autoreduction = use_autoreduction
        self.use_pol_exclusion = use_pol_exclusion

        if(use_matrix and not use_positions):
            raise RuntimeError("use_matrix needs use_positions enabled")
        if(order!="deglex" and order!="degrevlex"):
            raise RuntimeError("unknown order: use deglex or degrevlex")
        if((not use_positions) and order!="deglex"):
            raise RuntimeError("positions are necessary for all orders but deglex")

    def __cinit__(self):
        pass
    
    def calc_basis(self,generators,modPolynomial,reduce_monomials_but=None,keep_only=None):
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
        use_variable_exclusion = False
        if(keep_only != None):
            use_variable_exclusion = True

            variables = generators.variables()
            bool_keep_map = [False]*len(variables)
            for var in keep_only:
                for i in range(0,len(variables)):
                    if ("%s" % var) == ("%s" % (variables[i])):
                        bool_keep_map[i] = True
            keep_only = bool_keep_map
        else:
            keep_only = [False,]

        if(use_variable_exclusion and (not self.use_autoreduction or not modPolynomial==2)):
            raise RuntimeError("Variable exclusion is currently only enabled in GF(2) with enabled autoreduction.")

        if(reduce_monomials_but != None):
            generators = self.shrink_system(generators,reduce_monomials_but,self.use_autoreduction and modPolynomial==2)
        field = None
        matrix = None
        if(not self.use_matrix):
            field = PyFieldFn(False,modPolynomial)
            matrix = PyMatrixFactory_Fn_uint64(True,modPolynomial)
        else:
            field = PyFieldFn(True,modPolynomial)
            matrix = PyMatrixFactory_Fn_uint64(False,modPolynomial)
        polynomialFactory = PyPolynomialFactory_uint64(modPolynomial==2 and self.use_autoreduction)
        monFactory = PyMonomialFactory(self.use_positions,generators.nvariables(),modPolynomial==2 and self.use_autoreduction,self.order)
        bbt = PyBorderBasisTools_uint64(field,matrix,polynomialFactory,monFactory,generators.nvariables(),self.optimization,self.use_pol_exclusion,use_variable_exclusion,keep_only)

        basis,orderIdeal = bbt.calculate_basis(generators)
        statistics = bbt.get_statistics()
        
        del field
        del polynomialFactory
        del monFactory
        del bbt

        return (basis,orderIdeal,statistics)

    def shrink_system(self,plist,keyvars,polybori=False,max_degree=0x7fffffff):
        # if we should autoreduce the system, we port it to BooleanPolynomial, where this stuff gets done automatically
        if polybori and type(plist[0])!=BooleanPolynomial:
            newRing = BooleanPolynomialRing(len(plist.ring().variable_names()),plist.ring().variable_names())
            plistNew = PS([],newRing)
            for p in plist:
                pNew = BooleanPolynomial(newRing)
                for m in p.monomials():
                    mNew = BooleanPolynomial(newRing) + 1
                    for v in m.variables():
                        mNew = mNew * v
                    pNew = pNew + m
                plistNew.append(pNew)
            plist = plistNew
        # now we shrink the system
        variables = plist.ring().gens()
        for v in variables:
            useful = True
            for kv in keyvars:
                if kv == v:
                    useful = False
            if not useful:
                continue

            index = self._find_reduceable_index(plist,v)
            if index!=-1:
                 B = plist[index]+v
                 plist = self._reduce(plist,v,B,max_degree)
        return plist

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

    def _reduce(self,plist,A,B,max_degree):
        result = PS([],plist.ring())
        for pol in plist:
            pNew = plist[0]-plist[0]
            # subs() is buggy in Polybori - otherwise, just use that
            if type(pol)!=BooleanPolynomial:
                pNew = pol.subs({A.monomials()[0]: B})
                if pNew != plist[0]-plist[0] and pNew.degree()<=max_degree:
                    result.append(pNew)
                elif pNew != plist[0]-plist[0]:
                    result.append(pol)
                continue
            for monomial in pol.monomials():
                reducable = False
                reducable = monomial.reducible_by(A.monomials()[0])
                if reducable:
                    monNew = 1
                    try:
                        monNew = monomial/A.monomials()[0]
                    except:
                        for v in monomial.variables():
                            if v!=A.monomials()[0].variables()[0]:
                                monNew = monNew * v
                    pNew = pNew + monNew*B
                else:
                    pNew = pNew + monomial
            if pNew != plist[0]-plist[0] and pNew.degree()<=max_degree:
                result.append(pNew)
            elif pNew != plist[0]-plist[0]:
                result.append(pol)
        return result

    def _find_reduceable_index(self,plist,var):
        i = -1
        for pol in plist:
            i = i + 1
            pol = pol-var
            useable = True
            for mon in pol.monomials():
                reducible = False
                try:
                    reducible = mon.reducible_by(var.monomials()[0])
                except:
                    reducible = var.monomials()[0].divides(mon)
                if reducible:
                    useable = False
                    break
            if useable and pol.degree()<3:
                return i
        return -1


