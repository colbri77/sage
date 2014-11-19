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
    sage: gen = BBGenerator(optimization='none',use_positions=True,use_autoreduction=False)
    sage: basis,orderIdeal,statistics = gen.calc_basis(generators=F,modPolynomial=2)

    sage: basis
    [y^2 + x, x*y, x^2]
    sage: orderIdeal
    x + y + 1
    sage: statistics
    {'maxComparisons': 0L, 'maxMatrix': {'columns': 7L, 'rows': 6L}}

From mini AES:

    sage: from sage.borderbasis.generator import BBGenerator
    sage: sr = mq.SR(2,1,1,4,gf2=True,polybori=False)
    sage: F,s = sr.polynomial_system()
    sage: gen = BBGenerator(optimization='optimistic')
    sage: basis,orderIdeal,statistics = gen.calc_basis(F,2)

    sage: basis
     [k003*k200,
      k003*k201,
      k003*k202,
      k003*k203 + k003,
      k003*x200,
      k003*x201,
      k003*x202 + k003,
      k003*x203,
      k003*w200 + k003,
      k003*w201,
      k003*w202,
      k003*w203 + k003,
      k003*s100 + k003,
      k003*s101,
      k003*s102 + k003,
      k003*s103,
      k003*k100 + k003,
      k003*k101 + k003,
      k003*k102,
      k003*k103,
      k003*x100 + k003,
      k003*x101,
      k003*x102,
      k003*x103 + k003,
      k003*w100,
      k003*w101,
      k003*w102 + k003,
      k003*w103,
      k003*s000,
      k003*s001,
      k003*s002 + k003,
      k003*s003,
      k000*k003 + k003,
      k001*k003,
      k002*k003,
      k200,
      k201,
      k003 + k202 + 1,
      k003 + k203,
      k003 + x200 + 1,
      x201,
      x202 + 1,
      k003 + x203 + 1,
      k003 + w200,
      k003 + w201 + 1,
      w202,
      w203 + 1,
      k003 + s100,
      s101,
      s102 + 1,
      k003 + s103 + 1,
      k100 + 1,
      k101 + 1,
      k003 + k102 + 1,
      k103,
      k003 + x100,
      x101,
      x102,
      x103 + 1,
      w100,
      w101,
      k003 + w102,
      k003 + w103 + 1,
      k003 + s000 + 1,
      k003 + s001 + 1,
      k003 + s002,
      s003,
      k000 + 1,
      k001,
      k002 + k003 + 1]
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
    def __init__(self, optimization="enhanced", use_positions=True, use_matrix=True, use_autoreduction=True, use_pol_exclusion=False, order="deglex",min_mutants_limit=0):
        r"""
        Generates a ``BBGenerator`` and initializes it with the chosen optimization level

        INPUT::

            - ``optimization`` -- (default 'enhanced') which algorithm variant should be used
            - ``use_positions`` -- (default True) whether to use calculated DegLex/DegRevLex-Positions in the algorithm  
            - ``use_matrix`` -- (default True) whether to use a matrix during the calculation or a own gauss-variant instead
            - ``use_autoreduction`` (default True) if caluclations are executed in GF(2), defines whether the algorithm should autoreduce exponents >1 to 1
            - ``use_pol_exclusion`` (default False) whether minimal polynomials should affect the computational universe in future rounds
            - ``order`` (default 'deglex') which term ordering to use during the calculation
            - ``min_mutants_limit`` (default 0) percentage 0.0-1.0, that describes how many new polynomials are extended at once in the mutant algorithms. The final value is calculated as current_amount_of_leading_terms * min_mutants_limit / number_of_variables.      

        EXAMPLES::

            sage: from sage.borderbasis.generator import BBGenerator

            sage: BBGenerator()
            BBGenerator(optimization='enhanced')
            sage: BBGenerator(optimization='improved_mutant',order='degrevlex')
            BBGenerator(optimization='improved_mutant')
        """
        self.order = order
        self.optimization = optimization
        self.use_positions = use_positions
        self.use_matrix = use_matrix
        self.use_autoreduction = use_autoreduction
        self.use_pol_exclusion = use_pol_exclusion
        self.min_mutants_limit = min_mutants_limit

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
            - ``reduce_monomials_but`` (default None) if this value is a list of variables, the algorithm performs a (slow) quadratic substitution of the variables not provided in the list (as far as possible)
            - ``keep_only`` (default None) if this value is a list of variables, the algorithm substitutes variables not in the list with linear substitution during the calculation.

        OUTPUT::

            A 3-touple consisting of objects in this order:

            1. The calculated borderbasis as list of polynomials
            2. The according order ideal as polynomial
            3. A map containing statistics and benchmarks, currently either the biggest matrix, or the maximum amount of comparisons in one run is logged (depending on the configuration of the object)

        EXAMPLES:

            sage: from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence
            sage: from sage.borderbasis.generator import BBGenerator

            sage: R.<x,y> = PolynomialRing(GF(2),2)
            sage: F = PolynomialSequence([x*y,y**2+x],R)
            sage: gen = BBGenerator(optimization='none',use_autoreduction=False)
            sage: basis,orderIdeal,statistics = gen.calc_basis(generators=F,modPolynomial=2)

            sage: basis
            [y^2 + x, x*y, x^2]
            sage: orderIdeal
            x + y + 1
            sage: statistics
            {'maxComparisons': 0L, 'maxMatrix': {'columns': 7L, 'rows': 6L}}
            
            sage: gen = BBGenerator(optimization='experimental',use_autoreduction=False,use_matrix=False)
            sage: F = PolynomialSequence([x**2,y**2+1],R)
            sage: basis,orderIdeal,statistics = gen.calc_basis(F,2)
            sage: basis
            [y^2 + 1, x^2, x*y^2 + x, x^2*y]
            sage: statistics
            {'maxComparisons': 2L, 'maxMatrix': {'columns': 0L, 'rows': 0L}}

        WARNING::

            Currently, it is only possible to calculate border bases of polynomials in the galois field.
        """
        use_variable_exclusion = False
        if(keep_only != None):
            use_variable_exclusion = True

            variables = self._get_variables(generators)
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
        bbt = PyBorderBasisTools_uint64(field,matrix,polynomialFactory,monFactory,generators.nvariables(),self.optimization,self.use_pol_exclusion,use_variable_exclusion,keep_only,modPolynomial==2 and self.use_autoreduction,self.min_mutants_limit)

        basis,orderIdeal = bbt.calculate_basis(generators)
        statistics = bbt.get_statistics()
        
        del field
        del polynomialFactory
        del monFactory
        del bbt

        return (basis,orderIdeal,statistics)

    def shrink_system(self,plist,keyvars,polybori=False,max_degree=0x7fffffff):
        r"""
        Shrinks the provided polynomial system through quadratic substitution

        INPUT::
            - ``plist`` -- a PolynomialSequence containing the polynomials that are to be shrinked
            - ``keyvars`` -- a list of variables that should not be substituted
            - ``polybori`` (default False) whether polynomials sare in GF(2) and should be autoreduced (if exponents >1 are automatically 1)
            - ``max_degree`` (default 0x7fffffff) the degree that the resulting polynomials are allowed to have. If above, the substitution is reversed for this polynomial.

        OUTPUT::
            
            A PolynomialSequence containing the shrinked polynomials
        
        EXAMPLES:
        
            sage: from sage.borderbasis.generator import BBGenerator

            sage: sr = mq.SR(2,1,1,4,gf2=True,polybori=False)
            sage: F,s = sr.polynomial_system()
            sage: gen = BBGenerator('none')
            sage: key_vars = F.ring().gens()[len(F.ring().gens())-4:len(F.ring().gens())]
            sage: G = gen.shrink_system(F,key_vars,polybori=True)
            sage: F
            Polynomial Sequence with 104 Polynomials in 36 Variables
            sage: G
            Polynomial Sequence with 48 Polynomials in 16 Variables
        """
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
                 try:
                     for k in range(0,len(plist[index].monomials())):
                         m = plist[index].monomials()[k]
                         try:
                             reducible = m.reducible_by(v.monomials()[0])
                         except:
                             reducible = v.monomials()[0].divides(m)
                         if(reducible):
                             factor = plist[index].coefficients()[k]
                             break
                     B = plist[index]/factor
                     B = v - B
                 except:
                     # when this is a boolean polynomial
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
        r"""
        Reduce all appearences of A with B in plist, up to a result of degree max_degree
        """
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
        r"""
        returns an index for a polynomial in plist that enables the substitution of the variable var, or -1 if none such can be found.
        """
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

    def _get_variables(self,pythonList):
        r"""
        returns a list of all actually appearing variables in pythonList, ordered by the variable ordering.
        """
        existing = pythonList.variables()
        ordered = pythonList.ring().gens()
        result = []
        for i in ordered:
            for e in existing:
                if ("%s" % (i)) == ("%s" % (e)):
                    result.append(i)
                    break
        return result


