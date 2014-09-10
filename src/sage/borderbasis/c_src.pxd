cdef extern from "stdint.h":
    ctypedef unsigned long long uint64_t

cdef extern from "src/include/i_owningList.h" namespace "base":
    cdef cppclass IOwningList[T]:
        T& at(unsigned int pos)
        void clear()
        void push_back(T val)
        unsigned int size()

cdef extern from "src/include/i_monomial.h" namespace "polynomial":
    cdef cppclass IMonomial:
        void set(unsigned int index, unsigned int value)
        unsigned int& at(unsigned int & index)
        unsigned int getIndet()

cdef extern from "src/include/term.h" namespace "polynomial":
    cdef cppclass Term[T]:
        Term(T coef, IMonomial* monomial) except +
        Term() except +
        T getCoef() 
        void setCoef(T coef)
        IMonomial* getMonomial()
        void setMonomial(IMonomial* monomial)

cdef extern from "src/include/statistics.h" namespace "borderbasis":
    cdef struct MatrixProp:
        uint64_t rows
        uint64_t columns
    cdef cppclass Statistics:
        Statistics() except +
        uint64_t max_comparisons_in_reduction
        MatrixProp maxMatrix

cdef extern from "src/include/monomialFactory.h" namespace "polynomial":
    cdef cppclass IMonomialFactory:
        IMonomialFactory() except + 
        IMonomial* create()
    cdef cppclass MonomialFactoryNoOrderPos:
        MonomialFactoryNoOrderPos(unsigned int indet) except +
        IMonomial* create()
    cdef cppclass MonomialFactoryDegLex:
        MonomialFactoryDegLex(unsigned int indet) except +
        IMonomial* create()
    cdef cppclass MonomialFactoryNoOrderPosGF2:
        MonomialFactoryNoOrderPosGF2(unsigned int indet) except +
        IMonomial* create()
    cdef cppclass MonomialFactoryDegLexGF2:
        MonomialFactoryDegLexGF2(unsigned int indet) except +
        IMonomial* create()

cdef extern from "src/include/i_polynomial.h" namespace "polynomial":
    cdef cppclass IPolynomial_uint64:
        void push(Term[uint64_t]* term)
        unsigned int size()
        void clear()
        Term[uint64_t]* at(unsigned int index)
    cdef cppclass IPolynomial[T]:
        void push(Term[T]* term)
        unsigned int size()
        void clear()
        Term[T]* at(unsigned int index)

cdef extern from "src/include/owningVector.h" namespace "base":
    cdef cppclass OwningVector[T]:
        OwningVector() except +
        OwningVector(unsigned int n) except +
        T& at(unsigned int pos)
        void clear()
        void push_back(T val)
        unsigned int size()

cdef extern from "src/include/polynomialFactory.h" namespace "polynomial":
    cdef enum PolType:
        POLTYPE_VECTOR
        POLTYPE_VECTOR_GF2
    cdef cppclass PolynomialFactory[T]:
        PolynomialFactory(PolType t) except +
        IPolynomial[T]* create(unsigned int indet) 

cdef extern from "src/include/field.h" namespace "math":
    cdef cppclass IField[T]:
        IField() except +
    cdef cppclass FieldFn:
        FieldFn(uint64_t pol) except +

cdef extern from "src/include/i_matrixFactory.h" namespace "math":
    cdef cppclass IMatrixFactory[T]:
        IMatrixFactory() except +

cdef extern from "src/include/matrixFactory_fn.h" namespace "math":
    cdef cppclass MatrixFactory_Fn:
        MatrixFactory_Fn(uint64_t minPolynomial) except +

cdef extern from "src/include/borderBasisTools.h" namespace "borderbasis":
    cdef enum OptLevel:
        NONE
        ENHANCED
        OPTIMISTIC
        EXPERIMENTAL
        MUTANT
        IMPROVED_MUTANT
        IMPROVED_MUTANT_LINEAR
        IMPROVED_MUTANT_OPTIMISTIC
    cdef cppclass BorderBasisTools[T]:
        BorderBasisTools(IField[T]* field,
                         PolynomialFactory[T]* polFactory,
                         IMonomialFactory* monFactory,
                         unsigned int indeterminates,
                         OptLevel optimizations) except +
        BorderBasisTools(int dummyNecessaryForPython,
                         IMatrixFactory[T]* matrixFactory,
                         PolynomialFactory[T]* polFactory,
                         IMonomialFactory* monFactory,
                         unsigned int indeterminates,
                         OptLevel optimizations) except +
        void calculateBasis(void* i,
                         void* o,
                         void* orderIdeal)
        void getStatistics(Statistics* o) 
        
