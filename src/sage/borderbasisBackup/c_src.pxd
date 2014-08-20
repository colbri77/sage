cdef extern from "src/include/i_owningList.h" namespace "base":
    cdef cppclass IOwningList[T]:
        T& at(unsigned int pos)
        void clear()
        void push_back(T val)
        unsigned int size()

cdef extern from "src/include/i_monomial.h" namespace "polynomial":
    cdef cppclass IMonomial[T]:
        void set(unsigned int index, unsigned int value)
        unsigned int& at(unsigned int & index)
        unsigned int getIndet()

cdef extern from "src/include/term.h" namespace "polynomial":
    cdef cppclass Term[T]:
        Term(T coef, IMonomial[T]* monomial) except +
        Term() except +
        T getCoef() 
        void setCoef(T coef)
        IMonomial[T]* getMonomial()
        void setMonomial(IMonomial[T]* monomial)

cdef extern from "src/include/statistics.h" namespace "borderbasis":
    cdef struct MatrixProp:
        unsigned int rows
        unsigned int columns
    cdef cppclass Statistics:
        Statistics() except +
        MatrixProp maxMatrix

cdef extern from "src/include/monomialFactory.h" namespace "polynomial":
    cdef enum MonomialType:
        MONOMIALTYPE_DEGLEX
    cdef cppclass MonomialFactory[T]:
        MonomialFactory(MonomialType t) except + 
        IMonomial[T]* create(unsigned int indet) 

cdef extern from "src/include/i_polynomial.h" namespace "polynomial":
    cdef cppclass IPolynomial_uint64:
        void push(Term[long long unsigned int]* term)
        unsigned int size()
        void clear()
        Term[long long unsigned int]* at(unsigned int index)
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
    cdef cppclass PolynomialFactory[T]:
        PolynomialFactory(PolType t) except +
        IPolynomial[T]* create(unsigned int indet) 

cdef extern from "src/include/i_matrixFactory.h" namespace "matrix":
    cdef cppclass IMatrixFactory[T]:
        IMatrixFactory() except +

cdef extern from "src/include/matrixFactory_fn.h" namespace "matrix":
    cdef cppclass MatrixFactory_Fn:
        MatrixFactory_Fn(long long unsigned int minPolynomial)

cdef extern from "src/include/borderBasisTools.h" namespace "borderbasis":
    cdef enum OptLevel:
        NONE
        ENHANCED
        OPTIMISTIC
        EXPERIMENTAL
    cdef cppclass BorderBasisTools[T]:
        BorderBasisTools(IMatrixFactory[T]* matrixFactory, 
                         PolynomialFactory[T]* polFactory,
                         MonomialFactory[T]* monFactory,
                         unsigned int indeterminates,
                         OptLevel optimizations)
        void calculateBasis(void* i,
                         void* o,
                         void* orderIdeal)
        void getStatistics(Statistics* o) 
        
