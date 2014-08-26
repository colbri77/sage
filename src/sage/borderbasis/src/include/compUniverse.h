#ifndef __COMPUNIVERSE_H__
#define __COMPUNIVERSE_H__

#include "owningVector.h"
#include "polynomial.h"

namespace borderbasis {

template<typename T>
class ICompUniverse
{
public:
    ICompUniverse(uint indet);
    virtual ~ICompUniverse();

    virtual void clear() = 0;
    virtual void extend(uint64_t limitDegree) = 0;
    virtual uint getMaxDegree() const = 0;
    virtual uint64_t getMaxPos() const = 0;
    virtual bool contains(uint64_t pos) const = 0;
    virtual void add(uint64_t pos) = 0;
    virtual void addBorder() = 0;

    virtual bool contains(IPolynomial<T>* pol) const;
    virtual bool contains(IMonomial* monomial) const;
    virtual void add(IOwningList<IPolynomial<T>*>* additions,uint start);
    virtual void add(IOwningList<IPolynomial<T>*>* additions);
    virtual void add(IMonomial* monomial);

protected:
    uint indet;
};

template<typename T>
class LinearCompUniverse : public ICompUniverse<T>
{
public:
    LinearCompUniverse(uint indet);
    virtual ~LinearCompUniverse();

    virtual void clear() OVERRIDE;
    virtual void extend(uint64_t limitDegree) OVERRIDE;
    virtual uint getMaxDegree() const OVERRIDE;
    virtual uint64_t getMaxPos() const OVERRIDE;
    virtual bool contains(uint64_t pos) const OVERRIDE;
    virtual bool contains(IMonomial* monomial) const OVERRIDE;
    virtual void add(uint64_t pos) OVERRIDE;
    virtual void addBorder() OVERRIDE;
    virtual void add(IOwningList<IPolynomial<T>*>* additions,uint start) OVERRIDE;

private:
    uint limit;
};

template<typename T>
class SpecificCompUniverse : public ICompUniverse<T>
{
public:
    SpecificCompUniverse(uint indet);
    virtual ~SpecificCompUniverse();

    virtual void clear() OVERRIDE;
    virtual void extend(uint64_t limitDegree) OVERRIDE;
    virtual uint getMaxDegree() const OVERRIDE;
    virtual uint64_t getMaxPos() const OVERRIDE;
    virtual bool contains(uint64_t pos) const OVERRIDE;
    virtual bool contains(IMonomial* monomial) const OVERRIDE;
    virtual void add(uint64_t pos) OVERRIDE;
    virtual void addBorder() OVERRIDE;
    virtual void add(IMonomial* monomial) OVERRIDE;

private:
    Polynomial<T>* U;
};

} // namespace borderbasis

#endif // __COMPUNIVERSE_H__
