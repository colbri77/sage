#ifndef __COMPUNIVERSE_H__
#define __COMPUNIVERSE_H__

#include "owningVector.h"
#include "i_polynomial.h"
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
    virtual void addBorder() = 0;
    virtual bool beyondLastElement(IMonomial* monomial) const = 0;
    virtual bool contains(IPolynomial<T>* pol) const;
    virtual bool contains(IMonomial* monomial) const;
    virtual void add(IOwningList<IPolynomial<T>*>* additions,uint start);
    virtual void add(IOwningList<IPolynomial<T>*>* additions);
    virtual void add(IMonomial* monomial);
    virtual bool exclude(IMonomial* monomial);

protected:
    uint indet;
    vector<IMonomial*>* exclusions;

    virtual void add(uint64_t pos) = 0;
    virtual bool contains(uint64_t pos) const = 0;
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
    virtual bool contains(IMonomial* monomial) const OVERRIDE;
    virtual void addBorder() OVERRIDE;
    virtual void add(IOwningList<IPolynomial<T>*>* additions,uint start) OVERRIDE;
    virtual bool beyondLastElement(IMonomial* monomial) const OVERRIDE;

protected:
    uint limit;

    virtual void add(uint64_t pos) OVERRIDE;
    virtual bool contains(uint64_t pos) const OVERRIDE;
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
    virtual void addBorder() OVERRIDE;
    virtual void add(IOwningList<IPolynomial<T>*>* additions,uint start) OVERRIDE;
    virtual bool beyondLastElement(IMonomial* monomial) const OVERRIDE;
    virtual void add(IMonomial* monomial) OVERRIDE;

protected:
    uint8_t* U;
    uint64_t uLen;
    uint uBlocks;
    uint maxDegree;
    vector<IMonomial*>* lastUBorderCandidates;

    virtual void add(uint64_t pos) OVERRIDE;
    virtual bool contains(uint64_t pos) const OVERRIDE;
};

template<typename T>
class SpecificCompUniverseNoBorderLog : public SpecificCompUniverse<T>
{
public:
    SpecificCompUniverseNoBorderLog(uint indet);
    virtual ~SpecificCompUniverseNoBorderLog();

    virtual void addBorder() OVERRIDE DONT_USE;
    virtual void add(IOwningList<IPolynomial<T>*>* additions,uint start) OVERRIDE;
};

template<typename T>
class SpecificCompUniverseNoOrderPos : public ICompUniverse<T>
{
public:
    SpecificCompUniverseNoOrderPos(uint indet);
    virtual ~SpecificCompUniverseNoOrderPos();

    virtual void clear() OVERRIDE;
    virtual void extend(uint64_t limitDegree) OVERRIDE;
    virtual uint getMaxDegree() const OVERRIDE;
    virtual bool contains(IMonomial* monomial) const OVERRIDE;
    virtual void addBorder() OVERRIDE;
    virtual void add(IMonomial* monomial) OVERRIDE;
    virtual bool beyondLastElement(IMonomial* monomial) const OVERRIDE;

private:
    Polynomial<T>* U;

    virtual void add(uint64_t pos) OVERRIDE;
    virtual bool contains(uint64_t pos) const OVERRIDE;
};


} // namespace borderbasis

#endif // __COMPUNIVERSE_H__
