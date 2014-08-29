#ifndef __MONOMIALFACTORY_H__
#define __MONOMIALFACTORY_H__

#include "i_monomial.h"
#include "fastFlexibleArray.h"

namespace polynomial {

class IMonomialFactory
{
public:
    IMonomialFactory(){}
    virtual ~IMonomialFactory(){}

    virtual bool supportsGetPos() const = 0;
    virtual TAKE_OWN IMonomial* create() const = 0;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const = 0;
};

class MonomialFactoryNoOrderPos : public IMonomialFactory
{
public:
    MonomialFactoryNoOrderPos(uint indet);
    virtual ~MonomialFactoryNoOrderPos();

    virtual bool supportsGetPos() const OVERRIDE;

    virtual TAKE_OWN IMonomial* create() const OVERRIDE;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;

private:
    uint indet;
};

class MonomialFactoryDegLex : public IMonomialFactory
{
public:
    MonomialFactoryDegLex(uint indet);
    virtual ~MonomialFactoryDegLex();

    virtual bool supportsGetPos() const OVERRIDE;

    virtual TAKE_OWN IMonomial* create() const OVERRIDE;
    virtual TAKE_OWN IMonomial* create(uint64_t pos) const OVERRIDE;

private:
    FastFlexibleArray* monomials;
    uint indet;
};

} // namespace polynomial

#endif // __MONOMIALFACTORY_H__
