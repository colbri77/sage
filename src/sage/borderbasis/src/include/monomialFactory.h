#ifndef __MONOMIALFACTORY_H__
#define __MONOMIALFACTORY_H__

#include "i_monomial.h"

namespace polynomial {

enum MonomialType {
	MONOMIALTYPE_DEGLEX,
	MONOMIALTYPE_DEGLEX_NO_ORDER_POS
};

class MonomialFactory
{
public:
    MonomialFactory(MonomialType type);
    virtual ~MonomialFactory();

    virtual bool supportsGetPos() const;

    virtual TAKE_OWN IMonomial* create(uint indet) const;
    virtual TAKE_OWN IMonomial* create(uint64_t pos, uint indet) const;

private:
    MonomialType type;
};

} // namespace polynomial

#endif // __MONOMIALFACTORY_H__
