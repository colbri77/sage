#ifndef __MONOMIALFACTORY_H__
#define __MONOMIALFACTORY_H__

#include "i_monomial.h"

namespace polynomial {

enum MonomialType {
	MONOMIALTYPE_DEGLEX
};

template<typename T>
class MonomialFactory
{
public:
    MonomialFactory(MonomialType type);
    virtual ~MonomialFactory();

    virtual TAKE_OWN IMonomial<T>* create(uint indet) const;
    virtual TAKE_OWN IMonomial<T>* create(uint64_t pos, uint indet) const;

private:
    MonomialType type;
};

} // namespace polynomial

#endif // __MONOMIALFACTORY_H__
