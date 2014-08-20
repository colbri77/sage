#ifndef __STATISTICS_H__
#define __STATISTICS_H__

#include "definitions.h"

namespace borderbasis {

struct MatrixProp {
	uint rows;
	uint columns;
};

class Statistics {
public:
    Statistics();
    virtual ~Statistics();

    void start();
    void stop();
    void logMatrix(uint r,uint c);

    MatrixProp maxMatrix;
};

} // namespace borderbasis

#endif // __STATISTICS_H__
