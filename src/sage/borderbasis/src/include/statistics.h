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
    void logMatrix(uint64_t r,uint64_t c);

    uint64_t max_comparisons_in_reduction;
    MatrixProp maxMatrix;
};

} // namespace borderbasis

#endif // __STATISTICS_H__
