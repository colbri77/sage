#include "include/statistics.h"

namespace borderbasis {

Statistics::Statistics()
: max_comparisons_in_reduction(0)
{

}

Statistics::~Statistics()
{

}

void Statistics::start()
{
    max_comparisons_in_reduction = 0;
    maxMatrix.rows = 0;
    maxMatrix.columns = 0;
}

void Statistics::stop()
{

}

void Statistics::logMatrix(uint64_t r,uint64_t c)
{
    if(maxMatrix.rows*maxMatrix.columns < r*c) {
        maxMatrix.rows = r;
        maxMatrix.columns = c;
    }
}

} // namespace borderbasis
