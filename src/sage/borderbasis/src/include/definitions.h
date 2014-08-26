#ifndef __DEFINITIONS_H__
#define __DEFINITIONS_H__

#include <stdint.h>
#include <stdexcept>
//#include <cassert>

typedef unsigned int uint;

#define TAKE_OWN /*takes ownership of the reference*/
#define KEEP_REF /*keeps reference on this item, don't delete it during this classes lifetime*/
#define DONT_USE /*calling this method will result in invalid results or an exception*/

#define DEL_SAFE(x) if(x) delete (x)
#define ASSERT_NOT_REACHED throw runtime_error("assert(not reached) failed")
#define NOT_IMPLEMENTED throw runtime_error("not implemented")
#define ENSURE(x,y) if(!(x)) throw runtime_error(y)

#if __cplusplus <= 199711L
    #define OVERRIDE /*override*/
#else
    #define OVERRIDE override
#endif

// namespace definitions:
namespace base {using namespace std;}
namespace math {using namespace std;}
namespace polynomial {using namespace std;
                    using namespace base;
                    using namespace math;}
namespace borderbasis {using namespace std;
                    using namespace base;
                    using namespace polynomial;
                    using namespace math;}

#endif // __DEFINITIONS_H__
