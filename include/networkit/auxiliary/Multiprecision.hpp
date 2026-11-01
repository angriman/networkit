/*
 * Multiprecision.hpp
 *
 *  Created on: 21.03.2013
 *      Author: cls
 */

#ifndef NETWORKIT_AUXILIARY_MULTIPRECISION_HPP_
#define NETWORKIT_AUXILIARY_MULTIPRECISION_HPP_

// GCC reports a false out-of-bounds access inside of ttmath's division.
#if defined(__GNUC__) && (__GNUC__ >= 12)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
#include <ttmath/ttmath.h>
#if defined(__GNUC__) && (__GNUC__ >= 12)
#pragma GCC diagnostic pop
#endif

namespace NetworKit {
using bigfloat = ttmath::Big<TTMATH_BITS(64), TTMATH_BITS(64)>; ///< big floating point number
}

#endif // NETWORKIT_AUXILIARY_MULTIPRECISION_HPP_
