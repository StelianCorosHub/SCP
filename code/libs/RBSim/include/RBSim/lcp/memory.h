#pragma once

/* round an integer up to a multiple of 4, except that 0 and 1 are unmodified
 * (used to compute matrix leading dimensions)
 */
#define dPAD(a) (((a) > 1) ? (((a) + 3) & (int)(~3)) : (a))

#define EFFICIENT_ALIGNMENT 16

/* round something up to be a multiple of the EFFICIENT_ALIGNMENT */
#define dEFFICIENT_SIZE(x) ((((x)-1)|(EFFICIENT_ALIGNMENT-1))+1)

/* alloca aligned to the EFFICIENT_ALIGNMENT. note that this can waste
 * up to 15 bytes per allocation, depending on what alloca() returns.
 */
#define dALLOCA16(n) \
  ((char*)dEFFICIENT_SIZE(((size_t)(alloca((n)+(EFFICIENT_ALIGNMENT-1))))))

#define ALLOCA dALLOCA16
#define STACK_ALLOC_MAX 8192U



