#include <limits.h>
.globl fixedAddNoINTMIN
fixedAddNoINTMIN:
#define a r0
#define b r1
#define sum r0
#define result r0

	add sum, a, b, nov, return_sum
	jlts b, 0, ret_min

// ret_max:
	move result, INT_MAX
	jump r23

ret_min:
	move result, INT_MIN
	jump r23

return_sum:
	jump r23
