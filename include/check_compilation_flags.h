#include <stdio.h>
#include <stdbool.h>

#ifndef __check_comp_flags_h__
#define __check_comp_flags_h__

#ifdef __cplusplus
  extern "C" {
#endif

/**
 @brief Returns True if the code is compiled with the -DEVERYBEAM flag set, False otherwise.
 @return True if the code is compiled with the -DEVERYBEAM flag set, False otherwise.
*/
bool check_for_everybeam_compilation(void);

#ifdef __cplusplus
}
#endif

#endif // __check_comp_flags_h__