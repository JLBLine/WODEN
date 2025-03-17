#include "check_compilation_flags.h"

bool check_for_everybeam_compilation(void) {
  #ifdef HAVE_EVERYBEAM
    return true;
  #else
    return false;
  #endif
}