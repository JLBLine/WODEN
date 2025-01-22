#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "logger.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void fancier_woden_logger(const char *message) {
  printf("OH LOOKY AT ME: %s\n", message);
}

/*
Just call the logger to make sure it doesn't crash
*/
void test_woden_logger(void) {
  log_message("This is default log message");

  //Now set the callback to the fancier one
  //This is functionality to set the callback to the Python logger.
  //Test it here for the sake of the coverage report
  set_log_callback(fancier_woden_logger);
  log_message("This is a fancier log message");
  
}



//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_woden_logger);

    return UNITY_END();
}
