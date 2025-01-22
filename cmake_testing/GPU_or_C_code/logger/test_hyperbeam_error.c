#include <unity.h>
#include "logger.h"
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#include <string.h>

static jmp_buf jump_buffer; // Define the jump buffer.
static int mock_exit_called = 0; // Track whether exit was called.
static int mock_exit_status = 0;

void exit(int status) {
    mock_exit_called = 1;
    mock_exit_status = status;
    longjmp(jump_buffer, 1); // Use setjmp/longjmp to bypass actual exit.
}


#include "hyperbeam_error.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

void write_log_to_text(const char *message) {
  FILE *f = fopen("test_log.txt", "a");
  if (f == NULL) {
      printf("Error opening file!\n");
      exit(1);
  }

  fprintf(f, "log message: %s\n", message);
  fclose(f);
}


/*
Check that the hyperbeam error function works as intended for a normal error
failure of 1
*/
void test_hyperbeam_error(void) {

  const char *file = "test_log.txt";

  // Attempt to delete the file
  if (remove(file) == 0) {
      printf("Old test_log.txt deleted successfully.\n");
  } else {
      printf("No old test_log.txt found\n");
  }

  set_log_callback(write_log_to_text);
  log_message("Does this work?");

  if (setjmp(jump_buffer) == 0) {
    int status;
    struct FEEBeam *fee_beam;
    status = new_fee_beam("bad_hdf5_file_path", &fee_beam);

    if (status != 0) {
      handle_hyperbeam_error(__FILE__, __LINE__, "new_fee_beam");
    }
    log_message("WE BE HERE\n");

    TEST_FAIL_MESSAGE("exit was not called.");
  } else {
      // This block is executed after mock_exit calls longjmp.
      TEST_ASSERT_EQUAL_INT(1, mock_exit_called);
      TEST_ASSERT_EQUAL_INT(1, mock_exit_status); // Validate exit code.
      TEST_MESSAGE("exit was called as expected.");
 
      //read in the log file
      FILE *f = fopen("test_log.txt", "r");
      if (f == NULL) {
          printf("Error opening file!\n");
          exit(1);
      }

      char buffer[1024];
      int i = 0;
      while (fgets(buffer, 1024, f)) {
          
          i++;
          //Exit error should be in second log message
          if (i == 2) {

              TEST_ASSERT_NOT_NULL_MESSAGE(
                  strstr(buffer, "bad_hdf5_file_path"),
                  "Did not find bad_hdf5_file_path in log message, so error message is missing"
              );
              TEST_MESSAGE("Found bad_hdf5_file_path in log file, so error message was captured.");
          }
      }
  }
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_hyperbeam_error);

    return UNITY_END();
}
