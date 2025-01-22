#include <unity.h>
#include "logger.h"
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>
#include <string.h>

extern void try_to_memcpy_to_null_pointer_gpu(double *cpu_array, int num_elements);

// static jmp_buf jump_buffer; // Define the jump buffer.
// static int mock_exit_called = 0; // Track whether exit was called.
// static int mock_exit_status = 0;

// void exit(int status) {
//     mock_exit_called = 1;
//     mock_exit_status = status;
//     longjmp(jump_buffer, 1); // Use setjmp/longjmp to bypass actual exit.
// }



// #include "hyperbeam_error.h"

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
void test_log_gpu_errors(void) {

  const char *file = "test_log.txt";

  // Attempt to delete the file
  if (remove(file) == 0) {
      printf("Old test_log.txt deleted successfully.\n");
  } else {
      printf("No old test_log.txt found\n");
  }

  set_log_callback(write_log_to_text);
  log_message("Does this work?");

  try_to_memcpy_to_null_pointer_gpu(NULL, 10);
  
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
              strstr(buffer, "GPU ERROR cudaMemcpy"),
              "Did not find 'GPU ERROR cudaMemcpy' in log message, so error message is missing"
          );
          TEST_MESSAGE("Found 'GPU ERROR cudaMemcpy' in log file, so error message was captured.");
      }
  }
 
 
  // }
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_log_gpu_errors);

    return UNITY_END();
}