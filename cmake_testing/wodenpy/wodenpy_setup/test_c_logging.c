#include <stdio.h>
#include <stdlib.h>
#include "logger.h"

void check_the_logger(int number)
{

  char msg_buffer[50];
  
  int max_len = sizeof msg_buffer;

  int j = snprintf(msg_buffer, max_len,
                "We be running in C and the number is %d", number);

  if (j >= max_len) {
      fputs("Buffer length exceeded; string truncated", stderr);
  }

  log_message(msg_buffer);
    
}
