#include "logger.h"

// Default to printing to stdout
static log_callback_t current_log_callback = NULL;

void set_log_callback(log_callback_t callback) {
  current_log_callback = callback;
}

void log_message(const char *message) {
  if (current_log_callback) {
      current_log_callback(message);
  } else {
      printf("libwoden logger: %s\n", message);
  }
}