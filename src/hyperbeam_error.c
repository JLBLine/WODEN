#include "hyperbeam_error.h"

void handle_hyperbeam_error(const char file[], int line_num, const char function_name[]) {
    int err_length = hb_last_error_length();
    char *err = (char *)malloc(err_length * sizeof(char));
    int err_status = hb_last_error_message(err, err_length);
    if (err_status == -1) {
        log_message("Something really bad happened with hyperbeam");
        exit(EXIT_FAILURE);
    }
    //Print error message to stderr
    fprintf(stderr,"File %s:%d: hyperbeam error in %s: %s",
            file, line_num, function_name, err);

    //Also log the error message
    char log_buffer[1024];
    int log_len = sizeof log_buffer;
    snprintf(log_buffer, log_len,
             "File %s:%d: hyperbeam error in %s: %s", file, line_num, function_name, err);
    log_message(log_buffer);
    //Do an exit
    exit(EXIT_FAILURE);
}