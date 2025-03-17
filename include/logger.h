#include <stdio.h>

#ifndef __logger_h__
#define __logger_h__

#ifdef __cplusplus
  extern "C" {
#endif

/**
 @typedef log_callback_t
 @brief A function pointer type for log callback functions.

 This type defines a function pointer for callback functions that
 handle log messages. The callback function takes a single parameter,
 which is a constant character pointer to the log message.

 Initially sets the callback function to NULL. If the callback function
 is NULL, then `log message` defaults to logging all messages prefaced by
 "libwoden logger: " to stdout.

 @param[in] message A constant character pointer to the log message.
*/
typedef void (*log_callback_t)(const char *message);


/**
 * @brief Sets the callback function for logging.
 *
 * This function allows you to set a custom callback function that will be
 * used for logging messages. The callback function should match the 
 * signature defined by `log_callback_t`. Can be used from Python to
 * redirect logging to a Python logger.
 *
 * @param callback The callback function to be used for logging.
 */
void set_log_callback(log_callback_t callback);


/**
 * @brief Logs a message to the logging system.
 *
 * This function takes a string message as input and logs it to the 
 * appropriate logging system. The exact behavior of the logging 
 * if controlled by `set_log_callback`.
 *
 * @param message The message to be logged. It should be a null-terminated 
 *                string.
 */
void log_message(const char *message);

#ifdef __cplusplus
}
#endif

#endif // __logger_h__