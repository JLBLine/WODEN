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

 @param[in] message A constant character pointer to the log message.
*/
typedef void (*log_callback_t)(const char *message);

/**
 @brief Add a message to the log file

 @details TODO add how you call from Python once finalised

 @param[in] message The message to add to the log

*/
void set_log_callback(log_callback_t callback);

/**
@brief Add a message to the log file

@details TODO add how you call from Python once finalised

@param[in] *message The message to add to the log

*/
void log_message(const char *message);

#ifdef __cplusplus
}
#endif

#endif // __logger_h__