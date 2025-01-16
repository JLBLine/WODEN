import logging
import logging.handlers
from ctypes import CDLL, CFUNCTYPE, c_char_p
import ctypes
import importlib_resources
import wodenpy

class MultiLineFormatter(logging.Formatter):
    """Multi-line formatter for logging messages. Means that the indentation
    comes out inline with end of time and level message, if you happen submit
    a long multi-line string to a single log call.
    
    Taken from https://stackoverflow.com/questions/58590731/how-to-indent-multiline-message-printed-by-python-logger"""
    def get_header_length(self, record):
        """Get the header length of a given record."""
        return len(super().format(logging.LogRecord(
            name=record.name,
            level=record.levelno,
            pathname=record.pathname,
            lineno=record.lineno,
            msg='', args=(), exc_info=None
        )))

    def format(self, record):
        """Format a record with added indentation."""
        indent = ' ' * self.get_header_length(record)
        head, *trailing = super().format(record).splitlines(True)
        return head + ''.join(indent + line for line in trailing)


# Main logging configuration
def main_logging_config(queue, gitlabel, logging_level = logging.DEBUG, log_file = False):
    # Configure the main logger to consume messages from the queue
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging_level)
    
    # formatter = MultiLineFormatter("%(message)s")
    formatter = MultiLineFormatter("%(asctime)s - %(levelname)s - %(message)s",
                                          '%Y-%m-%d %H:%M:%S')
    
    stream_handler.setFormatter(formatter)
    handles = [stream_handler]
    
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging_level)
        file_handler.setFormatter(formatter)
        handles.append(file_handler)
        
    listener = logging.handlers.QueueListener(queue, *handles)
        
    listener.start()
    
    logger = get_logger_from_queue(queue)
    
    logo_string = rf"""
                 )  (              )  
     (  (     ( /(  )\ )        ( /(  
     )\))(   ')\())(()/(   (    )\()) 
    ((_)()\ )((_)\  /(_))  )\  ((_)\  
    _(())\_)() ((_)(_))_  ((_)  _((_) 
    \ \((_)/ // _ \ |   \ | __|| \| | 
     \ \/\/ /| (_) || |) || _| | .` | 
      \_/\_/  \___/ |___/ |___||_|\_| 
      
    You are using wodenpy version/git hash: {gitlabel}
    """
    
    logger.info(logo_string)
    
    return listener, handles

def get_logger_from_queue(queue, logging_level = logging.DEBUG):
    # Configure the worker logger to send messages to the queue
    handler = logging.handlers.QueueHandler(queue)
    root_logger = logging.getLogger()
    root_logger.setLevel(logging_level)
    root_logger.handlers.clear()  # Clear existing handlers to avoid duplication
    root_logger.addHandler(handler)
    
    return logging.getLogger(__name__)


def get_log_callback(logger, logging_level = logging.DEBUG):
    LogCallbackType = CFUNCTYPE(None, c_char_p)

    # Define the Python logging function for libwoden C library
    def log_callback(message):
        if logging_level == logging.DEBUG:
            logger.debug(f"libwoden: {message.decode('utf-8')}")
        else:
            logger.info(f"libwoden: {message.decode('utf-8')}")

    # Wrap the Python function as a C callback
    c_log_callback = LogCallbackType(log_callback)
    
    return c_log_callback

def simple_logger(logging_level = logging.DEBUG):
    """Use this to set a default logger for wodenpy functions when
    the user doesn't provide one. Basically useful for unit tests."""
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=logging_level,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    return logger