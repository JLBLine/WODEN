import logging
import logging.handlers
from ctypes import CDLL, CFUNCTYPE, c_char_p
import ctypes
import importlib_resources
import wodenpy
from multiprocessing import Process
from logging.handlers import QueueHandler, QueueListener

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
    
    
def listener_configurer(log_file):
    root = logging.getLogger()
    
    formatter = MultiLineFormatter("%(asctime)s - %(levelname)s - %(message)s",
                                   '%Y-%m-%d %H:%M:%S')
    
    stream_handler = logging.StreamHandler()
    # stream_handler.setLevel(logging_level)
    stream_handler.setFormatter(formatter)
    handles = [stream_handler]
    
    if log_file:
        file_handler = logging.FileHandler(log_file, mode='w')
        # file_handler.setLevel(logging_level)
        file_handler.setFormatter(formatter)
        handles.append(file_handler)
    
    for handle in handles:
        root.addHandler(handle)
    

# This is the listener process top-level loop: wait for logging events
# (LogRecords)on the queue and handle them, quit when you get a None for a
# LogRecord.
def listener_process(queue, configurer, log_file=False):
    configurer(log_file)
    while True:
        try:
            record = queue.get()
            if record is None:  # We send this as a sentinel to tell the listener to quit.
                break
            logger = logging.getLogger(record.name)
            logger.handle(record)  # No level or filter logic applied - just do it!
        except Exception as e:
            # import sys, traceback
            # print('Whoops! Problem:', file=sys.stderr)
            # traceback.print_exc(file=sys.stderr)
            
            exit(e)


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
    
    if queue == False:
        logger = simple_logger(logging_level)
    
    else:
        logger = logging.getLogger(__name__)
        logger.setLevel(logging_level)
        logger.handlers.clear()
        handler = QueueHandler(queue)
        logger.addHandler(handler)
    
    return logger
    
def set_logger_header(logger, gitlabel):
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