from sys import path
import os
import unittest
import numpy as np
from time import sleep
import logging
import logging.handlers
from ctypes import CDLL, CFUNCTYPE, c_char_p
import ctypes
from wodenpy.wodenpy_setup.woden_logger import get_logger_from_queue, get_log_callback, listener_process, listener_configurer
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager, set_start_method
from multiprocessing import Process, Queue
import multiprocessing
import concurrent.futures

code_dir = os.path.realpath(__file__)
code_dir = ('/').join(code_dir.split('/')[:-1])

##annoying path hack to find where the C library is
test_dir = code_dir + "/../../../build/cmake_testing/wodenpy/wodenpy_setup/"

testing_lib = f"{test_dir}/libc_logger_funcs.so"

class fake_args:
    def __init__(self):
        self.precision = "double"
        self.log_file = None
        
# Worker function
def python_worker_function(n, queue):
    logger = get_logger_from_queue(queue, logging.DEBUG)
    logger.debug(f"Thread number: {n}")
    
    return n**2, n


def get_check_the_logger(loaded_lib):
    check_the_logger = loaded_lib.check_the_logger
    check_the_logger.argtypes = [ctypes.c_int]
    return check_the_logger

def libwoden_worker_function(n, inputs, queue, testing_lib):
    logger = get_logger_from_queue(queue, logging.DEBUG)
    
    c_log_callback = get_log_callback(logger)
    woden_lib = ctypes.cdll.LoadLibrary(testing_lib)
    woden_lib.set_log_callback(c_log_callback)
    
    check_the_logger = get_check_the_logger(woden_lib)
    check_the_logger(inputs[n])
    logger.info(f"C thread {n} was fed {inputs[n]}")
    
    return 


def run_cpu_style_parallel(log_queue, testing_lib, num_threads,
                           outputs):
    
    logger = get_logger_from_queue(log_queue, logging.DEBUG)
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as python_exec, concurrent.futures.ProcessPoolExecutor(max_workers=num_threads) as c_exec:
            
        for round_num in range(2):
            logger.info(f"Reading set {round_num} things")
        
            python_data = [python_exec.submit(python_worker_function, i, log_queue)
                                for i in range(num_threads)]
            results_python = [future_python.result() for future_python in python_data]
            
            all_results_python = [0]*num_threads
            # all_orders_python = []
            
            for future_python in as_completed(python_data):
                result, order = future_python.result()
                all_results_python[order] = result
                # all_orders_python.append(order)
                
            c_data = [c_exec.submit(libwoden_worker_function, 
                                    i, all_results_python, log_queue,
                                    testing_lib)
                                for i in range(num_threads)]
            results = [future.result() for future in c_data]
            for future in as_completed(c_data):
                future.result()
                
    return outputs

##Vehicle for running tests
class Test(unittest.TestCase):
    """Test whether the args collected by the argument parser are read in
    correctly, and that sanity checks on certain combinations work such that
    we don't feed WODEN arguments that won't work"""

        
    def test_make_and_run_logger_multithread(self):
        """Check we can create the logger and call it from both Python and C
        when doing multithreading. This function is basically a paired down
        version of what happens in run_woden.py"""
        
        
        # for meh in range(10):
        args = fake_args()
        args.log_file = 'test.log'
        
        num_threads = 3
        
        ##RIGHTO, so unittest throws in some of it's own threading (I think)
        ##which will throw an error, unless we explicitly set the start method
        try:
            set_start_method("forkserver")  # Alternatives: "spawn" or "forkserver"
        except RuntimeError:
            # Start method is already set in some environments like Jupyter
            pass
        
        log_queue = multiprocessing.Manager().Queue(-1)
        listener = multiprocessing.Process(target=listener_process,
                                        args=(log_queue, listener_configurer,
                                              args.log_file))
        listener.start()
            
        main_logger = get_logger_from_queue(log_queue, logging.INFO)
        
        main_logger.info("Before we be doing parallel")
        main_logger.warning("He's not the messiah, he's a very naughty boy")
        main_logger.error("RUN AWAY!")
        
        libwoden_worker_function(0, [10], log_queue, testing_lib)
        # sleep(1.0)
        
        outputs = []
        
        outputs = run_cpu_style_parallel(log_queue, testing_lib, num_threads,
                        outputs)

                
        # sleep(0.5)
        main_logger.info("After we be doing parallel")
        
        log_queue.put_nowait(None)
        listener.join()
                    
        outputs = []
        ##The start of each line is a time, so strip that out to check
        ##the outputs against stuff we can predict
        ##We also have a silly WODEN logo at the start of the log
        ##Just skip testing that as it doesn't really matter what it looks
        ##like, we care about the actual information     
        with open('test.log', 'r') as f:
            lines = f.read().split('\n')
            
            first_line = False
            
            for line in lines:
                if "Before we be" in line:
                    first_line = True
                    
                if first_line:
                
                    bits = line.split(' - ')
                    
                    if len(bits) == 3:
                        outputs.append(' - '.join(bits[1:]))
            
        ##First four things are run in parallel, so we know the order
        expected_start = ["INFO - Before we be doing parallel",
                          "WARNING - He's not the messiah, he's a very naughty boy",
                          "ERROR - RUN AWAY!",
                          "DEBUG - libwoden: We be running in C and the number is 10",
                          ]
            
        self.assertEqual(outputs[:4], expected_start)
        
        ##The next three things are run in parallel, so might plop
        ##out in any order. Just assert they exist in the outputs
        
        expected_parallel = ["DEBUG - Thread number: 0",
                             "DEBUG - Thread number: 1",
                             "DEBUG - Thread number: 2",
                             "INFO - C thread 0 was fed 0",
                             "INFO - C thread 1 was fed 1",
                             "INFO - C thread 2 was fed 4",
                             "DEBUG - libwoden: We be running in C and the number is 0",
                             "DEBUG - libwoden: We be running in C and the number is 1",
                             "DEBUG - libwoden: We be running in C and the number is 4",
                             ]
        
        for line in expected_parallel:
            self.assertIn(line, outputs)
            
        ##This should be the last thing to be logged
        self.assertEqual("INFO - After we be doing parallel", outputs[-1])
            

##Run the test
if __name__ == '__main__':
    unittest.main()