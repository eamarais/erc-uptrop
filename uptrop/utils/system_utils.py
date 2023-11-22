import datetime as dt
import os
import psutil
import gc
import numpy as np

def print_to_log(*args, **kwargs):
    """
    Prints the given arguments to the console and flushes the output to ensure it is immediately visible.
    This function is intended for use in long-running jobs where progress updates need to be logged in real-time.
    """
    print(*args, flush=True, **kwargs)
    
def get_current_datetime(as_str=True):
    """
    Returns the current date and time.
    Parameters:
    as_str (bool): If True, returns the date and time as a formatted string ('%Y-%m-%d %H:%M:%S').
                   Otherwise, returns a dt.datetime object.
    Returns:
    dt.datetime or str: The current date and time, as either a dt.datetime object or a formatted string.
    """
    if as_str:
        return dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    else:
        return dt.datetime.now()
    
def measure_run_time(func):
    """
    This is a decorator that measures the running time when a decorated function is called.
    """
    def wrapper(*args, **kwargs):
        start_time = get_current_datetime(as_str=False)
        result = func(*args, **kwargs)
        end_time = get_current_datetime(as_str=False)
        running_time = end_time - start_time
        total_seconds = running_time.total_seconds()
        minutes, seconds = divmod(total_seconds, 60)
        minutes = int(minutes)
        seconds = int(seconds)
        print_to_log(f'{func.__name__} run time: {minutes} minutes and {seconds} seconds.')
        print_to_log(f'#'*50)
        return result
    return wrapper

def get_current_mem_usage():
    """
    Returns the current memory usage of the program in gigabytes (GB).
    This function uses the psutil library to get the current memory usage of the program.
    This function assumes there is only a single process and no subprocesses in a job.
    This function needs to be updated if there are subprocesses.
    """
    # Get the current process ID
    pid = os.getpid()
    # Create a Process object using the process ID
    process = psutil.Process(pid)
    # Get the memory usage in bytes
    mem_usage = process.memory_info().rss
    # Convert bytes to gigabytes
    mem_usage_gb = mem_usage/(1024**3)
    mem_usage_gb = np.round(mem_usage_gb,2)
    return mem_usage_gb

def force_garbage_collection():
    """
    Force garbage collection to release the memory immediately.
    """
    mem_usage_before = get_current_mem_usage()
    gc.collect()
    mem_usage_after = get_current_mem_usage()
    mem_released = abs(mem_usage_before - mem_usage_after)
    mem_released = np.round(mem_released,2)
    print_to_log(f'Memory released from cleaning garbage: {mem_released} GB')
    print_to_log(f'Current memory usage: {mem_usage_after} GB')