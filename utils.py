import time






def func_timer(func):
    
    def wrap(*args, **kwargs):
        t0 = time.perf_counter()
        
        result = func(*args, **kwargs)
        
        tf = time.perf_counter()
        
        print(func.__name__, tf - t0)
        
        return result
        
    return wrap
        
