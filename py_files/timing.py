import functools
import pdb
import time

class Timing(object):
    time_dict = dict()
    count = 0
    def __init__(self, func):
        functools.update_wrapper(self, func)
        self.func = func
        #if n_dict != None:
        #    time_dict = n_dict
        self.__class__.count += 1
        print(self.__class__.count, self.func.__name__)

    def __call__(self, *args, **kwargs):
        pdb.set_trace()
        t0 = time.perf_counter()
        returned = self.func(*args, **kwargs)
        t1 = time.perf_counter()
        key_ = self.func.__name__
        _key = ','.join(map(lambda x: str(x), args))
        self.__class__.time_dict["{} ({})".format(key_, _key)] = t1-t0
        return returned
