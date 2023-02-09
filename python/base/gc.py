import gc

# refenrence: https://stackify.com/python-garbage-collection/


"""
    To see which Python you're using.
"""
def which_python():
    import platform
    print(platform.python_implementation)


"""
    Create a variablle and check its refrence count.
"""
def check_var_ref_count():
    import sys
    a = [0]
    count = sys.getrefcount(a)
    print(count)


"""
    check the configured thresholds of your garbage collector
"""
def check_config_thresholds():
    import gc
    print(gc.get_threshold())


"""
    check the number of objects in each of your generations
"""
def check_objects_nums():
    import gc
    print(gc.get_count())


"""
    trigger a manual garbage collection process
"""
def garbage_collect():
    import gc
    print(gc.get_count())
    print(gc.collect())
    print(gc.get_count())


garbage_collect()