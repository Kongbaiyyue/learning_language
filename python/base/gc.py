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


"""
    search object referrers
"""
def object_referrers():
    import gc


    def foo():
        a = [2, 4, 6]
        b = [1, 4, 7]

        l = [a, b]
        d = dict(a=a)
        return l, d

    l, d = foo()
    r1 = gc.get_referrers(l[0])
    r2 = gc.get_referrers(l[1])

    print(r1)
    print(r2)


object_referrers()

