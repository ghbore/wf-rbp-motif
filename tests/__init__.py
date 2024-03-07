from subprocess import run
import functools


def shell(func):
    @functools.wraps(func)
    def wrapper(*kargs, **kwargs):
        cmd = func.__doc__
        run(cmd, shell=True)

    return wrapper
