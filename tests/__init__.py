from subprocess import run
import functools
import pytest


def shell(func):
    @functools.wraps(func)
    def wrapper(*kargs, **kwargs):
        cmd = func.__doc__
        run(cmd, shell=True)

    return wrapper


@pytest.fixture(scope="session")
@shell
def docker():
    """
    docker buildx build -t wf-rbp-motif:dev .
    """
