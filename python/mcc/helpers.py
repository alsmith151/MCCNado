import numpy as np

def define_memory_requested(attempts: int = 1, initial_value: int  = 1, scale: float = 1) -> str:
    """
    Define the memory requested for the job.
    """
    memory = int(initial_value) * 2 ** (int(attempts) - 1)
    memory = memory * float(scale)
    return f"{memory}G"

def define_time_requested(attempts: int = 1, initial_value: int = 1, scale: float = 1) -> str:
    """
    Define the time requested for the job.

    Base time is 1 hour.
    """
    time = int(initial_value) * 2 ** (int(attempts) - 1)
    time = time * float(scale)
    return f"{time}h"

def is_on(param: str) -> bool:
    """
    Returns True if parameter in "on" values
    On values:
        - true
        - t
        - on
        - yes
        - y
        - 1
    """
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_off(param: str):
    """Returns True if parameter in "off" values"""
    values = ["", "none", "f", "n", "no", "false", "0"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_none(param: str) -> bool:
    """Returns True if parameter is none"""
    values = ["", "none"]
    if str(param).lower() in values:
        return True
    else:
        return False


def check_options(value: object):
    if value in [None, np.nan, ""]:
        return ""
    elif is_off(value):
        return ""
    else:
        return value