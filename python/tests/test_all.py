import pytest
import mcc


def test_sum_as_string():
    assert mcc.sum_as_string(1, 1) == "2"
