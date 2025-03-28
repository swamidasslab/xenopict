import pytest

from xenopict.monkey import Patcher


class _TestClass:
    clsattr = "class_value"

    def __init__(self):
        self.attr = "init_value"


def test_patcher_absent():
    P = Patcher()
    obj = _TestClass()

    P.replace(obj, "the_answer", 42)
    P.install()
    assert obj.the_answer == 42  # type: ignore
    assert "the_answer" in obj.__dict__

    P.uninstall()

    with pytest.raises(AttributeError):
        obj.the_answer  # type: ignore


def test_patcher_attr():
    P = Patcher()
    obj = _TestClass()

    P.replace(obj, "attr", 100)
    P.replace(obj, "new_attr", 42)

    with pytest.raises(AttributeError):
        obj.new_attr  # type: ignore

    with P:
        assert obj.new_attr == 42  # type: ignore
        assert obj.attr == 100

    assert obj.attr == "init_value"

    with pytest.raises(AttributeError):
        obj.new_attr  # type: ignore


def test_patcher_context():
    P = Patcher()
    obj = _TestClass()

    P.replace(obj, "attr", 42)
    with P:
        assert obj.attr == 42

    assert obj.attr == "init_value"


def test_patcher_clsattr():
    P = Patcher()
    obj = _TestClass()

    P.replace(_TestClass, "clsattr", 42)
    P.replace(_TestClass, "new_clsattr", 100)

    with pytest.raises(AttributeError):
        obj.new_clsattr  # type: ignore

    with P:
        assert obj.clsattr == 42
        assert obj.new_clsattr == 100  # type: ignore
        assert "new_clsattr" not in obj.__dict__
        assert hasattr(obj, "new_clsattr")

    P.uninstall()
    assert obj.clsattr == "class_value"
    assert "clsattr" not in obj.__dict__

    with pytest.raises(AttributeError):
        obj.new_clsattr  # type: ignore
