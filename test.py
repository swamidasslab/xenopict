from xenopict.monkey import Patch, Patcher


def test_patcher_on_absent():
    P = Patcher()
    obj = []
    P.append(Patch(obj, "the_answer", 42, False))

    P.install()
    assert obj.the_answer == 42  # type: ignore
    assert "the_answer" in obj.__dict__

    P.uninstall()
    assert "the_answer" not in obj.__dict__
