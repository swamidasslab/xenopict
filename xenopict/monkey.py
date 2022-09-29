import sys
import inspect
from decorator import decorator
from typing import Callable, Optional


class Patch:
    _installed = False

    def __enter__(self):
        self.install()

    def __exit__(self, *args, **kwargs):
        self.uninstall()

    def install(self):
        pass

    def uninstall(self):
        pass


class Patcher(Patch):
    def __init__(self, patches: Optional[list[Patch]] = None):
        self.patches: list[Patch] = [] if patches is None else patches

    def append(self, patch: Patch) -> None:
        self.patches.append(patch)

    def replace(self, obj, attr, new):
        self.append(Replace(obj, attr, new))

    def install(self):
        if not self._installed:
            for p in self.patches:
                p.install()
        self._installed = True

    def uninstall(self):
        if self._installed:
            for p in reversed(self.patches):
                p.uninstall()
        self._installed = False


class Replace(Patch):
    _ABSENT = []

    def __init__(self, obj, attr, replace_with):
        self.obj = obj
        self.attr = attr
        self.replace_with = replace_with
        self.absent = None

    def uninstall(self):
        if not self._installed:
            return

        if self.orig is self._ABSENT:
            if self.attr in self.obj.__dict__:
                delattr(self.obj, self.attr)
        else:
            setattr(self.obj, self.attr, self.orig)

        del self.orig
        self._installed = False

    def __repr__(self):

        name = self.obj.__name__

        if hasattr(self.obj, "__module__"):
            name = f"{self.obj.__module__}.{name}"

        return f"{self.__class__.__name__}({name}.{self.attr})"

    def wrap(self, replace_with):
        return replace_with

    def install(self):
        if self._installed:
            return

        self.orig = self.obj.__dict__.get(self.attr, self._ABSENT)

        wrap = self.wrap(self.replace_with)
        setattr(self.obj, self.attr, wrap)
        self._installed = True


class Wrap(Replace):
    def wrap(self, replace_with):
        return replace_with(self.orig)


class BoostError(RuntimeError):
    pass


class BoostTypeError(TypeError):
    pass


class BoostModulePatcher(Patcher):
    def __init__(self, mod=None, wrappers=None, depth=15):
        super().__init__()
        self.wrappers = [] if wrappers is None else wrappers
        self._crawled_objects: dict = {}

        self._before_call: list[
            Callable[[tuple, dict, Callable], tuple[tuple, dict]]
        ] = []
        self._after_call = []
        self._on_throw = []

        if mod:
            self._crawl(mod, mod.__name__, depth)

    def _crawl(self, mod, base, depth=3):
        self._crawled_objects[id(mod)] = None

        if depth < 1:
            return

        for name in dir(mod):
            if name.startswith("__"):
                continue

            obj = getattr(mod, name)

            if id(obj) in self._crawled_objects:
                if wrapped := self._crawled_objects[id(obj)]:
                    self.replace(mod, obj.__name__, wrapped)
                continue

            if inspect.ismodule(obj):
                if obj.__name__.startswith(base):
                    self._crawl(obj, base, depth - 1)
                continue

            if inspect.isclass(obj):
                continue
            # if obj.__module__.startswith(base):
            #     self._crawl(obj, base, depth - 1)

            if inspect.isbuiltin(obj):
                wrapped = self.boost_wrapper(obj)
                for w in self.wrappers:
                    wrapped = w(wrapped)

                self._crawled_objects[id(obj)] = wrapped
                self.replace(mod, obj.__name__, wrapped)
        return

    def _handle_boost_error(self, func, args, kwargs, exec_info):
        exe_info = self.on_throw(func, args, kwargs, exec_info)
        exc, info, tb = exec_info

        boost_error_name = exc.__name__ if exc else "Unknown"
        if boost_error_name == "ArgumentError":
            return BoostTypeError(info).with_traceback(tb)

        return BoostError(f"{boost_error_name}: {info}").with_traceback(tb)

    def boost_wrapper(self, func):
        def wrap(*args, **kwargs):
            try:
                args, kwargs = self.before_call(args, kwargs, func)
                result = func(*args, **kwargs)

                result = self.after_call(result, args, kwargs, func)
            except Exception:
                raise self._handle_boost_error(
                    func, args, kwargs, sys.exc_info()
                ) from None

        wrap.__name__ = func.__name__
        wrap.__qualname__ = func.__name__
        wrap.__module__ = func.__module__
        wrap.__doc__ = func.__doc__
        #    wrap.__wrapped__ = func  # breaks wrapping on boost functions

        return wrap

    def register_before_call(
        self, callback: Callable[[tuple, dict, Callable], tuple[tuple, dict]]
    ):
        self._before_call.append(callback)

    def before_call(
        self, args: tuple, kwargs: dict, func: Callable
    ) -> tuple[tuple, dict]:
        for cb in self._before_call:
            args, kwargs = cb(args, kwargs, func)
        return args, kwargs

    def after_call(self, result, args: tuple, kwargs: dict, func: Callable):
        return result

    def on_throw(self, func: Callable, args: tuple, kwargs: dict, exec_info):
        return exec_info
