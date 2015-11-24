from functools import wraps

from openalea.stat_tool._stat_tool import *
from openalea.stat_tool import __stat_tool

from types import ModuleType, ClassType

def wrapper(f):
    @wraps(f)
    def __str__(self):
        sstream = __stat_tool.std.Ostringstream(__stat_tool.std.ios_openmode.S__OUT)
        f(self, sstream)
        return sstream.str()
    return __str__

def get_lshift(obj):
    if isinstance(obj, ModuleType):
        for obj in obj.__dict__.itervalues():
            get_lshift(obj)
    elif hasattr(obj, '__lshift__'):
        try:
            obj.__str__ = wrapper(getattr(obj, '__lshift__'))
        except:
            pass

get_lshift(__stat_tool.stat_tool)
