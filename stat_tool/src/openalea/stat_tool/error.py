""" Error classes """
__revision__ = "$Id$"


from _stat_tool import _FormatError

__all__ = ["_FormatError",]

class StatToolError(Exception):
    def __init__(self, arg):
        Exception.__init__(self, str(arg))


