"""Sets the version number. READTHEDOCS seems to struggle with this,
so whack in a try/except block to catch the error and set the version
number manually"""
try:
    from importlib.metadata import version as _version
    __version__ = _version("wodenpy")
except:
    __version__ = "2.6.0-alpha"