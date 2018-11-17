import sys
if sys.version_info[0] < 3:
    from NAFFlib2_c import *
else:
    from .NAFFlib_c import *

