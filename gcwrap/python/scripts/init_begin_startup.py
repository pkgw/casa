## import original sys.path to insure that
## it is initialized before any adjustments
from casa_system import original_path
import sys

if len(sys.path[0]) == 0:
    sys.path = sys.path[1:]

