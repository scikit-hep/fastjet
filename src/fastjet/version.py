# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import re

__version__ = "3.4.0.2"
version = __version__
version_info = tuple(re.split(r"[-\.]", __version__))

del re
