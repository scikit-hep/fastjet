import fastjet._ext  # noqa: F401, E402
import numpy as np
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class cluster:
    def __init__(self, data, R, algor):
        data = (
            data.astype(
                [
                    ("px", np.float32),
                    ("py", np.float32),
                    ("pz", np.float32),
                    ("E", np.float32),
                ]
            )
            .view(np.float32)
            .reshape(-1, 4)
        )
        data = data.tolist()
        self.out = fastjet._ext.interface(data, R, algor)
