import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401


def test_listoffset_inclusive():
    array = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    inclusive_jets_out = [
        [[[{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}]]],
        [[[{"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12}]]],
    ]
    assert inclusive_jets_out == cluster.inclusive_jets().to_list()
