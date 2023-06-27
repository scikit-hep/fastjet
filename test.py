import awkward as ak
import fastjet
import vector

vector.register_awkward()

array1 = ak.Array(
    [
        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 23.5},
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 755.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 835.56},
    ],
)

array2 = ak.Array(
    [
        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 23.5},
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 755.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 835.56},
    ],
    with_name="Momentum4D",
)

array3 = ak.Array(
    [
        {"pt": 3.42, "eta": 1.24, "phi": 1.21, "M": 22.6},
        {"pt": 71.8, "eta": 2.72, "phi": 1.11, "M": 519},
        {"pt": 71.1, "eta": 2.73, "phi": 1.1, "M": 631},
    ],
    with_name="Momentum4D",
)

"""
array4 = ak.Array(
    [
        {"pt": 3.42, "eta": 1.24, "phi": 1.21, "M": 22.6},
        {"pt": 71.8, "eta": 2.72, "phi": 1.11, "M": 519},
        {"pt": 71.1, "eta": 2.73, "phi": 1.1, "M": 631},
    ],
)
"""

jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)

cluster1 = fastjet.ClusterSequence(array1, jetdef)
cluster2 = fastjet.ClusterSequence(array2, jetdef)
cluster3 = fastjet.ClusterSequence(array3, jetdef)
# cluster4 = fastjet.ClusterSequence(array4, jetdef)

inclusive_jets1 = cluster1.inclusive_jets().to_list()
inclusive_jets2 = cluster2.inclusive_jets().to_list()
inclusive_jets3 = cluster3.inclusive_jets().to_list()
# inclusive_jets4 = cluster4.inclusive_jets()

print(inclusive_jets1)
print(inclusive_jets2)
print(inclusive_jets3)

assert inclusive_jets1 == inclusive_jets2
assert inclusive_jets2 == inclusive_jets3
# assert inclusive_jets3 == inclusive_jets4
