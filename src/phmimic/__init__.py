"""ph-mimic: EBNA1 molecular-mimicry candidate screen.

The pipeline is deliberately built around one experimentally anchored seed
(the EBNA1 386-440 cross-reactive region) and always reports what a *random
peptide of the same kind* would have found. Any number that cannot be
recomputed from a delivered table is not considered a result.
"""

__version__ = "1.0.0"
