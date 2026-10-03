"""Cap Louvain sweeps in ComSeg 1.8.5, which can loop forever on single patches (sopa issue #292)."""

import comseg.utils.custom_louvain as cl

p = cl.__file__
s = open(p).read()
s = s.replace("    while nb_moves > 0:", "    while nb_moves > 0 and len(list_move) < 100:")
s = s.replace(
    "    partition = list(filter(len, partition))\n    inner_partition",
    "    if len(list_move) >= 100: print('COMSEG_CAP_HIT', list_move[:12], file=__import__('sys').__stderr__, flush=True)\n"
    "    partition = list(filter(len, partition))\n    inner_partition",
)
assert "len(list_move) < 100" in s, "ComSeg cap patch did not apply"
open(p, "w").write(s)
print("comseg cap patch applied: True")
