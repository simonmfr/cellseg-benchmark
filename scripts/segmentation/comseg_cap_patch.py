"""Cap Louvain sweeps in ComSeg 1.8.5, which can loop forever on single patches (sopa issue #292)."""

import comseg.utils.custom_louvain as cl

with open(cl.__file__) as f:
    src = f.read()

if "len(list_move) < 100" not in src:
    src = src.replace("import random\n", "import random\nimport sys\n", 1)
    src = src.replace(
        "    while nb_moves > 0:", "    while nb_moves > 0 and len(list_move) < 100:"
    )
    src = src.replace(
        "    partition = list(filter(len, partition))\n    inner_partition",
        "    if len(list_move) >= 100:\n"
        "        print('COMSEG_CAP_HIT', list_move[:12], file=sys.stderr, flush=True)\n"
        "    partition = list(filter(len, partition))\n    inner_partition",
    )
    assert "import sys" in src and "len(list_move) < 100" in src, (
        "ComSeg cap patch did not apply"
    )
    with open(cl.__file__, "w") as f:
        f.write(src)
print("comseg cap patch applied")
