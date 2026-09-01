import json
from pathlib import Path
from typing import Optional, Dict, List

import pandas as pd


def flatten_ontology(obj):
    """
    Recursively extract every ontology node that has an Allen 'id'.

    Handles standard Allen ontology JSON trees, including:
    - a root node dictionary with 'children'
    - a list of root nodes
    - wrapper dictionaries such as {'msg': [...]} or {'structures': [...]}
    """
    nodes = []

    def walk(x) -> None:
        if isinstance(x, dict):
            # An actual Allen structure / ontology node
            if "id" in x and ("name" in x or "acronym" in x):
                nodes.append(
                    {
                        "allen_id": str(x["id"]),
                        "allen_name": x.get("name"),
                        "allen_acronym": x.get("acronym"),
                        "parent_structure_id": (
                            str(x["parent_structure_id"])
                            if x.get("parent_structure_id") is not None
                            else None
                        ),
                        "structure_id_path": x.get("structure_id_path"),
                        "rgb_triplet": x.get("rgb_triplet"),
                    }
                )
            for value in x.values():
                if isinstance(value, (dict, list)):
                    walk(value)

        elif isinstance(x, list):
            for item in x:
                walk(item)

    walk(obj)
    unique = {}
    for node in nodes:
        unique[node["allen_id"]] = node
    return list(unique.values())

def load_allen_ontology(ontology_json: str | Path) -> pd.DataFrame:
    """Load Ontology.json and create an Allen-ID-to-name lookup table."""
    with open(ontology_json, "r", encoding="utf-8") as handle:
        ontology = json.load(handle)

    lookup = pd.DataFrame(flatten_ontology(ontology))

    if lookup.empty:
        raise ValueError(
            "No Allen structures were detected in Ontology.json. "
            "Inspect the JSON structure and adapt flatten_ontology()."
        )
    lookup["path_list"] = lookup["structure_id_path"].apply(parse_structure_id_path)
    parent_map = dict(
        zip(
            lookup["allen_id"].astype(str),
            lookup["parent_structure_id"].astype("string")
        )
    )
    parent_map = {
        str(k): (None if pd.isna(v) else str(v))
        for k, v in parent_map.items()
    }

    cache: Dict[str, List[str]] = {}

    missing = lookup["path_list"].map(len).eq(0)
    if missing.any():
        lookup.loc[missing, "path_list"] = lookup.loc[missing, "allen_id"].map(
            lambda x: compute_path_from_parents(str(x), parent_map, cache)
        )

    return lookup

def parse_structure_id_path(x):
    """
    Convert Allen structure_id_path to a list of IDs.
    Typical format: "/997/8/567/688/695/"
    """
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return []
    if isinstance(x, list):
        return [str(v) for v in x]
    x = str(x).strip()
    if not x:
        return []
    return [part for part in x.strip("/").split("/") if part]

def compute_path_from_parents(
    node_id: str,
    parent_map: Dict[str, Optional[str]],
    cache: Dict[str, List[str]],
) -> List[str]:
    """
    Compute root-to-node ancestry path from parent_structure_id.
    """
    node_id = str(node_id)

    if node_id in cache:
        return cache[node_id]

    seen = set()
    lineage = []
    cur = node_id

    while cur is not None:
        if cur in seen:
            raise ValueError(f"Cycle detected in ontology at node {cur}")
        seen.add(cur)
        lineage.append(cur)

        parent = parent_map.get(cur)
        if parent is None or parent == cur or parent not in parent_map:
            break

        cur = parent

    path = list(reversed(lineage))
    cache[node_id] = path
    return path

def build_lookup_dicts(ontology_df: pd.DataFrame):
    by_id = ontology_df.set_index("allen_id").to_dict(orient="index")
    return by_id

def restrict_region_to_depth(
    allen_id: Optional[str],
    ontology_by_id: Dict[str, Dict],
    target_depth: int,
    use_acronym: bool = False,
) -> Optional[str]:
    """
    Return the ancestor label at target_depth along structure_id_path.

    target_depth:
      0 = root
      1 = first level below root
      ...
      n = nth node in the Allen hierarchy path

    If the structure is shallower than target_depth, return the deepest available node.
    """
    if allen_id is None or pd.isna(allen_id):
        return None

    allen_id = str(allen_id)
    node = ontology_by_id.get(allen_id)
    if node is None:
        return None

    path_list = node.get("path_list", [])
    if not path_list:
        return node.get("allen_acronym" if use_acronym else "allen_name")

    idx = min(target_depth, len(path_list) - 1)
    ancestor_id = str(path_list[idx])

    ancestor = ontology_by_id.get(ancestor_id)
    if ancestor is None:
        return None

    return ancestor.get("allen_acronym" if use_acronym else "allen_name")

def add_depth_restricted_annotation(
    df: pd.DataFrame,
    ontology_df: pd.DataFrame,
    id_col: str = "atlas_region_id",
    target_depth: int = 3,
    use_acronym: bool = False,
) -> pd.DataFrame:
    ontology_by_id = build_lookup_dicts(ontology_df)

    out = df.copy()

    out["allen_name_full"] = out[id_col].map(
        lambda x: ontology_by_id.get(str(x), {}).get("allen_name")
        if pd.notna(x) else None
    )
    out["allen_acronym_full"] = out[id_col].map(
            lambda x: ontology_by_id.get(str(x), {}).get("allen_acronym")
            if pd.notna(x) else None
        )

    label_col = "allen_acronym_depth" if use_acronym else "allen_name_depth"
    out[label_col] = out[id_col].map(
            lambda x: restrict_region_to_depth(
                allen_id=x,
                ontology_by_id=ontology_by_id,
                target_depth=target_depth,
                use_acronym=False,
            )
        )

    return out