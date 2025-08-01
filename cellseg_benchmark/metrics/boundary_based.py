import geopandas as gpd
from spatialdata.models import PointsModel


def count_assigned_transcripts(
    sdata, sdata_transcripts_key, boundaries_to_process=None
):
    """Calculate transcript assignment to cells across boundaries, handling 3D z-planes."""
    # Load transcripts
    transcripts_df = PointsModel.parse(
        sdata[sdata_transcripts_key], sort=True
    ).compute()[["x", "y", "gene", "transcript_id", "global_z"]]
    transcripts_gdf = gpd.GeoDataFrame(
        transcripts_df, geometry=gpd.points_from_xy(transcripts_df.x, transcripts_df.y)
    )
    print(f"Loaded {len(transcripts_gdf)} transcripts")

    # Filter boundaries
    available_boundaries = {
        k: v for k, v in sdata.shapes.items() if k.startswith("boundaries_")
    }
    if boundaries_to_process:
        boundary_keys = [
            name if name.startswith("boundaries_") else f"boundaries_{name}"
            for name in boundaries_to_process
        ]
        available_boundaries = {
            k: v for k, v in available_boundaries.items() if k in boundary_keys
        }

    results = {}

    for boundaries_name, boundaries in available_boundaries.items():
        dataset_name = boundaries_name.replace("boundaries_", "")
        print(f"Processing {dataset_name}")

        if (
            boundaries is None
            or boundaries.empty
            or "geometry" not in boundaries.columns
        ):
            print("  Skipping - invalid geometry")
            continue

        try:
            if "ZIndex" not in boundaries.columns:
                # 2D case
                assigned, total = process_boundaries(boundaries, transcripts_gdf)
                dimension = "2d"
            else:
                unique_z = boundaries["ZIndex"].unique()
                if len(unique_z) == 1:
                    # Single z-plane
                    z_boundaries = boundaries[boundaries["ZIndex"] == unique_z[0]]
                    z_transcripts = transcripts_gdf[
                        transcripts_gdf["global_z"] == unique_z[0]
                    ]
                    assigned, total = _process_boundaries(z_boundaries, z_transcripts)
                    dimension = "2d"
                else:
                    # Multiple z-planes
                    print(f"  Processing {len(unique_z)} z-planes")
                    total_assigned, total_transcripts = 0, 0

                    for z_index in sorted(unique_z):
                        z_boundaries = boundaries[boundaries["ZIndex"] == z_index]
                        z_transcripts = transcripts_gdf[
                            transcripts_gdf["global_z"] == z_index
                        ]

                        z_assigned, z_total = _process_boundaries(
                            z_boundaries, z_transcripts
                        )
                        total_assigned += z_assigned
                        total_transcripts += z_total

                    assigned, total = total_assigned, total_transcripts
                    dimension = "3d"

            results[dataset_name] = {
                "assigned_count": assigned,
                "unassigned_count": total - assigned,
                "total_count": total,
                "pct_assigned": round((assigned / total) * 100, 2) if total > 0 else 0,
                "dimension": dimension,
            }
            print(
                f"  {assigned}/{total} assigned ({results[dataset_name]['pct_assigned']}%) - {dimension}"
            )

        except Exception as e:
            print(f"  Error: {e}")

    return results


def _process_boundaries(boundaries_gdf, transcripts_subset):
    """Helper to process boundaries and return assigned count."""
    if boundaries_gdf.empty or transcripts_subset.empty:
        return 0, len(transcripts_subset)

    boundaries_unified = gpd.GeoDataFrame(
        geometry=[boundaries_gdf["geometry"].union_all(method="unary")]
    )
    joined = gpd.sjoin(
        transcripts_subset, boundaries_unified, how="left", predicate="within"
    )
    assigned = joined["index_right"].notna().sum()
    return assigned, len(transcripts_subset)
