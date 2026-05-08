# sdata.zarr requirements

* boundaries have an unnamed index and a separate column for cell identifier called `cell_id`. If 3D boundaries are given, consider providing the levels converted to float
* the table has `spatial_attrs` setup in `.uns`. Cell identifier are also called `cell_id`.