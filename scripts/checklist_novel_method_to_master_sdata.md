# Checklist for novel methods added to the master sdata:

1. Has the method run successfully?
2. Does cell annotation exist?
3. **For 2D methods:** Is Ovrlpy computed?
4. **For 3D methods:** Are channel-intensities computed for each plane and averaged per cell?
5. Determine the scaling of the segmented boundaries -> micron or pixel space?
6. Double check the transformation setting within `cellseg_benchmark.sdata_utils.assign_transformations`
7. **For 3D methods:** Do you need to separately load 3D boundaries? -> check implementation in `cellseg_benchmark.sdata_utils.build_shapes`
8. Is `adata.uns.spatialdata_attrs` correctly set? Does the boundary data have a matching index?
