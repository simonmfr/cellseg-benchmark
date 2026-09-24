# Checklist for novel methods added to the master sdata:

1. Has the method run successfully?
2. Does it contain the location of the cells in sdata['table'].obsm['spatial']?
3. Does cell annotation exist?
4. **For 2D methods:** Is Ovrlpy computed?
5. **For 3D methods:** Are channel-intensities computed for each plane and averaged per cell?
6. Determine the scaling of the segmented boundaries -> micron or pixel space?
7. Double check the transformation setting within `cellseg_benchmark.sdata_utils.assign_transformations`
8. **For 3D methods:** Do you need to separately load 3D boundaries? -> check implementation in `cellseg_benchmark.sdata_utils.build_shapes`
9. Is `adata.uns.spatialdata_attrs` correctly set? Does the boundary data have a matching index?
