# Patched drop-in for vpt_plugin_cellpose/predict.py
#
# Purpose: stop a single degenerate (off-tissue / near-empty) tile from crashing
# the whole `vpt run-segmentation` job with:
#     OpenCV cv2.resize (-215:Assertion failed) inv_scale_x > 0
#
# Mechanism: the plugin already skips z-planes whose std < 0.1 and returns an
# empty mask when ALL z-planes are empty. The crash is the in-between case: a
# near-empty tile keeps 1-2 noisy z-planes, so cellpose runs do_3D on a volume
# that is too thin in z, and its 3D reslice resizes an axis to size 0.
#
# The only additions over the upstream file are marked "# PATCH". They fire ONLY
# on tiles that would otherwise crash (real tissue tiles have signal in every
# z-plane, so they are untouched and produce identical results). Cellpose version,
# model weights and all real-tile output are unchanged -> consistent with samples
# already segmented.
#
# Install: replace the file at
#   <vpt-env>/lib/python3.10/site-packages/vpt_plugin_cellpose/predict.py
# (in the traceback: /home/ubuntu/miniforge3/envs/vpt/lib/python3.10/site-packages/
#  vpt_plugin_cellpose/predict.py). Back up the original first. Dask workers import
# the module, so the patched file is picked up automatically; no code injection.

import warnings

import numpy as np
import cv2  # PATCH: used by the safety-net except below

from cellpose import models

from vpt_core.io.image import ImageSet
from vpt_plugin_cellpose import CellposeSegProperties, CellposeSegParameters


def run(images: ImageSet, properties: CellposeSegProperties, parameters: CellposeSegParameters) -> np.ndarray:
    warnings.filterwarnings("ignore", message=".*the `scipy.ndimage.filters` namespace is deprecated.*")
    is_valid_channels = parameters.nuclear_channel and parameters.entity_fill_channel
    image = (
        images.as_stack([parameters.nuclear_channel, parameters.entity_fill_channel])
        if is_valid_channels
        else images.as_stack()
    )
    empty_z_levels = set()
    for z_i, z_plane in enumerate(image):
        for channel_i in range(z_plane.shape[-1]):
            if z_plane[..., channel_i].std() < 0.1:
                empty_z_levels.add(z_i)

    empty_mask = np.zeros((image.shape[0],) + image.shape[1:-1])  # PATCH: shared empty-mask helper

    if len(empty_z_levels) == image.shape[0]:
        return empty_mask

    to_segment_z = list(set(range(image.shape[0])).difference(empty_z_levels))

    # PATCH: with only 1 non-empty z-plane, cellpose's 3D reslice rescales the
    # z-axis (~0.31x) to size 0 -> cv2.resize "inv_scale_x > 0". Such a tile holds
    # no real 3D nucleus, so return an empty mask, like the all-empty branch above.
    # Cut at < 2 (only the certain nz=1 crash); the 2-plane borderline case is left
    # to run and is covered by the except below if it ever also rounds to 0. Real
    # tiles have signal in every z-plane and never hit this.
    if properties.model_dimensions == "3D" and len(to_segment_z) < 2:
        return empty_mask

    if properties.custom_weights:
        model = models.CellposeModel(gpu=False, pretrained_model=properties.custom_weights, net_avg=False)
    else:
        model = models.Cellpose(gpu=False, model_type=properties.model, net_avg=False)

    try:  # PATCH: belt-and-suspenders for any remaining degenerate-resize crash
        mask = model.eval(
            image[to_segment_z, ...],
            z_axis=0,
            channel_axis=len(image.shape) - 1,
            diameter=parameters.diameter,
            flow_threshold=parameters.flow_threshold,
            mask_threshold=parameters.mask_threshold,
            resample=False,
            min_size=parameters.minimum_mask_size,
            tile=True,
            do_3D=(properties.model_dimensions == "3D"),
        )[0]
    except cv2.error:  # PATCH: degenerate tile -> no cells, same as empty branch
        return empty_mask

    mask = mask.reshape((len(to_segment_z),) + image.shape[1:-1])
    for i in empty_z_levels:
        mask = np.insert(mask, i, np.zeros(image.shape[1:-1]), axis=0)
    return mask
