BASE_PATH = "/dss/dssfs03/pn52re/pn52re-dss-0001/cellseg-benchmark"

# WT samples from the aging cohort used as additional controls for the htra1 cohort
htra1_aging_controls = [
    "aging_s1_r1",
    "aging_s5_r1",
    "aging_s6_r0",
    "aging_s7_r1",
    "aging_s7_r2",
    "aging_s8_r2",
    "aging_s11_r0",
]

# matched as prefixes against method names
methods_3D = [
    "Baysor_3D",
    "Proseg_3D",
    "SIS",
    "Watershed_Merlin",
    "vpt_3D",
]

image_based = ["Cellpose", "Negative_Control"]

method_colors = {
    # Baysor variants (red palette)
    "Baysor_2D_denovo": "#5c0010",
    "Baysor_3D_denovo": "#7a0012",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.2": "#9a0013",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.8": "#a71423",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.2": "#b32833",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.8": "#c03b43",
    "Baysor_2D_Cellpose_1_nuclei_model_1.0": "#cc4f53",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.2": "#d96363",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.8": "#e67673",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.2": "#f28a83",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.8": "#ff9e93",
    # vpt variants (blue palette)
    "vpt_2D_DAPI_PolyT": "#08306b",
    "vpt_2D_DAPI_PolyT_nuclei": "#08519c",
    "vpt_2D_DAPI_nuclei": "#2171b5",
    "vpt_3D_DAPI_PolyT": "#6baed6",
    "vpt_3D_DAPI_PolyT_nuclei": "#9ecae1",
    "vpt_3D_DAPI_nuclei": "#a9d5f1",
    # Cellpose core methods (green palette)
    "Cellpose_1_nuclei_model": "#00441b",
    "Cellpose_1_DAPI_PolyT": "#006d2c",
    "Cellpose_1_DAPI_Transcripts": "#217a37",
    "Cellpose_1_Merlin": "#31a354",
    "Cellpose_2_DAPI_PolyT": "#74c476",
    "Cellpose_2_DAPI_Transcripts": "#9cd8a2",
    # Proseg variants (outdated)
#    "Proseg_pure": "#4a1486",
#    "Proseg_Cellpose_1_DAPI_PolyT": "#5a0876",
#    "Proseg_Cellpose_1_DAPI_Transcripts": "#6a1894",
#    "Proseg_Cellpose_1_nuclei_model": "#7d3db3",
#    "Proseg_Cellpose_2_DAPI_PolyT": "#927ac6",
#    "Proseg_Cellpose_2_DAPI_Transcripts": "#b2a4db",
    # Proseg 3D variants (yellow palette)
    "Proseg_3D_Cellpose_1_DAPI_PolyT": "#f5f0df",
    "Proseg_3D_Cellpose_1_DAPI_Transcripts": "#f3e4bf",
    "Proseg_3D_Cellpose_1_nuclei_model": "#f1d9a1",
    "Proseg_3D_Cellpose_2_DAPI_PolyT": "#efce82",
    "Proseg_3D_Cellpose_2_DAPI_Transcripts": "#eec364",
    "Proseg_3D_vpt3D_DAPI_nuclei": "#edb846",
    "Proseg_3D_vpt3D_DAPI_PolyT": "#eaad28",
    "Proseg_3D_vpt3D_DAPI_PolyT_nuclei": "#e8a20a",
    # Negative controls (grey palette)
    "Negative_Control_Rastered_5": "#101010",
    "Negative_Control_Rastered_10": "#424141",
    "Negative_Control_Rastered_25": "#5A5A5A",
    "Negative_Control_Voronoi": "#797878",
    "Negative_Control_Visium": "#969696",
    # ComSeg standalone (dark purple)
    "ComSeg": "#d7f035",
    "SIS_DAPI_total_mrna": "#9cb01c",
    "Watershed_Merlin": "#8a9159",
    # Ficture based
    'Ficture_segments_dapi': "#5a0876",
    'Ficture_segments': "#7d3db3"
}

# Figure labels: <Algorithm> [3D] <input> [framework] [(prior, confidence)].
# P = PolyT, T = transcripts, N = nuclei; all image-based methods also use DAPI.
# Framework (Sopa, VPT, Merlin) only where a method ran inside one.
method_names = {
    "Baysor_2D_denovo": "Baysor",
    "Baysor_3D_denovo": "Baysor 3D",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.2": "Baysor Sopa (CP1 P, 0.2)",
    "Baysor_2D_Cellpose_1_DAPI_PolyT_0.8": "Baysor Sopa (CP1 P, 0.8)",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.2": "Baysor Sopa (CP1 T, 0.2)",
    "Baysor_2D_Cellpose_1_DAPI_Transcripts_0.8": "Baysor Sopa (CP1 T, 0.8)",
    "Baysor_2D_Cellpose_1_nuclei_model_1.0": "Baysor Sopa (CP1 N, 1.0)",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.2": "Baysor Sopa (CP2 P, 0.2)",
    "Baysor_2D_Cellpose_2_DAPI_PolyT_0.8": "Baysor Sopa (CP2 P, 0.8)",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.2": "Baysor Sopa (CP2 T, 0.2)",
    "Baysor_2D_Cellpose_2_DAPI_Transcripts_0.8": "Baysor Sopa (CP2 T, 0.8)",
    "vpt_2D_DAPI_PolyT": "Cellpose1 P VPT",
    "vpt_2D_DAPI_PolyT_nuclei": "Cellpose1 P+N VPT",
    "vpt_2D_DAPI_nuclei": "Cellpose1 N VPT",
    "vpt_3D_DAPI_PolyT": "Cellpose1 3D P VPT",
    "vpt_3D_DAPI_PolyT_nuclei": "Cellpose1 3D P+N VPT",
    "vpt_3D_DAPI_nuclei": "Cellpose1 3D N VPT",
    "Cellpose_1_nuclei_model": "Cellpose1 N Sopa",
    "Cellpose_1_DAPI_PolyT": "Cellpose1 P Sopa",
    "Cellpose_1_DAPI_Transcripts": "Cellpose1 T Sopa",
    "Cellpose_1_Merlin": "Cellpose1 P Merlin",
    "Cellpose_2_DAPI_PolyT": "Cellpose2 P Sopa",
    "Cellpose_2_DAPI_Transcripts": "Cellpose2 T Sopa",
    "Proseg_3D_Cellpose_1_DAPI_PolyT": "Proseg 3D (CP1 P)",
    "Proseg_3D_Cellpose_1_DAPI_Transcripts": "Proseg 3D (CP1 T)",
    "Proseg_3D_Cellpose_1_nuclei_model": "Proseg 3D (CP1 N)",
    "Proseg_3D_Cellpose_2_DAPI_PolyT": "Proseg 3D (CP2 P)",
    "Proseg_3D_Cellpose_2_DAPI_Transcripts": "Proseg 3D (CP2 T)",
    "Proseg_3D_vpt3D_DAPI_nuclei": "Proseg 3D (CP1 3D N)",
    "Proseg_3D_vpt3D_DAPI_PolyT": "Proseg 3D (CP1 3D P)",
    "Proseg_3D_vpt3D_DAPI_PolyT_nuclei": "Proseg 3D (CP1 3D P+N)",
    "Negative_Control_Rastered_5": "Raster 5µm",
    "Negative_Control_Rastered_10": "Raster 10µm",
    "Negative_Control_Rastered_25": "Raster 25µm",
    "Negative_Control_Voronoi": "Voronoi",
    "Negative_Control_Visium": "Visium",
    "ComSeg": "ComSeg",
    "SIS_DAPI_total_mrna": "Spots-In-Space 3D",
    "Watershed_Merlin": "Watershed 3D P Merlin",
    "Ficture_segments_dapi": "Ficture Seg (CP1 3D N)",
    "Ficture_segments": "Ficture Seg",
}

if set(method_names) != set(method_colors):
    raise ValueError("method_names and method_colors must cover the same methods")

if len(set(method_names.values())) != len(method_names):
    raise ValueError("duplicate method labels")

label_colors = {method_names[k]: c for k, c in method_colors.items()}

label_order = list(label_colors)

ficture_factor_to_celltype = {
    "0": "ABCs",
    "1": "Astrocytes",
    "2": "BAMs",
    "3": "Bergmann",
    "4": "ECs",
    "5": "Ependymal",
    "6": "Immune-Other",
    "7": "Microglia",
    "8": "Neurons-Dopa",
    "9": "Neurons-Dopa-Gaba",  # renamed to Neurons-Dopa, see true_cluster dict
    "10": "Neurons-Gaba",
    "11": "Neurons-Glut",
    "12": "Neurons-Glyc-Gaba",
    "13": "Neurons-Immature",  # renamed to Neurons-Granule-Immature, see true_cluster dict
    "14": "Neurons-Other",
    "15": "OECs",
    "16": "OPCs",
    "17": "Oligodendrocytes",
    "18": "Pericytes",
    "19": "SMCs",
    "20": "VLMCs",
}

true_cluster = {  # fine label -> canonical cell type
    "Neurons-Immature": "Neurons-Granule-Immature",  # legacy FICTURE factor label
    "Astrocytes": "Astrocytes",
    "Astroependymal": "Astrocytes",
    "BAMs": "BAMs",
    "Choroid-Plexus": "Ependymal",
    "ECs": "ECs",
    "Ependymal": "Ependymal",
    "Immune-Other": "Immune-Other",
    "Microglia": "Microglia",
    "Neurons-Dopa-Gaba": "Neurons-Dopa",
    "Neurons-Gaba": "Neurons-Gaba",
    "Neurons-Glut": "Neurons-Glut",
    "Neurons-Glyc-Gaba": "Neurons-Glyc-Gaba",
    "Neurons-Granule-Immature": "Neurons-Granule-Immature",
    "Neurons-Other": "Neurons-Other",
    "OECs": "OECs",
    "OPCs": "OPCs",
    "Oligodendrocytes": "Oligodendrocytes",
    "Pericytes": "Pericytes",
    "SMCs": "SMCs",
    "VLMCs": "VLMCs",
    "ABCs": "VLMCs",
    "Bergmann": "Astrocytes",
    "Neurons-Dopa": "Neurons-Dopa",
    "Tanycytes": "Ependymal",
}

# merges neurons for marker-gene-based metrics
merged_celltypes = {
    "Neurons-Dopa": "Neurons",
    "Neurons-Dopa-Gaba": "Neurons",
    "Neurons-Gaba": "Neurons",
    "Neurons-Glut": "Neurons",
    "Neurons-Glyc-Gaba": "Neurons",
}

# Vascular cell types, used for the secondary macro-F1 (see metrics.ficture).
vascular_celltypes = ["ECs", "Pericytes", "SMCs", "VLMCs"]

# Cell types without clear marker genes: their boundary annotation is unreliable, so
# transcripts touching them are dropped from the FICTURE F1 comparison (metrics.ficture).
unreliable_celltypes = ["Neurons-Other", "Immune-Other", "OECs"]

index_order = [
    "Astrocytes",
    "BAMs",
    "ECs",
    "Ependymal",
    "Immune-Other",
    "Microglia",
    "Neurons-Dopa",
    "Neurons-Gaba",
    "Neurons-Glut",
    "Neurons-Glyc-Gaba",
    "Neurons-Other",
    "OPCs",
    "Oligodendrocytes",
    "Pericytes",
    "SMCs",
    "VLMCs",
    "Astroependymal",
    "Choroid-Plexus",
    "Neurons-Granule-Immature",
    "Tanycytes",
    "Undefined",
    "Unknown",
    "Low-Read-Cells",
]

column_order = [
    "Astrocytes",
    "BAMs",
    "ECs",
    "Ependymal",
    "Immune-Other",
    "Microglia",
    "Neurons-Dopa",
    "Neurons-Gaba",
    "Neurons-Glut",
    "Neurons-Glyc-Gaba",
    "Neurons-Other",
    "OPCs",
    "Oligodendrocytes",
    "Pericytes",
    "SMCs",
    "VLMCs",
    "ABCs",
    "Bergmann",
    "Neurons-Dopa-Gaba",
    "Neurons-Immature",
    "OECs",
]

cell_type_colors = {
    "ECs": "#FF6464",
    "aECs": "#FF7700",  # for EC subtyping script
    "capECs": "#FF6464",  # for EC subtyping script
    "vECs": "#9966CC",  # for EC subtyping script
    "otherECs": "#FFBFBF",  # for EC subtyping script
    "Pericytes": "#F6EC2A",
    "SMCs": "#29FBA7",
    "VLMCs": "#85B0F9",
    "ABCs": "#AEC9F5",
    "Ependymal": "#FDC000",
    "Tanycytes": "#FFE180",
    "Choroid-Plexus": "#BF9800",
    "Astrocytes": "#FE9A30",
    "Astroependymal": "#FE9A30",
    "Bergmann": "#FFD1A3",
    "Oligodendrocytes": "#4564FF",
    "OPCs": "#0095FF",
    "Microglia": "#00C088",
    "BAMs": "#20B2AA",
    "Immune-Other": "#98DF8A",
    "Neurons-Gaba": "#B449F8",
    "Neurons-Glut": "#CEB3FF",
    "Neurons-Glyc-Gaba": "#DCAEFF",
    "Neurons-Dopa": "#FCA0FF",
    # "Neurons-Immature": "#FF50E5",
    "Neurons-Granule-Immature": "#FF50E5",
    "Neurons-Other": "#FCA0FF",
    "Neurons": "#B449F8",  # for marker f1 score where Neurons are merged in one celltype
    "OECs": "#9EDAE5",
    "Unknown": "#D9D9D9",  # = not found in mmc dict, see process_mapmycells_output()
    "Undefined": "#D9D9D9",  # = below QC threshold
    "Mixed": "#D9D9D9",
    "Low-Read-Cells": "#D9D9D9",
}  # Updated method-to-color mapping with distinct, moderately saturated shades

brain_regions_colors = {
    "BS": "FF7080",
    "CA3sp": "66A83D",
    "CTX": "B0FFB8",
    "DG-sg": "66A83D",
    "HIP": "7ED04B",
    "STR": "98D6F9",
    "VS": "AAAAAA",
    "fiber tracts": "CCCCCC",
    "BS/STR": "CCA3BC",
    "STR/CTX": "A4EAD8",
    "Meninges": "480091",
}
