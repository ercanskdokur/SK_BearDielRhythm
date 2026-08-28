#!/usr/bin/env python
# Multi-panel figure of ALL environmental predictor rasters (except precip/temp).
# Style mirrors the reference "Environmental predictors" panel. Python/geoenv, no container.
# Produces: Fig_EnvPredictors.png -> Figure S2 (SI)
import os, numpy as np
import rasterio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

B   = "/path/to/project"
ras = f"{B}/data/rasters"
OUT = f"{B}/BrownBearDielAct/figures/Fig_EnvPredictors.png"

# panel order: 8 distance rasters (rows 1-2) then 4 topographic (row 3)
# (name, file, cmap, cbar_label, to_km)
PANELS = [
    ("d2Builtareas",     "d2Builtareas.tif",     "viridis", "km", True),
    ("d2Crops",          "d2Crops.tif",          "viridis", "km", True),
    ("d2Forest",         "d2Forest.tif",         "viridis", "km", True),
    ("d2GarbageDump",    "d2GarbageDump.tif",    "viridis", "km", True),
    ("d2ProtectedAreas", "d2ProtectedAreas.tif", "viridis", "km", True),
    ("d2Rangeland",      "d2Rangeland.tif",      "viridis", "km", True),
    ("d2Roads",          "d2Roads.tif",          "viridis", "km", True),
    ("d2Water",          "d2Water.tif",          "viridis", "km", True),
    ("Elevation",        "Elevation.tif",        "terrain", "m a.s.l.", False),
    ("Slope",            "Slope.tif",            "magma",   "degrees", False),
    ("Roughness",        "Roughness.tif",        "cividis", "m", False),
    ("Aspect",           "Aspect.tif",           "twilight","degrees", False),
]

ncol, nrow = 4, 3
fig, axes = plt.subplots(nrow, ncol, figsize=(19, 13.5))
axes = axes.ravel()

crs_txt, res_txt = None, None
for ax, (name, fn, cmap, clab, to_km) in zip(axes, PANELS):
    with rasterio.open(f"{ras}/{fn}") as src:
        arr = src.read(1, masked=True)
        b = src.bounds
        if crs_txt is None:
            crs_txt = src.crs.to_string()
            res_txt = f"{abs(src.transform.a):.0f}"
    data = np.ma.masked_invalid(arr.astype(float))
    if to_km:
        data = data / 1000.0
    im = ax.imshow(data, cmap=cmap, extent=[b.left, b.right, b.bottom, b.top],
                   origin="upper", interpolation="nearest")
    ax.set_title(name, fontsize=15, fontweight="bold", pad=6)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.ax.tick_params(labelsize=10)
    cb.set_label(clab, fontsize=11, rotation=90, labelpad=2)

# hide any leftover axes (none here — 12 panels fill 12 slots)
for ax in axes[len(PANELS):]:
    ax.set_visible(False)

fig.suptitle("Environmental predictors", fontsize=22, fontweight="bold", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.975])
fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
print("CRS:", crs_txt, "res:", res_txt)
print("SAVED", OUT)
