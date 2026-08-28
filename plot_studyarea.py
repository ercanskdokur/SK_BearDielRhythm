#!/usr/bin/env python
# Study-area land-cover map (ESRI 2022 LULC) — Python/matplotlib, no container.
# Produces: Fig_StudyArea_New.png -> Figure 2 (main text)
import os, math, numpy as np, pandas as pd
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.enums import Resampling
from rasterio.transform import from_origin
import geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

B   = "/path/to/project"
shp = f"{B}/data/shapefiles"; esri = f"{B}/data/esri_lulc"; wdpa = f"{B}/data/wdpa"
OUT = f"{B}/BrownBearDielAct/figures/Fig_StudyArea_New.png"
DST = "EPSG:3857"

# ESRI class -> target palette (sampled from Fig_studyarea_map.png legend)
COL = {1:"#3c8dde", 2:"#1b7837", 4:"#6fae8f", 5:"#e69f00",
       7:"#b0322b", 8:"#b0aba3", 9:"#ffffff", 11:"#e5dca3"}
LAB = {2:"Forest", 5:"Crops", 7:"Built-up areas", 8:"Bare ground",
       11:"Rangeland", 1:"Water"}

def rd(p): return gpd.read_file(p).to_crs(DST)

# ---- frame from study area ----
sa = rd(f"{shp}/StudyAreaSk.shp")
xmin, ymin, xmax, ymax = sa.total_bounds
px, py = (xmax-xmin)*0.012, (ymax-ymin)*0.012
xmin-=px; xmax+=px; ymin-=py; ymax+=py
res = 90.0
Wd = int((xmax-xmin)/res); Hd = int((ymax-ymin)/res)
dst_tr = from_origin(xmin, ymax, res, res)
print("frame 3857:", round(xmin), round(ymin), round(xmax), round(ymax), "grid", Wd, Hd)

# ---- ESRI mosaic reprojected to target grid (cached) ----
cache = f"{B}/BrownBearDielAct/figures/.mosaic_{Wd}x{Hd}.npy"
if os.path.exists(cache):
    mosaic = np.load(cache); print("mosaic loaded from cache", mosaic.shape)
else:
    tiles = [f"{esri}/{t}_20220101-20230101.tif" for t in ("37S","37T","38S","38T")]
    mosaic = np.zeros((Hd, Wd), dtype=np.uint8)
    for t in tiles:
        with rasterio.open(t) as src:
            with WarpedVRT(src, crs=DST, transform=dst_tr, width=Wd, height=Hd,
                           resampling=Resampling.nearest) as vrt:
                a = vrt.read(1)
        mosaic = np.where(mosaic == 0, a, mosaic)
        print("merged", os.path.basename(t), "nonzero%%=%.1f" % (100*(mosaic>0).mean()))
    np.save(cache, mosaic)

# ---- RGBA raster image ----
rgba = np.zeros((Hd, Wd, 4), dtype=float)
for v, hx in COL.items():
    rgba[mosaic == v] = to_rgba(hx)

# ---- vectors ----
roads = rd(f"{shp}/MainRoads.shp")
amnp  = rd(f"{shp}/SK_AMNP.shp")
dump  = rd(f"{shp}/GarbageDumpSites.shp")
arm   = rd(f"{wdpa}/ARM/WDPA_WDOECM_Aug2026_Public_ARM_shp-polygons.shp")
geo   = rd(f"{wdpa}/GEO/WDPA_WDOECM_Aug2026_Public_GEO_shp-polygons.shp")
wd = gpd.GeoDataFrame(pd.concat([arm, geo], ignore_index=True), crs=DST)
wd = gpd.clip(wd, sa)                       # clip WDPA to study area
prot = gpd.GeoDataFrame(geometry=list(amnp.geometry)+list(wd.geometry), crs=DST)  # AMNP + WDPA
borders = rd(f"{B}/data/ne/ne_10m_admin_0_boundary_lines_land.shp")  # country borders
from shapely.geometry import box
frame = box(xmin, ymin, xmax, ymax)
roads = gpd.clip(roads, frame)
borders = gpd.clip(borders, frame)
print("roads:", len(roads), "prot polys:", len(prot), "(amnp %d + wdpa %d)" % (len(amnp), len(wd)),
      "borders:", len(borders))

# ---- plot ----
fig, ax = plt.subplots(figsize=(12, 10.5))
ax.imshow(rgba, extent=[xmin, xmax, ymin, ymax], origin="upper", interpolation="nearest", zorder=1)
# protected areas: grey fill + prominent dark border
prot.plot(ax=ax, facecolor="#b7b7b7", edgecolor="#333333", alpha=0.42, linewidth=1.4, zorder=2)
prot.boundary.plot(ax=ax, color="#333333", linewidth=1.4, zorder=2.1)
roads.plot(ax=ax, color="#111111", linewidth=0.45, zorder=3)
# country borders: thick black, much heavier than roads
borders.plot(ax=ax, color="black", linewidth=3.4, zorder=3.5)
sa.boundary.plot(ax=ax, color="black", linewidth=2.4, zorder=4)

# dump symbol (small)
img = mpimg.imread(f"{B}/data/dump_symbol.png")
cen = dump.geometry.union_all().centroid
oi = OffsetImage(img, zoom=0.06)
ax.add_artist(AnnotationBbox(oi, (cen.x, cen.y), frameon=False, box_alignment=(0.5,0.5), zorder=6))

ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
ax.set_aspect("equal"); ax.set_axis_off()
ax.set_title("Study area and land cover - Northeastern Türkiye",
             fontsize=17, fontweight="bold", loc="left", pad=12)

# scale bar (correct for Web-Mercator distortion at central latitude)
yc = (ymin+ymax)/2
latc = math.degrees(2*math.atan(math.exp(yc/6378137.0)) - math.pi/2)
gkm = 50; bar = gkm*1000.0/math.cos(math.radians(latc))
xb = xmin+(xmax-xmin)*0.05; yb = ymin+(ymax-ymin)*0.045
ax.plot([xb, xb+bar], [yb, yb], color="black", lw=4, solid_capstyle="butt", zorder=7)
ax.plot([xb, xb+bar/2], [yb, yb], color="white", lw=4, solid_capstyle="butt", zorder=7)
ax.plot([xb, xb+bar], [yb, yb], color="black", lw=1.2, zorder=7)
ax.text(xb+bar/2, yb+(ymax-ymin)*0.012, f"{gkm} km", ha="center", va="bottom", fontsize=11)

# north arrow — top-right, just below the frame's top edge
xn = xmax-(xmax-xmin)*0.05
yb = ymax-(ymax-ymin)*0.135; yt = ymax-(ymax-ymin)*0.065
ax.annotate("", xy=(xn, yt), xytext=(xn, yb),
            arrowprops=dict(facecolor="black", edgecolor="black", width=5, headwidth=16), zorder=7)
ax.text(xn, ymax-(ymax-ymin)*0.038, "N", ha="center", va="center",
        fontsize=15, fontweight="bold", zorder=7)

# legend (with garbage-dump symbol rendered via a custom image handler)
from matplotlib.legend_handler import HandlerBase
class HandlerImage(HandlerBase):
    def __init__(self, im, zoom): super().__init__(); self.im=im; self.zoom=zoom
    def create_artists(self, legend, orig, xd, yd, width, height, fontsize, trans):
        oi = OffsetImage(self.im, zoom=self.zoom); oi.image.axes = legend.axes
        ab = AnnotationBbox(oi, (width/2.-xd, height/2.-yd), xycoords=trans,
                            frameon=False, box_alignment=(0.5, 0.5))
        return [ab]
class _Dump: pass
dh = _Dump()
handles = [Patch(fc=COL[2]), Patch(fc=COL[5]), Patch(fc=COL[7]), Patch(fc=COL[8]),
           Patch(fc=COL[11]), Patch(fc=COL[1]),
           Patch(fc="#b7b7b7", alpha=0.42, ec="#333333", lw=1.4),
           Line2D([0],[0], color="black", lw=3.4),
           Line2D([0],[0], color="#111111", lw=1.2),
           Line2D([0],[0], color="black", lw=2.4),
           dh]
labels = ["Forest","Crops","Built-up areas","Bare ground","Rangeland","Water",
          "Protected areas","Country border","Roads","Study area","Garbage dump"]
ax.legend(handles, labels, handler_map={dh: HandlerImage(img, zoom=0.10)},
          loc="center left", bbox_to_anchor=(1.005, 0.5),
          fontsize=12, frameon=True, borderpad=0.9, labelspacing=0.8, handleheight=1.6)

fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
print("SAVED", OUT)
