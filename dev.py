import h5py
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, LineString
import pandas as pd

# === CONFIGURATION ===
hdf_path = "path/to/your/output.hdf"
shapefile_path = "path/to/reference_line.shp"
flow_area_name = "Your2DAreaName"

# === LOAD REFERENCE LINE ===
gdf = gpd.read_file(shapefile_path)
line_geom = gdf.geometry.unary_union
if not isinstance(line_geom, LineString):
    raise ValueError("Expected LineString geometry for reference line.")

# === LOAD MESH FROM HDF ===
with h5py.File(hdf_path, "r") as f:
    base_geom = f[f"/Geometry/2D Flow Areas/{flow_area_name}"]
    points = base_geom["Points"][:]               # Nx2
    faces = base_geom["Faces"][:]                 # MxN

    # Optional: Read center X/Y if needed
    face_x = base_geom["Face X Coordinates"][:]
    face_y = base_geom["Face Y Coordinates"][:]

    # === LOAD DEPTH AND VELOCITY ===
    summary_path = f"/Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/2D Flow Areas/{flow_area_name}"
    depth = f[summary_path + "/Maximum Depth"][:]
    try:
        velocity = f[summary_path + "/Maximum Velocity"][:]
    except KeyError:
        # Fallback to face velocity if needed
        velocity = f[summary_path + "/Maximum Face Velocity"][:]

# === BUILD FACE POLYGONS ===
face_polygons = []
for face in faces:
    face_coords = points[face]
    face_coords = face_coords[face_coords[:, 0] != -1]  # remove -1 padding
    if len(face_coords) >= 3:
        poly = Polygon(face_coords)
        face_polygons.append(poly if poly.is_valid else None)
    else:
        face_polygons.append(None)

# === INTERSECT REFERENCE LINE WITH FACES ===
data = []
for i, poly in enumerate(face_polygons):
    if poly is None or not poly.is_valid:
        continue
    intersection = line_geom.intersection(poly)
    if intersection.is_empty:
        continue

    length = 0
    if intersection.geom_type == "MultiLineString":
        length = sum(seg.length for seg in intersection.geoms)
    elif intersection.geom_type == "LineString":
        length = intersection.length

    if length > 0:
        dep = depth[i]
        vel = velocity[i]
        unit_discharge = dep * vel
        seg_discharge = unit_discharge * length
        data.append({
            "Face Index": i,
            "Width (ft)": length,
            "Centroid X": poly.centroid.x,
            "Centroid Y": poly.centroid.y,
            "Depth (ft)": dep,
            "Velocity (ft/s)": vel,
            "Unit Discharge (cfs/ft)": unit_discharge,
            "Segment Discharge (cfs)": seg_discharge
        })

# === SAVE TO CSV ===
df = pd.DataFrame(data)
df.to_csv("reference_line_flow_data.csv", index=False)

print(f"Extracted flow data for {len(df)} intersected faces.")
