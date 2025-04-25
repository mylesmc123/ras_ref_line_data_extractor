# %%
import h5py
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon, LineString, MultiLineString
from shapely.ops import unary_union, split
import pandas as pd
import matplotlib.pyplot as plt

# %%
# === CONFIGURATION ===
hdf_path = r"C:\Users\MBMcmanus\OneDrive - Garver\Documents\Work\ArkansasRiverWest\RAS_ArkansasRiverWest\Arkansas2019.p17.hdf"
shapefile_path = r"L:\2021\21T01042 - ARDOT 080676 Arkansas River - West\Design\Calculations\H&H\Hydraulics\ArkansasRiverWest\Reference_Lines.shp"

# === LOAD REFERENCE LINES ===
ref_gdf = gpd.read_file(shapefile_path)
ref_gdf

# %%

def get_first_flow_area(hdf_path):
    # Open the HDF5 file and get the geometry group
    with h5py.File(hdf_path, "r") as f:
        # Known non-area keys to ignore
        ignore_keys = {"Attributes", "Cell Info", "Cell Points", 
                    "Polygon Info", "Polygon Parts", "Polygon Points"}
        # Get the first flow area name from the HDF5 file that is not in the ignore keys
        flow_areas = [k for k in f["/Geometry/2D Flow Areas"].keys() if k not in ignore_keys]
        if not flow_areas:
            raise ValueError("No flow areas found in the HDF5 file.")
        
        return flow_areas[0]  # Return the first flow area name

def load_mesh_from_ras_66(hdf_path, flow_area_name, crs):
    with h5py.File(hdf_path, "r") as f:
        group = f[f"/Geometry/2D Flow Areas/{flow_area_name}"]

        # Coordinates of all face points (Nx2 array)
        coords = group["FacePoints Coordinate"][:]  # shape (N, 2)

        # Face point indexes for each cell (pairs of indices)
        faces_indexes = group["Faces FacePoint Indexes"][:]  # shape (M, 2)

        # Create a list to hold the edges (line segments)
        lines = []
        face_ids = []

        # Loop through each face and create LineString for each edge
        for i, face in enumerate(faces_indexes):
            # Get the two indices (start and end points)
            valid_indices = face[face >= 0]  # Ignore invalid (-1) indices
            
            if len(valid_indices) == 2:
                # Extract the coordinates for the two points (start and end of the edge)
                start_point = coords[valid_indices[0]]
                end_point = coords[valid_indices[1]]

                # Create a LineString for the edge (face)
                line = LineString([start_point, end_point])
                lines.append(line)
                face_ids.append(f"Edge_{i}")

    # Create a GeoDataFrame from the edges
    mesh_gdf = gpd.GeoDataFrame({"FaceID": face_ids}, geometry=lines, crs=crs)
    return mesh_gdf


from shapely.ops import unary_union, split
from shapely.geometry import MultiLineString, LineString
import geopandas as gpd

def split_lines_by_polygon_boundaries(lines_gdf, polygon_gdf):
    """
    Split lines in lines_gdf using the boundaries of polygons in polygon_gdf (polygons' boundaries are used for splitting).
    Returns a new GeoDataFrame with the split lines.
    """
    results = []

    # Step 1: Union all polygons into a single geometry
    merged_polygon = unary_union(polygon_gdf.geometry)

    for line_idx, line_row in lines_gdf.iterrows():
        line_geom = line_row.geometry
        split_result = line_geom

        # Step 2: Split the line by the union of all polygon boundaries
        if merged_polygon.intersects(split_result):
            try:
                split_result = split(split_result, merged_polygon)
                # If the result is a MultiLineString, ensure it's properly processed
                if isinstance(split_result, MultiLineString):
                    split_result = split_result.geoms
                elif isinstance(split_result, LineString):
                    split_result = [split_result]
            except Exception as e:
                print(f"Error splitting line {line_idx}: {e}")
                continue

        # Step 3: Store the split lines as separate records
        for part in split_result:
            if isinstance(part, (LineString, MultiLineString)):
                results.append({**line_row.drop(labels='geometry'), 'geometry': part})

    # Step 4: Check if results are empty and handle it
    if results:
        return gpd.GeoDataFrame(results, crs=lines_gdf.crs)
    else:
        print("No valid splits found.")
        return gpd.GeoDataFrame(columns=lines_gdf.columns, crs=lines_gdf.crs)
    
# %%
# === LOAD MESH FROM HDF ===
try:
    with h5py.File(hdf_path, "r") as f:
        # get the first flow area name
        flow_area = get_first_flow_area(hdf_path)     
except OSError as e:
    print(f"Error opening HDF5 file: {e}")  

# %%
# === LOAD MESH FROM HDF FOR THE FLOW AREA===
mesh_gdf = load_mesh_from_ras_66(hdf_path, flow_area, crs='PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["false_easting",0.0],PARAMETER["false_northing",0.0],PARAMETER["central_meridian",-96.0],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["latitude_of_origin",23.0],UNIT["Foot_US",0.3048006096012192]]')
print (f"Loaded {len(mesh_gdf)} polygons from the mesh.")
# %%
# plot the polygons
mesh_gdf.plot()
# %%
# save the polygons to a gpkg file
mesh_gdf.to_file(r"C:\Users\MBMcmanus\OneDrive - Garver\Documents\Work\ArkansasRiverWest\GIS\mesh.gpkg", driver='gpkg', layer='mesh', index=False)

# %%
ref_gdf = ref_gdf.to_crs(mesh_gdf.crs)  # Reproject reference lines to match mesh CRS
ref_gdf

# %%
# === INTERSECT REFERENCE LINE WITH MESH ===
segments_gdf = split_lines_by_polygon_boundaries(ref_gdf, mesh_gdf)
segments_gdf
# %%
print("Before split:", len(ref_gdf))
print("After split:", len(segments_gdf))
# %%
# get the first line from ref_gdf
line = ref_gdf.iloc[0].geometry

# Initialize the result with the original line
split_result = [line]

# Loop through each polygon in mesh_gdf
for polygon in mesh_gdf.geometry:
    boundary = polygon.boundary  # Get the boundary of the polygon
    
    # Check if the line intersects the polygon boundary
    if split_result[0].intersects(boundary):  
        # If the line intersects the boundary, split it
        temp_split = []
        for geom in split_result:
            temp_split.extend(split(geom, boundary))  # Split each segment
        split_result = temp_split  # Update the result with the new segments

# Print the results
for geom in split_result:
    print(geom)
# %%
# plot the split lines
fig, ax = plt.subplots()
for geom in split_result:
    x, y = geom.xy
    ax.plot(x, y, color='blue')
ax.set_title("Split Lines")
ax.set_xlabel("X Coordinate")
ax.set_ylabel("Y Coordinate")
plt.show()
# %%
