"""
This script takes as input a .osm.pbf file of an area and returns a GeoDataFrame representing the
edges of the road network.
The GeoDataFrame has the following columns:
- `geometry` (EPSG:4326)
- `osm_id`
- `source`
- `target`
- `speed` (km/h, can be null)
- `length` (m)
- `lanes` (can be null)
- `urban` (bool, edge is in an urban area)
- `name` (can be null)
- `road_type` (id of the corresponding OSM highway tag)
"""
import os
import sys
import time
import tomllib

import numpy as np
import geopandas as gpd
import osmium
from osmium.geom import WKBFactory
import pyproj
from shapely.ops import transform
from shapely.geometry import LineString
from shapely.prepared import PreparedGeometry, prep


def valid_way(way, highways):
    """Returns True if the way is a valid way to consider."""
    has_access = not "access" in way.tags or way.tags["access"] in (
        "yes",
        "permissive",
        "destination",
    )
    return has_access and len(way.nodes) > 1 and way.tags.get("highway") in highways


def is_urban_area(area, urban_landuses):
    """Returns True if the area is an urban area."""
    return area.tags.get("landuse") in urban_landuses and (area.num_rings()[0] > 0)


class UrbanAreasReader(osmium.SimpleHandler):
    def __init__(self, urban_landuses):
        super().__init__()
        # Set of landuse tags to be considered as urban areas.
        self.urban_landuses = urban_landuses
        self.wkb_factory = WKBFactory()
        self.areas_wkb = list()

    def area(self, area):
        if not is_urban_area(area, self.urban_landuses):
            return
        self.handle_area(area)

    def handle_area(self, area):
        self.areas_wkb.append(self.wkb_factory.create_multipolygon(area))

    def get_urban_area(self):
        polygons = gpd.GeoSeries.from_wkb(self.areas_wkb)
        return polygons.unary_union


class NodeReader(osmium.SimpleHandler):
    def __init__(self, highways):
        super().__init__()
        # Set of highway tags to be considered.
        self.highways = set(highways)
        # Ids of the nodes explored.
        self.nodes_explored = set()
        # Ids of all the nodes in the final graph.
        self.node_ids = set()
        # Ids of all the edges in the final graph.
        self.edge_ids = set()

    def way(self, way):
        if not valid_way(way, self.highways):
            # Only consider valid highways.
            return
        self.edge_ids.add(way.id)
        # Always add source and target node to the final graph.
        self.node_ids.add(way.nodes[0].ref)
        self.node_ids.add(way.nodes[-1].ref)
        # Add the other nodes if they were already explored, i.e., they
        # intersect with another valid highway.
        for i in range(1, len(way.nodes) - 1):
            node = way.nodes[i]
            if node.ref in self.nodes_explored:
                self.node_ids.add(node.ref)
            else:
                self.nodes_explored.add(node.ref)


class EdgeReader(osmium.SimpleHandler):
    def __init__(self, node_ids, edge_ids):
        super().__init__()
        # Ids of the nodes in the final graph.
        self.node_ids = node_ids
        # Ids of the edges in the final graph.
        self.edge_ids = edge_ids
        # List of edges in the final graph, with their description.
        self.edges: list[dict] = list()

    def way(self, way):
        if not way.id in self.edge_ids:
            return

        road_type = way.tags["highway"]

        name = (
            way.tags.get("name", "") or way.tags.get("addr:street", "") or way.tags.get("ref", "")
        )
        if len(name) > 50:
            name = name[:47] + "..."

        oneway = (
            way.tags.get("oneway", "no") == "yes" or way.tags.get("junction", "") == "roundabout"
        )

        # Find maximum speed if available.
        maxspeed = way.tags.get("maxspeed", "")
        speed = None
        back_speed = None
        if maxspeed == "FR:walk":
            speed = 20
        elif maxspeed == "FR:urban":
            speed = 50
        elif maxspeed == "FR:rural":
            speed = 80
        else:
            try:
                speed = float(maxspeed)
            except ValueError:
                pass
        if not oneway:
            try:
                speed = float(way.tags.get("maxspeed:forward", "0")) or speed
            except ValueError:
                pass
            try:
                back_speed = float(way.tags.get("maxspeed:backward", "0")) or speed
            except ValueError:
                pass

        # Find number of lanes if available.
        lanes = None
        back_lanes = None
        if oneway:
            try:
                lanes = int(way.tags.get("lanes", ""))
            except ValueError:
                pass
            else:
                lanes = max(lanes, 1)
        else:
            try:
                lanes = (
                    int(way.tags.get("lanes:forward", "0")) or int(way.tags.get("lanes", "")) // 2
                )
            except ValueError:
                pass
            else:
                lanes = max(lanes, 1)
            try:
                back_lanes = (
                    int(way.tags.get("lanes:backward", "0")) or int(way.tags.get("lanes", "")) // 2
                )
            except ValueError:
                pass
            else:
                back_lanes = max(back_lanes, 1)

        for i, node in enumerate(way.nodes):
            if node.ref in self.node_ids:
                source = i
                break
        else:
            # No node of the way is in the nodes.
            print("Error: Found a valid way ({}) with no valid node".format(way.id))
            return

        j = source + 1
        for i, node in enumerate(list(way.nodes)[j:]):
            if node.ref in self.node_ids:
                target = j + i
                self.add_edge(
                    way,
                    source,
                    target,
                    oneway,
                    name,
                    road_type,
                    lanes,
                    back_lanes,
                    speed,
                    back_speed,
                )
                source = target

    def add_edge(
        self,
        way,
        source,
        target,
        oneway,
        name,
        road_type,
        lanes,
        back_lanes,
        speed,
        back_speed,
    ):
        source_id = way.nodes[source].ref
        target_id = way.nodes[target].ref
        if source_id == target_id:
            # Self-loop.
            return
        # Create a geometry of the road.
        coords = list()
        for i in range(source, target + 1):
            if way.nodes[i].location.valid():
                coords.append((way.nodes[i].lon, way.nodes[i].lat))
        geometry = LineString(coords)
        back_geometry = None
        if not oneway:
            back_geometry = LineString(coords[::-1])

        self.edges.append(
            {
                "geometry": geometry,
                "name": name,
                "road_type": road_type,
                "lanes": lanes,
                "speed": speed,
                "source": source_id,
                "target": target_id,
                "osm_id": way.id,
            }
        )

        if not oneway:
            self.edges.append(
                {
                    "geometry": back_geometry,
                    "name": name,
                    "road_type": road_type,
                    "lanes": back_lanes,
                    "speed": back_speed,
                    "source": target_id,
                    "target": source_id,
                    "osm_id": way.id,
                }
            )

    def post_process(self, urban_area: PreparedGeometry | None, metric_crs):
        edges = gpd.GeoDataFrame(self.edges, crs="EPSG:4326")

        if not urban_area is None:
            edges["urban"] = [urban_area.contains(geom) for geom in edges.geometry]

        # Compute length.
        edges["length"] = edges.to_crs(metric_crs).geometry.length

        print("Number of edges: {}".format(len(edges)))

        edges = edges[
            [
                "geometry",
                "source",
                "target",
                "length",
                "speed",
                "lanes",
                "urban",
                "osm_id",
                "name",
                "road_type",
            ]
        ].copy()

        edges["id"] = np.arange(len(edges))

        self.edges_df = edges

    def write_edges(self, output_file):
        directory = os.path.dirname(output_file)
        if not os.path.isdir(directory):
            os.makedirs(directory)
        if output_file.endswith("parquet"):
            self.edges_df.to_parquet(output_file)
        elif output_file.endswith("geojson"):
            self.edges_df.to_file(output_file, driver="GeoJSON")
        elif output_file.endswith("gpkg"):
            self.edges_df.to_file(output_file, driver="GPKG")
        elif output_file.endswith("fgb"):
            self.edges_df.to_file(output_file, driver="FlatGeobuf")
        elif output_file.endswith("shp"):
            self.edges_df.to_file(output_file, driver="Shapefile")
        else:
            raise Exception(f"Unsupported format for output file: `{output_file}`")


def buffer(geom, distance, metric_crs):
    wgs84 = pyproj.CRS("EPSG:4326")
    metric_crs = pyproj.CRS(metric_crs)
    project = pyproj.Transformer.from_crs(wgs84, metric_crs, always_xy=True).transform
    inverse_project = pyproj.Transformer.from_crs(metric_crs, wgs84, always_xy=True).transform
    metric_geom = transform(project, geom)
    buffered_geom = metric_geom.buffer(distance).simplify(0, preserve_topology=False)
    geom = transform(inverse_project, buffered_geom)
    return geom


if __name__ == "__main__":
    t0 = time.time()

    if len(sys.argv) == 1:
        # Read the config from the default path.
        config_path = "config.toml"
    else:
        if len(sys.argv) == 3 and sys.argv[1] == "--config":
            config_path = sys.argv[2]
        else:
            raise SystemExit(f"Usage: {sys.argv[0]} [--config <path_to_config.toml>]")

    if not os.path.exists(config_path):
        raise Exception(f"Cannot find config file `{config_path}`")

    with open(config_path, "rb") as f:
        try:
            config = tomllib.load(f)
        except Exception as e:
            raise Exception(f"Cannot parse config:\n{e}")

    if not "osm" in config:
        raise Exception("Missing key `osm` in config")
    if not "crs" in config:
        raise Exception("Missing key `crs` in config")

    input_file = config["osm"].get("input_file")
    if input_file is None:
        raise Exception("Missing key `osm.input_file` in config")
    if not os.path.exists(input_file):
        raise Exception(f"OSM input file not found:\n`{input_file}`")
    if not "raw_edges_file" in config:
        raise Exception("Missing key `raw_edges_file` in config")

    highways = config["osm"].get("highways")
    if not isinstance(highways, list):
        raise Exception("Invalid or missing array `osm.highways` in config")

    urban_landuses = config["osm"].get("urban_landuse")
    if not isinstance(urban_landuses, list):
        print("Warning: Invalid or missing array `osm.urban_landuse` in config")
        print("The urban column will not be created")

    print("Finding nodes...")
    node_reader = NodeReader(highways)
    node_reader.apply_file(input_file, locations=True, idx="flex_mem")

    print("Reading edges...")
    edge_reader = EdgeReader(node_reader.node_ids, node_reader.edge_ids)
    edge_reader.apply_file(input_file, locations=True, idx="flex_mem")

    if urban_landuses is None:
        urban_area = None
    else:
        print("Finding urban areas...")
        area_reader = UrbanAreasReader(urban_landuses)
        area_reader.apply_file(input_file, locations=True, idx="flex_mem")
        urban_area = area_reader.get_urban_area()

        # Buffer the urban areas by 50 meters to capture all nearby roads.
        urban_area = buffer(urban_area, 50, config["crs"])
        urban_area = prep(urban_area)

    print("Post-processing...")
    edge_reader.post_process(urban_area, config["crs"])

    print("Writing edges...")
    edge_reader.write_edges(config["raw_edges_file"])

    print("Done!")

    print("Total running time: {:.2f} seconds".format(time.time() - t0))
