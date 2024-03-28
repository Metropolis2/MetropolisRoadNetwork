"""
This script takes as input a GeoDataFrame of edges and performs various operations to make the data
ready to use with METROPOLIS2.
"""
import sys
import os
import time
import tomllib

import numpy as np
import networkx as nx
import pandas as pd
import geopandas as gpd


def read_edges(input_file):
    if input_file.endswith(".parquet"):
        gdf = gpd.read_parquet(input_file)
    else:
        gdf = gpd.read_file(input_file)
    columns = [
        "geometry",
        "source",
        "target",
        "length",
        "speed",
        "lanes",
        "road_type",
    ]
    for col in columns:
        if not col in gdf.columns:
            print("Error: Missing column {}".format(col))
    return gdf


def set_default_values(gdf, config):
    # Set default speeds.
    if "urban" in gdf.columns:
        # Set default speeds based on urban vs rural areas.
        urban_speeds = config["default_speed"].get("urban")
        if not isinstance(urban_speeds, dict):
            raise Exception(
                "Missing or invalid table `postprocess_network.default_speed.urban` in config"
            )
        rural_speeds = config["default_speed"].get("rural")
        if not isinstance(rural_speeds, dict):
            raise Exception(
                "Missing or invalid table `postprocess_network.default_speed.rural` in config"
            )
        urban_speeds = pd.DataFrame(
            list(urban_speeds.values()),
            index=list(urban_speeds.keys()),
            columns=["urban_speed"],
        )
        rural_speeds = pd.DataFrame(
            list(rural_speeds.values()),
            index=list(rural_speeds.keys()),
            columns=["rural_speed"],
        )
        default_speeds = pd.concat((urban_speeds, rural_speeds), axis=1)
        gdf = gdf.merge(default_speeds, left_on="road_type", right_index=True, how="left")
        gdf.loc[gdf["speed"].isna() & gdf["urban"], "speed"] = gdf["urban_speed"]
        gdf.loc[gdf["speed"].isna() & ~gdf["urban"], "speed"] = gdf["rural_speed"]
        gdf = gdf.drop(columns=["urban_speed", "rural_speed"])
    else:
        speeds = config["default_speed"]
        if not isinstance(speeds, dict):
            raise Exception("Invalid table `postprocess_network.default_speed` in config")
        gdf["speed"] = gdf["speed"].fillna(gdf["road_type"].map(speeds))
    # Set default number of lanes.
    nb_lanes = config["default_nb_lanes"]
    if not isinstance(nb_lanes, dict):
        raise Exception("Invalid table `postprocess_network.default_nb_lanes` in config")
    gdf["lanes"] = gdf["lanes"].fillna(gdf["road_type"].map(nb_lanes))
    # Set default bottleneck capacity.
    capacities = config["default_capacity"]
    if not isinstance(capacities, dict):
        raise Exception("Invalid table `postprocess_network.default_capacity` in config")
    if "capacity" in gdf.columns:
        gdf["capacity"] = gdf["capacity"].fillna(gdf["road_type"].map(capacities))
    else:
        gdf["capacity"] = gdf["road_type"].map(capacities)
    return gdf


def remove_duplicates(gdf):
    """Remove the duplicates edges, keeping in order of priority the one in the main graph, with the
    largest capacity and with smallest free-flow travel time."""
    print("Removing duplicate edges")
    n = len(gdf)
    # Sort the dataframe.
    gdf["tt"] = gdf["length"] / (gdf["speed"] / 3.6)
    gdf.sort_values(["capacity", "tt"], ascending=[False, True], inplace=True)
    gdf.drop(columns="tt", inplace=True)
    # Drop duplicates.
    gdf.drop_duplicates(subset=["source", "target"], inplace=True)
    if n > len(gdf):
        print("Warning: discarding {} edges duplicated".format(n - len(gdf)))
    return gdf


def select_connected(gdf):
    print("Building graph...")
    G = nx.DiGraph()
    G.add_edges_from(
        map(
            lambda v: (v[0], v[1]),
            gdf[["source", "target"]].values,
        )
    )
    # Keep only the nodes of the largest strongly connected component.
    nodes = max(nx.strongly_connected_components(G), key=len)
    if len(nodes) < G.number_of_nodes():
        print(
            "Warning: discarding {} nodes disconnected from the largest graph component".format(
                G.number_of_nodes() - len(nodes)
            )
        )
        gdf = gdf.loc[gdf["source"].isin(nodes) & gdf["target"].isin(nodes)].copy()
    #  # We now do the same for the main graph.
    #  G = nx.DiGraph()
    #  G.add_edges_from(
        #  map(
            #  lambda v: (v[0], v[1]),
            #  gdf.loc[gdf["main_graph"], ["source", "target"]].values,
        #  )
    #  )
    #  connected_nodes = max(nx.strongly_connected_components(G), key=len)
    #  if len(connected_nodes) < G.number_of_nodes():
        #  n0 = gdf["main_graph"].sum()
        #  gdf["main_graph"] = (
            #  gdf["main_graph"]
            #  & gdf["source"].isin(connected_nodes)
            #  & gdf["target"].isin(connected_nodes)
        #  )
        #  n1 = gdf["main_graph"].sum()
        #  print(
            #  "Warning: removing {} edges from the main graph as they are not connected".format(
                #  n0 - n1
            #  )
        #  )
    return gdf


def reindex(gdf):
    gdf["id"] = np.arange(len(gdf))
    return gdf


def check(gdf, config):
    gdf["lanes"] = gdf["lanes"].clip(config.get("min_nb_lanes", 1))
    gdf["length"] = gdf["length"].clip(config.get("min_length", 0.0))
    gdf["speed"] = gdf["speed"].clip(config.get("min_speed", 1e-4))
    gdf["target_count"] = gdf.groupby("target")["id"].transform("count")
    return gdf


def clean(gdf, config):
    gdf = set_default_values(gdf, config)
    gdf = remove_duplicates(gdf)
    if config.get("ensure_connected", True):
        gdf = select_connected(gdf)
    gdf = reindex(gdf)
    gdf = check(gdf, config)
    return gdf


def save(gdf, output_file):
    print("Saving post-processed edges")
    directory = os.path.dirname(output_file)
    if not os.path.isdir(directory):
        os.makedirs(directory)
    if output_file.endswith("parquet"):
        gdf.to_parquet(output_file)
    elif output_file.endswith("geojson"):
        gdf.to_file(output_file, driver="GeoJSON")
    elif output_file.endswith("fgb"):
        gdf.to_file(output_file, driver="FlatGeobuf")
    elif output_file.endswith("shp"):
        gdf.to_file(output_file, driver="Shapefile")
    else:
        raise Exception(f"Unsupported format for output file: `{output_file}`")


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

    input_file = config.get("raw_edges_file")
    if input_file is None:
        raise Exception("Missing key `raw_edges_file` in config")
    if not os.path.exists(input_file):
        raise Exception(f"Raw edges file not found:\n`{input_file}`")
    if not "clean_edges_file" in config:
        raise Exception("Missing key `clean_edges_file` in config")

    if not "postprocess_network" in config:
        raise Exception("Missing key `postprocess_network` in config")
    if not "default_nb_lanes" in config["postprocess_network"]:
        raise Exception("Missing key `default_nb_lanes` in config")
    if not "default_capacity" in config["postprocess_network"]:
        raise Exception("Missing key `default_capacity` in config")
    if not "default_speed" in config["postprocess_network"]:
        raise Exception("Missing key `default_speed` in config")

    gdf = read_edges(input_file)
    gdf = clean(gdf, config["postprocess_network"])
    save(gdf, config["clean_edges_file"])

    print("Total running time: {:.2f} seconds".format(time.time() - t0))
