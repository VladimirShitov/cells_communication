from typing import Optional

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.path import Path as MplPath
import matplotlib.patches as patches
import networkx as nx
import numpy as np


def normalize_vector(vector: np.array, normalize_to: float) -> np.array:
    """Make `vector` norm equal to `normalize_to`

    vector: np.array
        Vector with 2 coordinates
    normalize_to: float
        A norm of the new vector

    Returns
    -------
    Vector with the same direction, but with length normalized to `normalize_to`
    """
    vector_norm = np.linalg.norm(vector)

    return vector * normalize_to / vector_norm


def orthogonal_vector(point: np.array, width: float,
                      normalize_to: Optional[float] = None) -> np.array:
    """Get orthogonal vector to a `point`

    point: np.array
        Vector with x and y coordinates of a point
    width: float
        Distance of the x-coordinate of the new vector from the `point` (in orthogonal direction)
    normalize_to: Optional[float] = None
        If a number is provided, normalize a new vector length to this number

    Returns
    -------
    Array with x and y coordinates of the vector, which is orthogonal to the vector from (0, 0) to `point`
    """
    EPSILON = 0.000001

    x = width
    y = -x * point[0] / (point[1] + EPSILON)

    ort_vector = np.array([x, y])

    if normalize_to is not None:
        ort_vector = normalize_vector(ort_vector, normalize_to)

    return ort_vector


def draw_self_loop(
        point: np.array,
        ax: Optional[plt.Axes] = None,
        padding: float = 1.5,
        width: float = 0.3,
        plot_size: int = 10,
        linewidth: float = 0.2,
        color: str = "pink"
) -> plt.Axes:
    """Draw a loop from `point` to itself

    !Important! By "center" we assume a (0, 0) point. If your data is centered around a different points,
    it is strongly recommended to center it around zero. Otherwise, you will probably get ugly plots

    Parameters
    ----------
    point: np.array
        1D array with 2 coordinates of the point. Loop will be drawn from and to these coordinates.
    padding: float = 1.5
        Controls how the distance of the loop from the center. If `padding` > 1, the loop will be
        from the outside of the `point`. If `padding` < 1, the loop will be closer to the center
    width: float = 0.3
        Controls the width of the loop
    linewidth: float = 0.2
        Width of the line of the loop
    ax: Optional[matplotlib.pyplot.Axes]:
        Axis on which to draw a plot. If None, a new Axis is generated
    plot_size: int = 7
        Size of the plot sides in inches. Ignored if `ax` is provided
    color: str = "pink"
        Color of the arrow

    Returns
    -------
    Matplotlib axes with the self-loop drawn
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=(plot_size, plot_size))

    point_with_padding = padding * point
    ort_vector = orthogonal_vector(point, width, normalize_to=width)

    first_anchor = ort_vector + point_with_padding
    second_anchor = -ort_vector + point_with_padding

    verts = [point, first_anchor, second_anchor, point]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]

    path = MplPath(verts, codes)

    patch = patches.FancyArrowPatch(
        path=path,
        lw=linewidth,
        arrowstyle="-|>",
        color=color,
        alpha=0.5,
        mutation_scale=30  # arrowsize in nx.draw_networkx_edges()
    )
    ax.add_patch(patch)

    return ax


def draw_graph_edge(graph: nx.Graph, pos: dict, edge: tuple, edge_weight: float,
                    ax: plt.Axes, color: str, arc_radius: float = 0.2) -> plt.Axes:
    """Draw the given edge of the network

    Parameters
    ----------
    graph: nx.Graph
        A graph, edge of which you want to draw
    pos: dict
        A dictionary with the positions of the nodes. Keys are nodes of the `graph`,
        and values are their positions
    edge: tuple
        Edge in the `graph`. Tuple must contain 2 elements: the start and the end of the edge
    edge_weight: float
        Number describing how much weight does the edge have. It is used for the visualisation
    ax: plt.Axes
        Axis on which to draw a plot
    color: str
        Color of the node
    arc_radius: float = 0.2
        Edges are drawn as arcs, not as straight lines. This parameter controls the curvature of
        the `edge` arc

    Returns
    -------
    Matplotlib axes with the edge drawn
    """

    nx.draw_networkx_edges(
        graph,
        pos=pos,
        width=edge_weight,
        edgelist=[edge],
        alpha=0.5,
        edge_color=color,
        ax=ax,
        arrowsize=30,
        connectionstyle=f"arc3,rad={arc_radius}",
        node_size=1000
    )

    return ax


def graph_edges_weights(graph):
    return {edge: graph.edges[edge]["weight"] for edge in graph.edges}


def draw_graph_edges(graph: nx.Graph, pos: dict, ax: plt.Axes):
    """Draw graph edges so that edges in the opposite directions look differently

    Parameters
    ----------
    graph: nx.Graph
        Graph, edges of which you want to draw
    pos: dict
        Dictionary, where keys are nodes and values are their positions. Can be obtained
        through networkx layout algorithms (e. g. nx.circular_layout())
    ax: plt.Axes
        Axis on which the edges are drawn

    Returns
    -------
    Axis with the edges drawn
    """

    edge_weights = graph_edges_weights(graph)
    edges_to_draw = set(graph.edges)

    for edge in graph.edges:
        if edge not in edges_to_draw:
            continue

        if edge[0] == edge[1]:  # By default, networkx doesn't draw self loops correctly
            draw_self_loop(point=pos[edge[0]], ax=ax, linewidth=edge_weights[edge])
            continue

        draw_graph_edge(graph, pos, edge, edge_weight=edge_weights[edge], ax=ax, color="pink")
        edges_to_draw.remove(edge)

        # Edges between the same vertices look confusing, if they have the same style
        # So we draw such edges with different colors and curvature
        reverse_edge = (edge[1], edge[0])

        if reverse_edge in graph.edges and reverse_edge in edges_to_draw:
            draw_graph_edge(
                graph,
                pos,
                reverse_edge,
                edge_weight=edge_weights[edge],
                ax=ax,
                color="lightblue",
                arc_radius=0.1
            )
            edges_to_draw.remove(reverse_edge)

    return ax


def draw_chord_diagram(graph: Optional[nx.Graph] = None, matrix: Optional[pd.DataFrame] = None,
                       ax=None, plot_size=10) -> plt.Axes:
    """Draw weighted connections between the nodes of the `graph`"""

    if graph is None and matrix is None:
        raise ValueError("Please provide either `graph` or `matrix`")

    elif graph is None:
        graph = nx.DiGraph(matrix)

    if ax is None:
        fig, ax = plt.subplots(figsize=(plot_size, plot_size))

    pos = nx.circular_layout(graph, center=(0, 0))

    nx.draw_networkx_nodes(
        graph,
        pos=pos,
        node_color="springgreen",
        ax=ax,
        alpha=0.7,
        node_size=1000)

    draw_graph_edges(graph, pos, ax)

    nx.draw_networkx_labels(graph, pos=pos, font_weight="bold", ax=ax)

    ax.margins(0.1)  # Give some space for labels

    return ax


def draw_dotplot(df, figsize=(12, 6)):
    df["radius"] = np.sqrt(df["count"]) / 3
    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(df["index"], df["variable"], s=0)

    # plot data row-wise as text with circle radius according to count
    for _, row in df.iterrows():
        bbox_props = {"boxstyle": f'circle, pad = {row["radius"]}', "fc": "w", "ec": "r", "lw": 2}

        if row["count"]:
            ax.annotate(
                str(row["count"]),
                xy=(row["index"], row["variable"]),
                bbox=bbox_props,
                ha="center",
                va="center",
                zorder=2,
                clip_on=True
            )

    ax.grid(ls="--", zorder=1)
    fig.tight_layout()
    ax.set_title("Interactions counts", fontsize=18, fontweight="bold")

    return ax
