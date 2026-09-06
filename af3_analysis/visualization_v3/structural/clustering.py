"""
V3 Structural Clustering and MDS.

Uses pairwise structural distance matrix.
Clustering approaches:
- hierarchical clustering
- k-medoids

MDS/PCA for 2D structural-space representation.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.manifold import MDS
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import pairwise_distances


def hierarchical_clustering(
    distance_matrix: np.ndarray,
    *,
    n_clusters: Optional[int] = None,
    linkage: str = "average",
    metric: str = "precomputed",
) -> Dict[str, Any]:
    """
    Perform hierarchical clustering on distance matrix.

    Parameters
    ----------
    distance_matrix : (N, N) array
    n_clusters : int, optional
        Number of clusters. If None, determine automatically.
    linkage : str
        Linkage method: 'single', 'complete', 'average', 'ward'.
    metric : str
        Distance metric. Use 'precomputed' if matrix is already distances.

    Returns
    -------
    dict with:
        - labels: (N,) array of cluster assignments
        - n_clusters: int
        - method: str
    """
    n = distance_matrix.shape[0]

    if n == 0:
        return {
            "labels": np.array([]),
            "n_clusters": 0,
            "method": "hierarchical",
        }

    # Determine number of clusters
    if n_clusters is None:
        # Use elbow method or default to sqrt(N)
        n_clusters = max(1, int(np.sqrt(n)))

    # Ensure at least 2 clusters for visualization
    n_clusters = min(n_clusters, n)
    n_clusters = max(n_clusters, 1)

    # Perform clustering
    clustering = AgglomerativeClustering(
        n_clusters=n_clusters,
        linkage=linkage,
        metric=metric,
    )

    # AgglomerativeClustering with metric='precomputed' expects the full
    # (N, N) distance matrix, not the condensed vector (scipy's linkage
    # functions are the ones that consume condensed distances).
    clustering.fit(distance_matrix)

    labels = clustering.labels_

    return {
        "labels": labels,
        "n_clusters": n_clusters,
        "method": "hierarchical",
        "linkage": linkage,
    }


def kmeans_clustering(
    distance_matrix: np.ndarray,
    *,
    n_clusters: Optional[int] = None,
    n_init: int = 10,
    random_state: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Perform k-medoids-like clustering using KMeans on distance coordinates.

    Note: This is not true k-medoids. For true k-medoids, use sklearn_extra.

    Parameters
    ----------
    distance_matrix : (N, N) array
    n_clusters : int, optional
    n_init : int
    random_state : int, optional

    Returns
    -------
    dict with cluster assignments
    """
    n = distance_matrix.shape[0]

    if n == 0:
        return {
            "labels": np.array([]),
            "n_clusters": 0,
            "method": "kmeans",
        }

    if n_clusters is None:
        n_clusters = max(1, int(np.sqrt(n)))

    n_clusters = min(n_clusters, n)
    n_clusters = max(n_clusters, 1)

    # Convert distance matrix to coordinates using MDS
    mds = MDS(n_components=min(2, n - 1), dissimilarity="precomputed", random_state=random_state)
    coords = mds.fit_transform(distance_matrix)

    # Cluster in coordinate space
    kmeans = KMeans(
        n_clusters=n_clusters,
        n_init=n_init,
        random_state=random_state,
    )
    kmeans.fit(coords)
    labels = kmeans.labels_

    return {
        "labels": labels,
        "n_clusters": n_clusters,
        "method": "kmeans",
        "coordinates": coords,
    }


def mds_embedding(
    distance_matrix: np.ndarray,
    *,
    n_components: int = 2,
    dissimilarity: str = "precomputed",
    random_state: Optional[int] = None,
    normalized_stress: bool = True,
) -> Dict[str, Any]:
    """
    Perform MDS embedding of structural distance matrix.

    Parameters
    ----------
    distance_matrix : (N, N) array
    n_components : int
        Number of dimensions for embedding.
    dissimilarity : str
        'precomputed' if matrix is already distances.
    random_state : int, optional
    normalized_stress : bool
        Whether to use normalized stress.

    Returns
    -------
    dict with:
        - coordinates: (N, n_components) array
        - stress: float
        - n_components: int
    """
    n = distance_matrix.shape[0]

    if n == 0:
        return {
            "coordinates": np.array([]).reshape(0, n_components),
            "stress": np.nan,
            "n_components": n_components,
        }

    # Limit components
    n_components = min(n_components, n - 1)
    n_components = max(n_components, 1)

    # Perform MDS
    mds = MDS(
        n_components=n_components,
        dissimilarity=dissimilarity,
        random_state=random_state,
        normalized_stress=normalized_stress,
    )
    coordinates = mds.fit_transform(distance_matrix)

    return {
        "coordinates": coordinates,
        "stress": float(mds.stress_),
        "n_components": n_components,
        "method": "MDS",
    }


def pca_embedding(
    coordinates: np.ndarray,
    *,
    n_components: int = 2,
) -> Dict[str, Any]:
    """
    Perform PCA on coordinate data.

    Parameters
    ----------
    coordinates : (N, D) array
    n_components : int

    Returns
    -------
    dict with PCA results
    """
    from sklearn.decomposition import PCA

    n = coordinates.shape[0]
    n_components = min(n_components, min(n, coordinates.shape[1]))
    n_components = max(n_components, 1)

    pca = PCA(n_components=n_components)
    coords_transformed = pca.fit_transform(coordinates)

    return {
        "coordinates": coords_transformed,
        "explained_variance_ratio": pca.explained_variance_ratio_.tolist(),
        "components": pca.components_.tolist(),
        "n_components": n_components,
        "method": "PCA",
    }


def cluster_summary(
    labels: np.ndarray,
    conditions: List[str],
    seeds: List[int],
    predictions: List[str],
) -> Dict[str, Any]:
    """
    Generate summary of clustering results.

    Parameters
    ----------
    labels : (N,) array of cluster assignments
    conditions : list of condition_id
    seeds : list of seed
    predictions : list of prediction_id

    Returns
    -------
    dict with cluster composition
    """
    df = pd.DataFrame({
        "prediction_id": predictions,
        "condition": conditions,
        "seed": seeds,
        "cluster": labels,
    })

    # Count per cluster
    cluster_counts = df.groupby("cluster").size().to_dict()

    # Composition per cluster
    cluster_composition = {}
    for cluster_id in sorted(cluster_counts.keys()):
        cluster_df = df[df["cluster"] == cluster_id]
        composition = cluster_df.groupby("condition").size().to_dict()
        cluster_composition[cluster_id] = {
            "n": int(cluster_counts[cluster_id]),
            "composition": composition,
        }

    return {
        "n_clusters": len(cluster_counts),
        "cluster_counts": cluster_counts,
        "cluster_composition": cluster_composition,
        "n_observations": len(df),
    }
