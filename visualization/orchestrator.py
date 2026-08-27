# -*- coding: utf-8 -*-
"""
============================================================================
AF3 Analysis Visualization Orchestrator
============================================================================

This module orchestrates the loading of raw data and calling all 8 core plotting functions. 
It acts as the single entry point for visualization generation, ensuring proper 
data handling before passing necessary schema context to specialized plots like F8.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional
import pandas as pd
import os
import numpy as np

# Experiment design metadata (optional)
try:
    from af3_analysis.experiment_metadata import ExperimentDesign, load_experiment_design
except ImportError:
    ExperimentDesign = None  # type: ignore[assignment,misc]
    load_experiment_design = None  # type: ignore[assignment,misc]

# Import local modules that were just created/updated
from .core_plots import (
    plot_qc_completeness,
    plot_seed_distributions,
    plot_factorial_interaction,
    plot_effect_size_forest,
    plot_variability,
    plot_seed_trajectories,
    plot_ecdf_overlay,
    plot_metric_relationships
)

def load_data(tables_dir: Path):
    """Loads all required CSV data files from the specified directory."""
    print("Loading visualization prerequisites...")
    try:
        seed_aggregated_df = pd.read_csv(tables_dir / "seed_aggregated.csv")
        descriptive_stats_df = pd.read_csv(tables_dir / "descriptive_stats.csv")
        pairwise_comparisons_df = pd.read_csv(tables_dir / "pairwise_comparisons.csv")
        
        # Add condition_name to descriptive_stats_df if missing
        if 'condition_name' not in descriptive_stats_df.columns and 'condition_name' in seed_aggregated_df.columns:
            mapping = dict(zip(seed_aggregated_df['condition_id'], seed_aggregated_df['condition_name']))
            descriptive_stats_df['condition_name'] = descriptive_stats_df['condition_id'].map(mapping)
            
    except FileNotFoundError as e:
        print(f"ERROR: Required data file missing. Please ensure {e} exists in the project directory.")
        raise

    return {
        'seed_aggregated': seed_aggregated_df,
        'descriptive_stats': descriptive_stats_df,
        'pairwise_comparisons': pairwise_comparisons_df
    }


def generate_all_figures(
    tables_dir: str,
    figures_dir: str,
    *,
    metadata_path: Optional[str] = None,
    experiment_design: Optional["ExperimentDesign"] = None,
):
    """
    The main entry point function. Loads data and executes the 8 plotting routines in order.

    Parameters
    ----------
    tables_dir : str
        Directory containing the analysis CSVs.
    figures_dir : str
        Directory to write figures to.
    metadata_path : str, optional
        Path to ``experiment_metadata.json``.  Loaded automatically if
        *experiment_design* is not provided.
    experiment_design : ExperimentDesign, optional
        Pre-loaded experiment design metadata.  If *None* and
        *metadata_path* is given, the metadata is loaded from that path.
    """
    print("\n================ Starting Visualization Figure Generation ================")
    
    # 1. Load Data
    try:
        data = load_data(Path(tables_dir))
    except Exception as e:
        print(f"🚨 FATAL ERROR during data loading: {e}. Halting visualization process.")
        return None

    # 1b. Load experiment design metadata (optional)
    design = experiment_design
    if design is None and metadata_path is not None and load_experiment_design is not None:
        try:
            design = load_experiment_design(metadata_path)
            print(f"  Loaded experiment metadata: {design.experiment_id} "
                  f"({design.n_conditions} conditions, {design.n_attributes} attributes)")
        except Exception as e:
            print(f"  WARNING: Could not load experiment metadata from {metadata_path}: {e}")
            design = None

    fig_map = {} # Dictionary to store names/paths of generated figures for the report step
    all_figures_ready = True
    
    # Propagate condition_name from seed_aggregated to descriptive_stats if
    # missing.  The analysis stage groups by condition_id but may not carry
    # the human-readable name through.
    ds = data['descriptive_stats']
    sa = data['seed_aggregated']
    if 'condition_name' not in ds.columns and 'condition_name' in sa.columns:
        name_map = dict(zip(sa['condition_id'], sa['condition_name']))
        if 'condition_id' in ds.columns:
            ds['condition_name'] = ds['condition_id'].map(name_map)
        elif 'condition_id_' in ds.columns:
            # Flattened multi-index columns from agg(['mean','std'])
            ds['condition_name'] = ds['condition_id_'].map(name_map)
        # else: cannot map — plots that need it will skip gracefully
    
    # Infer schema dynamically
    seed_agg = data['seed_aggregated']
    
    # Curated global metrics (exclude redundant max/min/median and constant lengths)
    global_metrics = ['pLDDT_mean', 'pLDDT_min', 'pae_mean', 'contact_prob_mean', 'ranking_score']
    # Chain-level metrics are NOT included in the plotting schema because:
    # 1. Different conditions have different chains → not comparable
    # 2. Hundreds of chain names (A, AA, AB…, ZZ) × many conditions creates
    #    thousands of subplots that hang matplotlib during rendering.
    # Chain-level data remains in seed_aggregated.csv for downstream analysis.
    
    valid_metrics = [m for m in global_metrics if m in seed_agg.columns]
    
    # Build condition order from metadata if available, else empty (fallback in utils)
    condition_order: list = []
    if design is not None:
        condition_order = design.condition_names

    schema = {
        'metrics': valid_metrics,
        'factors': [],  # Factors now come from metadata, not hard-coded
        'condition_order': condition_order,
        'experiment_design': design,  # Pass metadata through for plots that need it
    }

    # 2. Execute Plots in Order (M1 through M8)
    try:
        plot_qc_completeness(data['seed_aggregated'], schema, Path(figures_dir))
        fig_map["QC_Completeness"] = "fig_qc_completeness.png"

        seed_dist_files = plot_seed_distributions(data['seed_aggregated'], schema, Path(figures_dir)) 
        fig_map["Seed_Distributions"] = seed_dist_files if seed_dist_files else ["fig_seed_distributions.png"]
        
        plot_factorial_interaction(data['descriptive_stats'], schema, Path(figures_dir))
        fig_map["Factorial_Interaction"] = "fig_factorial_interaction.png"

        plot_effect_size_forest(data['pairwise_comparisons'], schema, Path(figures_dir))
        fig_map["Effect_Size_Forest"] = "fig_effect_size_forest.png"
        
        plot_variability(data['descriptive_stats'], schema, Path(figures_dir)) 
        fig_map["Variability"] = "fig_variability.png"

        plot_seed_trajectories(data['seed_aggregated'], schema, Path(figures_dir)) 
        fig_map["Seed_Trajectories"] = "fig_seed_trajectories.png"

        plot_ecdf_overlay(data['seed_aggregated'], schema, Path(figures_dir))
        fig_map["ECDF_Overlay"] = "fig_ecdf_overlay.png"

        plot_metric_relationships(
            dataframe=data['seed_aggregated'], 
            schema=schema,
            save_path=Path(figures_dir)
        )
        fig_map["Scatter_Matrix"] = "fig_scatter_matrix.png"

    except Exception as e:
        print(f"\n❌ An error occurred during visualization generation: {e}")
        import traceback
        traceback.print_exc()
        all_figures_ready = False
    
    print("\n================ Visualization Generation Complete ================")
    if all_figures_ready:
        return fig_map # Return successfully generated figure names/paths
    else:
        return None

# Example usage if this file was run directly (should be called by orchestrator.py)
if __name__ == "__main__":
    print("--- Running visualization module dry-run ---")
    dummy_tables = "path/to/your/data/directory" 
    dummy_figures = "path/to/your/figure/directory"
    # Use the defined dummy schema for testing.
    if os.path.exists(Path(dummy_tables)):
         generate_all_figures(dummy_tables, dummy_figures)