# Detecting transitions between collective motion regimes using functional hypothesis test with persistent homology data
This study (https://arxiv.org/abs/2403.03428) demonstrates how the significant changes based on shape or pattern features correlate with the existing transition events seen from samples of trajectory observations of experimentally collected data and synthetic data from cell migration model.

## Description
This software will provide simulation code to generate the particle positions at discrete time steps based on a first-order discrete-time stochastic model describing the motion of “cells” (particles) in a two-dimensional rectangular domain with periodic boundary conditions as introduced in https://doi.org/10.1039/D1SM00072A
Based on the collection motion model simulations, we verify the proposed method of using functional hypothesis test with persistent landscape curves

## Getting Started
By following these instructions, you can setup a local copy of the project for use in testing and development
A list of prerequisites R packages used within the project:
*	TDA: Statistical Tools for Topological Data Analysis (https://cran.r-project.org/web/packages/TDA/index.html)
*	TDAkit: Toolkit for Topological Data Analysis (https://cran.r-project.org/web/packages/TDAkit/index.html)
*	igraph (https://cran.r-project.org/web/packages/igraph/vignettes/igraph.html)

### Codes
* The following R codes can be used to generate samples of simulations for specific parameters using the self-propelled particle model and to obtain persistent homology descriptors (persistence landscape) based on particle point cloud at discrete time steps 
- The main files is Main_code.R and required function files are Grad.R, Interaction_force.R, Prop_force_file.R, Compute_homology.R, Landscape_plot.R, Plot_cell_by_frame.R, Plotcomparisonmatrix.R, Lp.R.

* The Stat_test_code.R is used for conducting Functional hypothesis test with Westfall-Young correction to compare persistence landscape contours of random vs observed data at each time step 
