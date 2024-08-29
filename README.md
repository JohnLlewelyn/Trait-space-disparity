# Trait-space-disparity

Title: Trait-space disparity in fish communities spanning 380 million years from the Late Devonian to present

Authors: John Llewelyn, John A. Long, Richard Cloutier, Alice M. Clement, Giovanni Strona, Frédérik Saltré, Michael S. Y. Lee, Brian Choo, Kate Trinajstic, Olivia Vanhaesebroucke, Corey J. A. Bradshaw

In this study, we compare trait space of modern and Devonian fish communities to test whether species have dealt with their environment using similar strategies through time. Fish trait data was collected from databases, scientific literature, photographs, and expert opinion. Missing data was inferred using multiple imputation. The complete trait data set was then used to calculate Gower distances, quantifying differences between fish species in trait space. Principal coordinate analysis was applied to these Gower distances to ordinate species in lower dimensional space, and the species’ coordinates from the PCoA analysis were used to calculate functional diversity metrics, including: functional richness, functional nearest neighbour, functional specialisation, distance between community trait space centroids, and overlap between community trait spaces (Jaccard Index). The data and R code for replicating our analyses are provided in this repository. In the code, search for ‘###’ to find lines where file paths need to be specified/updated. The code files are numbered in the order they should be run.

Repository Structure

Data Folder
This folder contains five subfolders:

1.	taxonomy/
  - Contains an Excel file (.xlsx) with taxonomic information for the Devonian fish species.
2.	Devonian_fish_traits/
  - Contains an Excel file (.xlsx) of Devonian fish traits, excluding traits inferred through multiple imputation.
3.	modern_fish_traits_RDSs/
  - Contains RDS files of species traits for each modern community.
4.	Devonian_fish_traits_RDS/
  - Contains an RDS file with the Devonian species' traits, including imputed data.
5.	data_for_plots/
  - Contains RDS files produced by the code (below) and used to generate the study’s figures.

Code Folder
This folder contains five subfolders, organized by the specific steps of the analysis:

1.	1_multiple_imputation/
    Contains one file:
    - 1_missForest_taxonomyHiGitHub.R: Imputes missing trait data for Devonian fish.
2.	2_3_tidy_data/
    Contains two files:
    - 2_Devonian_combine_and_tidy_modern_RDS_files.R: Tidies the modern fish data.
    - 3_Devonian_tidy_Gogo_Miguasha.R: Tidies the Devonian fish data.
3.	4_5_6_mFD/
    Contains three files:
    - 4_gawdis_and_mFD_allSP.R: Calculates Gower distances, applies PCoA, and computes community functional diversity metrics (functional richness, nearest neighbour, specialization).
    - 5_gawdis_and_mFD_subsampling.R: Performs the same calculations as 4_gawdis_and_mFD_allSP.R, but subsamples communities to control for species diversity.
    - 6_gawdis_and_mFD_subsampling_plots.R: Plots the functional diversity metrics calculated in the previous two files.
4.	7_8_hypervolume/
    Contains two files:
    - 7_Hypervolumes_fixedBandwidths.R: Uses PCoA coordinates to fit hypervolumes describing each community’s trait space and calculates distances between trait space centroids and overlap (Jaccard Index).
    - 8_combine_heatmaps_and_within_vs_between.R: Plots the hypervolume results.
5.	9_traits_by_site/
    Contains one file:
    - 9_plot_traits_by_sites_updated2.R: Plots the distribution of traits across the two Devonian and six modern communities.

