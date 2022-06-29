#!/bin/bash

# prepare PhyCLIP input file and identify clusters in overview tree

# from the wiki:
# PhyCLIP’s user-defined parameters (S, gamma, FDR) can be calibrated across a 
# range of input values to optimise the statistical properties of the clustering
# results so as to select an optimal parameter set and its associated clustering 
# result.
rm -r phyclip_clusts
mkdir phyclip_clusts
cd phyclip_clusts

echo '../RAxML_bestTree.small_core' > phyclip_params.txt
echo '3-15(1),0.05-0.35(0.05),1-5(0.5)' >> phyclip_params.txt


# options to consider:
# Clean-up options:
#   --subsume_sensitivity_induced_clusters {0,1}
#                         Subsume cluster-size sensitivity-induced clusters into
#                         parent cluster (default = 1).
#   --sensitivity_percentile SENSITIVITY_PERCENTILE
#                         Percentile of cluster size distribution under which a
#                         cluster is considered to be sensitivity-induced
#                         (advanced option, default = 25%).
#   --subsume_subclusters {0,1}
#                         Subsume sub-clusters into their respective parent
#                         clusters (default = 0).
#   --force {0,1}         Force-cluster low information putative clusters.

# fine grained clusters
# phyclip.py --input_file phyclip_params.txt --optimise high --collapse_zero_branch_length 1 --equivalent_zero_length 1.000001e-06 --threads 40

# more broad clusters, not subsumed.
# phyclip.py --input_file phyclip_params.txt --optimise intermediate --collapse_zero_branch_length 1 --equivalent_zero_length --threads 40

# more broad clusters
phyclip.py --input_file phyclip_params.txt --optimise intermediate --collapse_zero_branch_length 1 --subsume_subclusters 1 --equivalent_zero_length 1.000001e-06 --threads 40 --tree_outgroup midpoint

  
