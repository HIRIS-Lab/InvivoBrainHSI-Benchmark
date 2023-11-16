# Utils 

This folder contains utility functions and scripts that provide additional functionality and support for the main components of the project.

## Contents

1. **[generateColorMap](generateColorMap.m):**
   - This code generate a color map based in the results obtained from the HELICoiD classification result. It assigns colors depending on a class list and the label in each pixel of the classification map.

2. **[hierclust2nmfMulti](hierclust2nmfMulti.m):**
   - Hierarchical custering based on rank-two nonnegative matrix factorization applied to hyperspectral image clustering.

3. **[knnFilter_window](knnFilter_window.m):**
   - Applies a k-Nearest Neighbors (kNN) filter to enhance an image.
     
4. **[majorityVoting](majorityVoting.m):**
   - Performs majority voting to generate an output map.

5. **[plotExample](plotExample.m):**
   - Plots various maps for visualization.
  
6. **[labelConsistentHierarchy](hkm/labelConsistentHierarchy.m):**
   - labels classes consistently with the same label across the hierarchy iteratively checks between parent and child if they share the same class and propogates the label till root.

7. **[splitclust](hkm/splitclust.m):**
   - Given a matrix M, split its columns into two subsets.
  
