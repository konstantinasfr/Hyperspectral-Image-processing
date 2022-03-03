# Clustering: Hyperspectral-Image-processing

## Aim of this project: Identification of homogeneous regions in the Salinas HSI.

The goal of this study is to compare the performance of (a) cost function optimization clustering
algorithms (the k-means, the fuzzy c-means, the possibilistic c-means and the
probabilistic (where each cluster is modelled by a normal distribution) clustering algorithms) on the one hand and (b) the hierarchical algorithms (Complete-link,
WPGMC, Ward algorithms) on the other hand, in finding homogeneous regions in the
Salinas HSI, focusing ONLY on the pixels for which the class label information is
available.

## Implementaton

For the purpose of this study we will use Hyperspectral images (HSIs) that depict a
specific scene at several (L) narrow continuous spectral bands. This images can be
represented by a MxNxL three-dimensional cube, where the first two dimensions correspond to the spatial information, while the third corresponds to the spectral information. In our case the dimensions are 150x150x204, so the specific scene is depicted
by 204 spectral bands. The true labels of every sub-region is given so we will be able to
extract further conclusions about the accuracy of every clustering algorithm that will
be used.
With the use of principal component analysis we will reduse the dimensions of our
problem so that it will be more fast and easy to handle.
### Clustering algorithms used:
* k-means
* Fuzzy
* Probabilistic
* Agglomeratives
     * Single Link
     * Complete Link
     * Weighted Pair Group Method with Arithmetic Mean
     * Unweighted pair group method with arithmetic mean
     * Weighted pair group method centroid
     * Ward or minimum variance algorithm agglomerative
 
 ## Conclusions
 In conclusion the best agglomerative algorithm appears to be the Ward algorithm as
it outperforms the others. However in general the CFO algorithms gave good results
too, as we saw in case of fuzzy with 10 clusters and k-means with 9 clusters. Probably
possibilistic is able to give even better results if we perform a most careful tuning. It
is worth noting that the CFO algorithms were more noisy. Finally, the CFO algorithms
gave more mediocre results in contrast to the agglomerative ones which gave quite
bad results, but in the end both categories of algorithms managed to solve the problem
efficiently.
   
