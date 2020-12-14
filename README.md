prim-perm-pca.R is an R function for performing a PRIM-based analysis of 3D chromatin relocalizations between aligned 3D reconstructions.
The referent (first named) structure is rotated to principal component axes prior to applying PRIM as has been recommended.
The relocalization (Euclidean) distance between the two 3D structures serves as outcome while the three (rotated) coordinates constitute the covariates. 
The entire PRIM procedure is then re-applied to permuted relocalization values in order to obtain a null referent distribution to assess the significance 
of the originally identified relocalization regions, facilitated by output graphics for the top 3 regions. 
