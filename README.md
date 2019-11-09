# StarMatch.jl

Solve the "lost in space" star tracker problem.

This package is a Julia implementation of Approach 3 (_Rotation Invariant Vector Frame_) from: [A reliable and fast lost-in-space mode star tracker](https://dr.ntu.edu.sg/handle/10220/47620) by Metha Deval Samirbhai. Matlab code is available at https://github.com/DevD1092/spr_riav.

## Description

The rotation invariant vector frame is a geometric technique that compares vector paths drawn between stars on the image with prebuilt vector paths drawn from a catalog and saved in a star pattern database (SPD). This technique does not make use of magnitude information and is robust to false and missing stars.

What are rotation invariant vectors? (picture)

Image solution overview:

  1. Locate the star nearest the image center.
  2. Sort neighbor stars by their distance from the center star.
  3. Build the rotation invariant vector path beginning from each of the `N_NEAREST` neighbor stars.
  4. Shortlist SPD by checking for entries that have a distance from center star to nearest neighbor star within `distancetolerance` from the corresponding distance in the image.
  5. For each of the candidate SPD entries, increment a vote counter for each SPD entry vector path component that is within `vectortolerance` from each of the image vector path components.
  6. The image star nearest the center is identified as the entry with the most votes.
  7. Vector path components for the SPD entry of the winning center star are compared to the image vector components to identify the image neighbor stars.

## Usage

### Generate Star Pattern Database (SPD)

Mention catalog file.

Demo generatespd

### Solve image

Demo solve


## Citing

Please cite Mehta's paper if used in scientific work:

  1. Mehta Deval Samirbhai, Shoushun Chen, and Kay-Soon Low, “A Rotation-Invariant Additive Vector Sequence Based Star Pattern Recognition”, IEEE Transactions on Aerospace and Electronic Systems, vol. 55, no. 2, pp. 689-705, April 2019. https://ieeexplore.ieee.org/document/8430566/

We also encourage you to cite this package. Refer to the Zenodo DOI at the top of the page or [CITATION.bib](CITATION.bib).
