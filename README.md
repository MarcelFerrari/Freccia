Here's a short Markdown Readme for the GitHub page of the Freccia Eigensolver library:

---

# Freccia Eigensolver Library

![Freccia Eigensolver Logo](logo.svg)

# The Freccia Eigensolver Library
The Freccia Eigensolver library is a high-performance computational tool designed for efficient diagonalization of special symmetric real matrices. 
Currently, the following matrix structures are supported:
- Diagonal Plus Rank One (DPR1)
- Diagonal Plus Rank k (DPRK)
- Banded
- Arrowhead
- Broad Arrowhead

Additionally, it features algorithms for efficient eigensystem updates for many of the listed cases.

## Status

The project is very early in development and currently only supports double precision computations relying on the Eigen3 and the Intel oneAPI-MKL libraries.
Currently, only the Eigen C++ interface is supported, however Python bindings to support numpy are in development.
