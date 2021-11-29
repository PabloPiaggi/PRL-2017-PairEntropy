# Enhancing Entropy and Enthalpy Fluctuations to Drive Crystallization in Atomistic Simulations
## Pablo Piaggi, Omar Valsson, and Michele Parrinello
### Phys. Rev. Lett. 119, 015701 (2017)

[![DOI](http://img.shields.io/badge/DOI-10.1103%2FPhysRevLett.119.015701-blue)](https://doi.org/10.1103/PhysRevLett.119.015701)
[![arXiv](http://img.shields.io/badge/arXiv-1612.03235-B31B1B.svg)](https://arxiv.org/abs/1612.03235)
[![plumID:21.047](https://www.plumed-nest.org/eggs/21/047/badge.svg)](https://www.plumed-nest.org/eggs/21/047/)

This repository contains the input files to reproduce the results of the paper mentioned above. 
The source code for the pair entropy CV is provided (PairEntropy.cpp) and can be used in PLUMED with the [LOAD keyword](https://www.plumed.org/doc-master/user-doc/html/_l_o_a_d.html).
This is a minimal input for the CV:
```
PAIRENTROPY ...
  LABEL=s2
  ATOMS=1-250
  MAXR=0.65
  SIGMA=0.01
... PAIRENTROPY
```

More information about the use of the CV can be found [here](https://sites.google.com/site/pablompiaggi/scripts/pair-entropy-cv).

Please e-mail me if you have trouble reproducing the results.
