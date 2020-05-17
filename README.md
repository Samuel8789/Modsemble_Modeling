# Modsemble
Darik A. O’Neil // Laboratory of Rafael Yuste // Columbia University

## Modsemble Analysis Pipeline 
Computationally-tractable and efficient conditional random field modeling of arbitrary graphical structures. Global optimum guarantee when max degree <= Olog(n), incredibly reliable in general up to at least max degree Olog(n^8). O(1/T) convergence rate guarantee.

![IMAGE of FS](https://github.com/darikoneil/Modsemble_Modeling/blob/master/FFC.jpg)

## REQUIREMENTS
- Linux
- MATLAB 2019b (+ Symbolic Links)
- Make (G++ & GCC) 
- Signal Processing Toolbox 
- Parallel Computing Toolbox
- Image Processing Toolbox
- Curve Fitting Toolbox
## Required Third-Party Dependencies:
- QPBO-v1.32 
  - Optimizing binary MRFs via extended roof duality. C. Rother, V. Kolmogorov, V. Lempitsky, and M. Szummer. 2007.
- GLMNet 
  - Regularization Paths for Generalized Linear Models via Coordinate Descent. Friedman J, Hastie T, and Tibshirani R. 2010. 
- MLE-STRUC
  - Bethe Learning of Graphical Models via MAP Decoding. K. Tang, N. Ruozzi, D. Belanger, and T. Jebara. 2016.
## For Further Reading
- First Modeling Efforts
  - Identification and Targeting of Cortical Ensembles. Luis Carrillo-Reid, Shuting Han, Ekaterina Taralova, Tony Jebara, Rafael Yuste. BioRxiv, 2017.
- Papers Utilizing Method
  - Controlling Visually Guided Behavior by Holographic Recalling of Cortical Ensembles. Luis Carrillo-Reid, Shuting Han, Weijan Yang, Alejandro Akrouh, Rafael Yuste. Cell, 2019
  - Triggering visually-guided behavior by holographic activation of pattern completion neurons in cortical ensembles. Luis Carrillo-Reid, Shuting Han, Weijan Yang, Alejandro Akrouh, Rafael Yuste. BioRxiv, 2018.
  - **Imprinting and Recalling Cortical Ensembles. Luis Carrillo-Reid, Weijan Yang, Yuki Bando, Darcy Peterka, Rafael Yuste. Science, 2017**
  
