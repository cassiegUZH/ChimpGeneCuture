# Social Network Analysis

This folder contains all the necessary code and data required to perform social network analysis on the chimps data.


## Data

These datasets can be found in the `Data` folder.

- graphs of IBD and NePRA genetic relationships.
- tool usage datasets


## Code
- `main_fig3_diagrams.ipynb`: code to generate the diagrams in Fig. 3 of the paper.
- `SM_figS9_diagrams.ipynb`: code to generate the diagrams in Fig. S9 of the paper.
- `SM_figS10_S11_diagrams.ipynb`: code to generate the diagrams in Fig. S10 and S11 of the paper.


### Network backboning and projection
- Apart from the standard `networkx` library, we have used the code from [Extracting the Multiscale Backbone of Complex Weighted Networks](1) to perform network backboning and projection. 

[1]: https://ieeexplore.ieee.org/document/9073533 "The Impact of Projection and Backboning on Network Topologies"
