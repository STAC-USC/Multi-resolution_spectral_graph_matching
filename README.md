# Multi-resolution spectral graph matching
This function performs a graph matching algorithm explained in:
"Spectral Graph Matching through a Multi-Resolution Approach" by
Victor González Navarro and Antonio Ortega

If you use the code, please include the following citation.

V. González and A. Ortega, "Multi-Resolution Spectral Graph Matching," 2019 IEEE International Conference on Image Processing (ICIP), Taipei, Taiwan, 2019, pp. 2319-2323.

Paper link: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8803306

To use this code you need to download and install the following library:
https://epfl-lts2.github.io/gspbox-html/ to plot the graphs. We have successfully run this code on Matlab 2017.

If graph data is available and loaded, the algorithm can be directly run in 'MAIN.m' file.
If there is no data, 'generate_graph.m' file can be used to generate a pair of graphs.
The rest of files are secondary function and must be place in the same fold as 'MAIN.m' file.
Other graph matching algorithms we use to compare in the paper are in 'Other S.G.M algorithms' foler.
