# Atyuha - Optimal Path Planning - Vishwakarma 2024
<i>This is the repository for the Vishwakarma Awards 2024 on Optimal Path Planning for Multipe UAV Area Coverage.</i>

The folder `area-coverage` on the `main` branch contains the code for our algorithm on generating polygons optimized for number of sides and least redunant area. The folder `Tests`, along with others contains our iterations on other algorithms connected to optimizing for dividing a given area into sub-regions.

The folder `circle-fitting` contains the algorithm to generate minimum number of circles (and their respective centers - as nodes) which can be fit into any given rectangle. This will soon be extended to the case of a general quadrilateral. The nodes will be used as waypoints to generate optimal trajectories in order to cover the assigned area by the UAV.
