# CLP-GRASP

This is the repository of the article titled "Integrating Multi-Drop, Load-Bearing, and Static-Mechanical-Equilibrium in Container Loading with a Multithreaded GRASP Algorithm".

## Instances

Three types of instances were used:
*ceschia_CS: Instances obtained from Ceschia, S., & Schaerf, A. (2013). Local search for a multi-drop multi-container loading problem. Journal of Heuristics, 19(2), 275–294. https://doi.org/10.1007/s10732-011-9162-6
*DATA: Instances obtained from Junqueira, L., Morabito, R., & Sato Yamashita, D. (2012). MIP-based approaches for the container loading problem with multi-drop constraints. Annals of Operations Research, 199(1), 51–75. https://doi.org/10.1007/s10479-011-0942-z
*MLBR: Instances obtained from Christensen, S. G., & Rousøe, D. M. (2009). Container loading with multi-drop constraints. International Transactions in Operational Research, 16(6), 727–743. https://doi.org/https://doi.org/10.1111/j.1475-3995.2009.00714.x

## Solutions

The packing pattern is composed by the position of the packed boxes. Each row indicates two corners of each packed box.

## Implementation

This is the C++ implementation of the GRASP algorithm used to solve the Container Loading Problem. This code has the following parameters:
* -ins: Name of the instance file
* -nThreads: Number of threads used for the execution
* -maxTime: Maximum execution time
* -estabilidad: Indicates how stability of the cargo is tackled. 0 for no stability, 1 for partial support, 2 for full support
* -juntar: Indicates how maximal spaces (MS) are treated. 0 for not merging nor expanding MS, 1 for merging but not expanding MS, 2 for merging and expanding MS
* -presion: Indicates wether the load-bearing constraint (LBC) is handled. 0 for not considering LBX, 1 for considering LBC.
* -multidrop: Indicates how the multi-drop constraint is handled. 0 for not considering it, 1 for considereing Visibility criterion, 2 for considering vertical separating walls.

Any question about this repository should be adress to Germán Pantoja via e-mail: germanpantojab@gmail.com
