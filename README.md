This project is reimplementing various monte carlo integrations that calculate the Jet distributions of particle collisions. The reimplementation aims to address integrable singularities contained in the equations being evaluated as to reduce error in the results, as well as eventually add basic multi-threading support and a higher level of compiler optimization, which the standard implementation found in the ARES framework does not do. Once this has been achieved, this code will be integrated into the ARES framework to be made avalible to any particle physics researches who need it. 

Currently, an arbitrary number of gluons can be calculated as well as their most significant correction, which have been combined into a total monte carlo integration for calculating any addative observable. The script will soon add all required corrections for the relevant collision, which after having corrected for Fcorrel should be trivial. After this the script will also be extended to non-addative observables.

To build each script with the appropiate optimizations, run `make <script name>`
Each script will write it's results to a file which can be read by `parseResults.py <file name>` in order to be plotted on a graph in matplotlib

Equations and methods based off the following papers:

https://arxiv.org/pdf/1412.2126.pdf

https://arxiv.org/pdf/1807.11487.pdf

Project supervised by Andrea Banfi
