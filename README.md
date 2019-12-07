These are the scripts for the calculations required for my final year project. Currently the project is still getting off the ground so there isn't a massive amount to talk about, but hopefully soon it should have some interesting scripts to do with calculating particle collisions, and possible even prototypes for contributions to the ARES project that CERN use for particle collision data processing.

To run, compile the C script as a shared library. This can be done with a simple execution of
`bash gcc`

From there you can run the python script directly which will write results to the `results` file and generate a graph at the end of computation. 

For now this contains a basic monte carlo integration script for numerically evaluating integrals from this paper: https://arxiv.org/pdf/1807.11487.pdf

The script depends on numpy and matplotlib.
