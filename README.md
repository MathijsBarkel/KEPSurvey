# KEPSurvey
This repository contains the instances and code for all model combinations discussed in the paper "Operational research approaches and mathematical models for kidney exchange: A literature survey and empirical evaluation" by Mathijs Barkel, Rachael Colley, Maxence Delorme, David Manlove and William Pettersson.

All algorithms are coded in C++ and require the commercial solver Gurobi (we used version 10.0.3). The files have the following contents:
- main.h/cpp             - example of front-end code for solving a given instance using any of the ILP models that we consider;
- main_functions.h/cpp   - code containing our main functions;
- RCVF.h/cpp             - code for solving combinations of CF-CYCLE, HCF-CYCLE or PIEF-CYCLE with CF-CHAIN or PIEF-CHAIN using reduced cost variable fixing; and
- helper_functions.h/cpp - code containing miscellaneous simple/supportive functions.

The file main.cpp is meant to be self-explanatory, one simply needs to make the following ten choices:
1. Choose which instance to solve.
2. Decide whether to solve the problem hierarchically (first maximize number of transplants, then maximize weight).
3. Decide whether to solve the unweighted or weighted variant of the instance.
4. Decide whether to count the weights to the terminal vertex tau (representing donations to the DDWL or donors becoming bridge donors).
5. Determine the maximum cycle and chain length.
6. Determine which cycle and chain models to apply and/or which hybrid models to apply.
7. Determine if reduced cost variable fixing should be used.
8. Determine if the LP relaxation should be solved seperately.
9. Determine the time limit per run (in seconds).
10. Determine what information to print.

Moreover, "Instances.zip" contains a JSON-file for each compatibility graph. Our instances are defined by the chosen compatibility graph together with choices 2-5.
Our compatibility graphs are generated using the instance generator created by Delorme et al. (2022), which is available at https://wpettersson.github.io/kidney-webapp/#/generator.
More specifically, we used their ``SplitPRA BandXMatch PRA0'' profile, which was shown by these authors to create instances with similar characteristics to those found in historical UKLKSS datasets.
Each JSON file has a name of the format "uk_2019_splitpra_bandxmatch_pra0_pdd\_{x}\_{y}\_{z}", where:
- x in {0.05, 0.10, 0.20} specifies the number of NDDs |N| relative to the number of RDPs |R|, namely x = |N|/|R|;
- y in {50, 100, 200, 500, 750, 1000} specifies the number of RDPs; and
- z in {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} specifies the instance index.

The results of our computational experiments based on this code can be found on https://www.optimalmatching.com/kep-survey-2025/.

For questions on the code, send an e-mail to mathijsbarkel3@gmail.com.

References:
- Delorme, M., García, S., Gondzio, J., Kalcsics, J., Manlove, D., Pettersson, W., and Trimble, J. (2022). Improved instance generation for kidney exchange programmes. Computers & Operations Research, 141:105707.
