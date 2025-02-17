# KEPSurvey
This repository contains the instances and code for all model combinations discussed in the paper "Operational research approaches and mathematical models for kidney exchange: A literature survey and empirical evaluation" by Mathijs Barkel, Rachael Colley, Maxence Delorme, David Manlove and William Pettersson.

All algorithms are coded in C++ and require the commercial solver Gurobi (we used version 10.0.3). The files have the following contents:
- main.h/cpp             - example of front-end code for solving a given instance using any of the ILP models that we consider;
- main_functions.h/cpp   - code containing our main functions;
- RCVF.h/cpp             - code for solving combinations of CF-CYCLE/HCF-CYCLE/PIEF-CYCLE and CF-CHAIN/PIEF-CHAIN using Reduced Cost Variable Fixing; and
- helper_functions.h/cpp - code containing miscellaneous simple/supportive functions.

Moreover, "Instances.zip" contains a JSON-file for each of our test instances. 
These files have a name of the format "uk_2019_splitpra_bandxmatch_pra0_pdd\_{x}\_{y}\_{z}", where:
- x in {0.05, 0.10, 0.20} specifies the number of NDDs |N| relative to the number of RDPs |R|, namely x = |N|/|R|;
- y in {50, 100, 200, 500, 750, 1000} specifies the number of RDPs; and
- z in {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} specifies the instance index.

For questions on the code, send an e-mail to mathijsbarkel3@gmail.com.
