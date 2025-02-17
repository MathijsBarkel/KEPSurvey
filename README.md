# KEPSurvey
This repository contains the instances and code for all model combinations discussed in the paper "Operational research approaches and mathematical models for kidney exchange: A literature survey and empirical evaluation" by Mathijs Barkel, Rachael Colley, Maxence Delorme, David Manlove and William Pettersson.

All algorithms are coded in C++ and part of our methods require the commercial solver Gurobi (we used version 10.0.3). The files have the following contents:

main.h/cpp             | example of front-end code for solving a given instance using any of the ILP models that we consider
main_functions.h/cpp   | code containing our main functions
RCVF.h/cpp             | code for solving combinations of CF-CYCLE/HCF-CYCLE/PIEF-CYCLE and CF-CHAIN/PIEF-CHAIN using Reduced Cost Variable Fixing 
helper_functions.h/cpp | code containing miscellaneous simple/supportive functions.

Moreover, "InstancesBPPMCF.zip" contains a txt-file for each of our test instances. These are spread over 5 folders, each corresponding to a different data set: "Dataset 1", "Dataset 2", "Dataset 3" and "Dataset 4" correspond to datasets D1, D2, D3 and D4 as introduced by Mehrani et al. (2022) (see https://github.com/saharnazmehrani/BPPMCF-IJOC), and "Triplets" corresponds to D5*. Instances D1*, D2*, D3* and D4* (introduced by Barkel et al. (2025) are obtained by setting the number of bins to the minimum number of required bins (see main.cpp), where the minimum number of required bins per instance is saved in the file "minNumberOfBinsPerInstance.txt".

Each instance file is structured as follows: For Datasets 1,2,3 and 4:

The first line is always a 1
The second line contains the number of bins (B)
The third line contains the bin capacity (W)
Then a blank line, followed by a 0-matrix of size BxW, followed by three more blank lines
The next line contains the number of colors (C)
The next line contains the number of items (I)
Then follows a blank line
Finally there is a matrix of size Ix2, where each row i gives the color c and the weight w of item i For Dataset 5 (Triplets):
The first line contains the number of items I
The remaining lines contain a matrix of size Ix2, where each row i gives the color c and the weight w of item i.

For questions on the code, send an e-mail to mathijsbarkel3@gmail.com.
