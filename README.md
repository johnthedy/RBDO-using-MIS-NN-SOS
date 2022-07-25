# RBDO-using-MIS-NN-SOS

This is RBDO (Reliability Based Design Optimization) code written in MATLAB language
MAIN.m is main code consist of SOS optimizer, random variable number and its properties can be edited in MAIN.m code
G.m consist of reliability failure function 
fobj1.m consist of SOS Metaheuristic optimization objective functionn 
fobj2.m consist of Multisphere Importance Sampling algorithm

Template reliability problem is short column problem subjected to axial load and moment can be seen in G.m consisting of 2 optimized parameter
and 6 random variable follow normal distribution for reliability analysis.

This recreate SOS optimization based on paper from:
Cheng M-Y, Prayogo D. Symbiotic Organisms Search: A new metaheuristic optimization algorithm. Computers & Structures. 2014;139:98-112.

M-IS reliability method is cited from:
Thedy J, Liao K-W. Multisphere-based importance sampling for structural reliability. Structural Safety. 2021;91:102099.

While Neural Network used in this code come from MATLAB toolbox:
