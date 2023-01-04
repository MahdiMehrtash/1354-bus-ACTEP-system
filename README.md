# ACTEP 1354-bus Test System
This repository includes a new large-scale test system for Transmission Expansion Planning with AC Networks Model (ACTEP) studies. The proposed 1354-bus ACTEP system is based on Case1354pegase available in MATPOWER. The system is modified to be useful for ACTEP studies, and the candidate lines and candidate generating units are added to the system based on their role in improving the operation of the system. Global-TEP (a global solver for the ACTEP problem) is used to find the ACTEP global solutions with a guaranteed optimality gap.
To better understand the proposed ACTEP 1354-bus test system, you can refer to the manuscript: M. Mehrtash, B. F. Hobbs, and Y. Cao, "A Large-Scale Test System for Transmission Expansion Planning with AC Networks Model," 2022 IEEE Texas Power and Energy Conference (TPEC), doi: 10.1109/TPEC54980.2022.9750848.

## Developers
This test system has been developed by:
    
    Mahdi Mehrtash, Johns Hopkins University
    Benjamin F. Hobbs, Johns Hopkins University
    Yankai Cao, University of British Columbia
    
## Citing Global-TEP
If you find ACTEP 1354-bus Test System useful for your work, you might cite the manuscript:

    @articlemehrtash2022,
    author={Mehrtash, Mahdi and Hobbs, Benjamin F and Cao, Yankai,
    title="A Large-Scale Test System for Transmission Expansion Planning with AC Networks Model",
    conference="2022 IEEE Texas Power and Energy Conference (TPEC)",
    year="2022",
    month="March",
    day="1",
    doi="10.1109/TPEC54980.2022.9750848",
    url="[https://ieeexplore.ieee.org/document/9445630](https://ieeexplore.ieee.org/abstract/document/9750848)"
    }

## Features
    Test system data are in .CSV format
    The ACTEP code is ompatible with Julia 1.6 and the lastest version of JuMP
