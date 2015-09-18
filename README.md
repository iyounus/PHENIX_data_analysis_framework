# [PHENIX](http://www.phenix.bnl.gov/) Data Analysis Framework

This analysis framework was created at University of New Mexico (UNM) to be used for physics analyses carried out by the PHENIX research group at UNM.

This code uses [ROOT](https://root.cern.ch/) classes extensively. The output data is in [ROOT Trees](https://root.cern.ch/doc/master/classTTree.html). ROOT Trees are commonly used as by physicists to store large experimental data.

This code allows to extract relevant data, needed by UNM research group, from PHENIX experiment data at RCF ([RHIC Computing Facility](https://www.racf.bnl.gov/)). The output data files are much smaller and can be easily transfered to local server at UNM.

## Description of Mudules

IContainers module provides the objects stored in the trees. IPdstReco is the main module that allows to read the data files at RCF, extract the data with some constraints, and store the output trees. This module cannot be compiled on local machine because it requires PHENIX environment.

The rest of the modules are used for different type of analyses. The IExample module provides a simple use case.
