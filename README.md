ddRAD-seq-Pipeline
==================

Created by Michael Sorenson and Jeffrey DaCosta

Copyright (c) 2011-2013 Boston University. All rights reserved.

This pipeline was created to process “double-digest” restriction-site associated DNA sequences
(ddRAD-seq) as described in DaCosta and Sorenson (PLoS One, submitted). The scripts are written for
execution in the python (version 3) language, and in some cases are dependent on other programs or
specification files. Some scripts run jobs in parallel using the multiprocessing module and the ssh
command to access nodes on a computer cluster. These scripts were optimized to run on our Diskless
Remote Boot Linux (DRBL) cluster and may need editing to run on other clusters. The software is free
and distributed WITHOUT warranty; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.


See MANUAL.txt for details on pipeline steps and scripts.
