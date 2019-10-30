#!/bin/bash

python2.7 ../scripts/surviv.py MDS-MPN-U_inc.txt MDS-MPN-U_surv.txt MDS-MPN-U_results.txt
python2.7 ../scripts/surviv.py MDS-MPN-RS-T_inc.txt MDS-MPN-RS-T_surv.txt MDS-MPN-RS-T_results.txt
python2.7 ../scripts/surviv.py CMML_inc.txt CMML_surv.txt CMML_results.txt
python2.7 ../scripts/surviv.py MDS_inc.txt MDS_surv.txt MDS_results.txt
python2.7 ../scripts/surviv.py AML_inc.txt AML_surv.txt AML_results.txt
