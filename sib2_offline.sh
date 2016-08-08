#!/bin/bash

#
# Script for run SIB2 (offline) - see 'sib2_offine.par' file for ordering variables
#

#echo "1960111601 1961123024 '/home/nelsonvn/folders/SiB2_HRR/data2' 3600.0 20 3 6 3 0 1 \
#	0 0.85 -45.0 -20.5 0.458 -0.2 3.5E-06 7.797 0.08 1.0 298.0 298.0 297.0 \
#	0.0 0.0 0.0 0.0 0.75 0.85 0.98 0.0 0.0 0.0001 20.0 0.9999 298.0 298.0 297.0 \
#	0.70 0.85 0.98 45.0 45.0 1.0" | time -p ./sib2_offline.out

echo " 1995010101 2008123123 '/home/roilan/Dropbox/Dissertacao/FROM_SIB2DIAG.txt' \
3600 20 5 2 3 0 1 0 0.85 -47.63 -21.61 0.379898 -0.4713136 2.2223463e-05 8.0804005 \
0.04185 1.0 298.0 298.0 297.0 0.0 0.0 0.0 0.0 0.7997442 0.8543309 0.9572984 0.0 \
6.0 0.0001 20.0 0.9999 298.0 298.0 297.0 0.7997442 0.8543309 0.9572984 45.0 45.0 \
1.0 0.98 999 '/home/roilan/Dropbox/Dissertacao/' " | time -p ./sib2_offline.out 

