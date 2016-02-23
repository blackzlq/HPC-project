#!/bin/sh
#BSUB -J Liqin-project
#BSUB -o out_file
#BSUB -e error_file
#BSUB -n 1
#BSUB -q ht-10g


./test 
