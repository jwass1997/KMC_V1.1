#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --mem=1gb
#SBATCH --export=NONE

module load compiler/gnu/14.1
module load lib/boost/1.80.0