#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-4 pipeline_executor_tra.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 5-8 pipeline_executor_tra.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 9-12 pipeline_executor_tra.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-16 pipeline_executor_tra.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 17-26 pipeline_executor_tra.pbs
qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 27-30 pipeline_executor_tra.pbs