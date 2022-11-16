#!/bin/bash

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 1-4 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 5-8 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 9-12 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 13-16 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 17-20 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 21-24 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 25-28 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 29-32 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 33-36 pipeline_executor.pbs

qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 37-40 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 41-44 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 45-48 pipeline_executor.pbs




qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 49-58 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 59-64 pipeline_executor.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 65-68 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 69-72 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 73-76 pipeline_executor.pbs

qsub -l select=1:ncpus=4:mem=124gb,walltime=4:00:00 -M a.vijayan@unsw.edu.au -m ae -J 77-82 pipeline_executor.pbs


qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 83-86 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 87-90 pipeline_executor.pbs
qsub -l select=1:ncpus=16:mem=124gb,walltime=12:00:00 -M a.vijayan@unsw.edu.au -m ae -J 91-94 pipeline_executor.pbs
