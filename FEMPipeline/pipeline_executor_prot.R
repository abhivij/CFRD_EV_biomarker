library(FEMPipeline)
source("FEMPipeline/dataset_pipeline_arguments_prot.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) > 1) {
  print(paste('Executing pipeline on dataset', args[2]))
  dparg = dataset_pipeline_arguments_prot[[strtoi(args[2])]]
  do.call(execute_pipeline, dparg)
} else {
  print('Executing pipeline on all datasets')
  for (dparg in dataset_pipeline_arguments_prot) {
    do.call(execute_pipeline, dparg)
  }
}