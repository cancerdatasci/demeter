DEMETER v2.20.2
------------

Scripts for deconvolving seed effects from gene effects in shRNA screens via the DEMETER method.

The contents are as follows:

README                 - this file

kubeque/*              - Files use to run hyperparameter sampling in parallel via kubeque

setup-r                - Script for installing all R dependencies

sample-hyperparams     - Randomly sample from space of hyperparameters and generate a table of combinations

merge-eval-outputs     - Produces a summary table from multiple runs of eval-demeter (as part of hyperparameter search)

preprocess-gcts        - Given multiple GCT files with fold change data and a file containing mapping from celline to batch, 
                         preprocess the GCT files to produce a file which eval-demeter can use as input.

eval-demeter           - Script for running DEMETER given a preprocessed file, and fit the model given the data and hyperparameters

*.R, *.cpp             - These files are the supporting R code and the stochastic gradient descent 
                         that was implemented in cpp which is used by the scripts above.

The process of running DEMETER is in three steps:

1. Preprocess the data from raw GCT files from Achilles into an R structure that DEMETER operates on

example: ./preprocess-gcts test/reformatted.Rdata test/testdata/batches.csv test/testdata/fc_matrix_1.gct,test/testdata/fc_matrix_2.gct

2. Find decent hyperparameters by performing a random sampling of the space and picking the best

example:
  # preprocess gct files 
  ./preprocess-gcts processed/reformatted.Rdata data/TableS1_SampleInfo_excluding_2lines.csv data/Achilles_v2.4.6.rnai.gct,data/Achilles_v2.19.2_mapped.rnai.gct

  # prepare for hyperparameter sampling
  ./sample-hyperparams --output.file=processed/hyperparameter-samples.csv
 
  # kick off batch job
  kubeque sub --fetch processed/hpsearch -u @kubeque/files \
      --params processed/hyperparameter-samples.csv \
      ./eval-demeter --holdout.fold='{fold}' --config.index='{config.index}' \
        --fold.count=5 --randseed=1 --G.S='{G.S}' --alpha.beta.gamma='{alpha.beta.gamma}' \
        --full.data.file '^processed/reformatted.Rdata' --output.file=out.rds \
        --learning.rate='{learning.rate}' --max.num.iter=500

  # collect results
  ./merge-eval-outputs processed/hpsearch

3. Run DEMETER with the hyperparameters chosen based on which lead to the best VSE (validation squared error)

  ./eval-demeter --full.data.file '^processed/reformatted.Rdata' 
      --output.file=processed/final.rds \
      --learning.rate=0.005 --randseed=1 \
      --G.S=4e-5 --alpha.beta.gamma=0.8
