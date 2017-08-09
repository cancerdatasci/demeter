DEMETER v2.20.2
------------

Scripts for deconvolving seed effects from gene effects in shRNA screens via the DEMETER method.

The contents are as follows:

README                 - this file

setup-r                - Script for installing all R dependencies

sample-hyperparams     - Randomly sample from space of hyperparameters and generate a table of combinations

merge-eval-outputs     - Produces a summary table from multiple runs of eval-demeter (as part of hyperparameter search)

preprocess-gcts        - Given multiple GCT files with fold change data and a file containing mapping from celline to batch, 
                         preprocess the GCT files to produce a file which eval-demeter can use as input.

eval-demeter           - Script for running DEMETER given a preprocessed 
                         file, and fit the model given the data and hyperparameters

*.R, *.cpp             - These files are the supporting R code and the stochastic gradient descent 
                         that was implemented in cpp which is used by the scripts above.

test/testdata          - The test/testdata directory contains some sample files
                         which do not represent any real data, but only
                         are examples of the required format for the input
                         data.  These files are used in the example steps
                         below.  The advantage of using these test files when
                         walking through the example is that these files are
                         tiny and hence each step quickly executes.

The process of running DEMETER is in three steps:

1. Preprocess the data from raw GCT files from Achilles into an R structure that DEMETER operates on

example: 
  # preprocess gct files
  ./preprocess-gcts test/reformatted.Rdata test/testdata/batches.csv test/testdata/fc_matrix_1.gct,test/testdata/fc_matrix_2.gct

The batches.csv file must have at least two columns in it: "Name" which
should be the cell line name as it appears in the GCT file and 
"DEMETER batch" which should be a number (the first batch should be denoted
as 1) which represents which batch the line belongs to.

2. Find decent hyperparameters by performing a random sampling of the space and picking the best

example:
  # sample hyperparameters to evaluate
  ./sample-hyperparams --output.file=test/hyperparameter-samples.csv
 
  # Generate commands to run eval-demeter on each hyperparameter.  
  mkdir test/hpsearch
  ./run-eval-for-each test/reformatted.Rdata test/hyperparameter-samples.csv \
     test/hpsearch/out > hyperparams-eval.sh
  # Evaluate each hyperparameter setting sequentially.  The
  # shell script (hyperparams-eval.sh) can be run sequentially,
  # probably will take a prohibitively long time.  However, each line can be
  # run in parallel, so each line could be run via submitting to a batch
  # engine such as SGE, slurm, PBS, etc to execute in parallel and get the
  # results back in a more timely manner.
  sh ./hyperparams-eval.sh

  # collect results and print the hyperparameters with the best performance
  ./merge-eval-outputs test/hpsearch test/hpsearch-merged.rds

3. Run DEMETER with the hyperparameters (G.S, alpha.beta.gamma,
   max.num.iter -- the values below are those chosen based running 
   the hyperparameter search on the achilles 2.19.2 and achilles 2.4.6 
   datasets) chosen based on which lead to the best VSE (validation 
   squared error)

example:
  ./eval-demeter-final --full.data.file=test/reformatted.Rdata \
      --dest.dir=test/final --learning.rate=0.005 --randseed=1 \
      --G.S=4e-5 --alpha.beta.gamma=0.9 --max.num.iter=10

-------------------------------------------------------------------------

The above were walking through the steps using the tiny (artificial) test 
data in test/testdata.  However, the steps are basically the same to
reproduce the DEMETER 2.20.2 dataset based on the real Achilles data.

1. Make a "data" directory and download the v2.4.6 and v2.19.2 datasets to that
directory.   Also, download supplemental table 1, which contains the
"SampleInfo" sheet, and save as a csv named "data/sampleinfo.csv"

2. Run the preprocessing:

  ./preprocess-gcts data/reformatted.Rdata data/sampleinfo.csv data/Achilles_v2.4.6.rnai.gct,data/achilles-v2-19-2-mapped-rnai_v1-data.gct

3. Run the process of selecting hyperparameters

  ./sample-hyperparams --output.file=data/hyperparameter-samples.csv
 
  # Generate commands to run eval-demeter on each hyperparameter.  
  mkdir data/hpsearch
  ./run-eval-for-each data/reformatted.Rdata data/hyperparameter-samples.csv \
     data/hpsearch/out > data/hyperparams-eval.sh
  sh data/hyperparams-eval.sh
  ./merge-eval-outputs data/hpsearch data/hpsearch-merged.rds

4. After reviewing performance of various hyperparameters, we chose G.S=4e-5 alpha.beta.gamma=0.9 max.num.iter=250

  ./eval-demeter-final --full.data.file=data/reformatted.Rdata \
      --dest.dir=data/final --learning.rate=0.005 --randseed=1 \
      --G.S=4e-5 --alpha.beta.gamma=0.9 --max.num.iter=250

once complete, the file data/final/ExpandedGeneZSolsCleaned.csv will contain
the final z-scored gene solutions reported as v2.20.2.

