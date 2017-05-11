kubeque sub -r memory=6G -u @kubeque/files --params processed/hyperparameter-samples-retry.csv \
    /usr/bin/time -v ./eval-demeter --holdout.fold='{fold}' --config.index='{config.index}' \
      --fold.count=5 --randseed=1 --G.S='{G.S}' --alpha.beta.gamma='{alpha.beta.gamma}' \
      --full.data.file '^processed/reformatted.Rdata' --output.file=out.rds \
      --learning.rate='{learning.rate}' --max.num.iter=350
      
