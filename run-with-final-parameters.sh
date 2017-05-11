mkdir -p processed/final
./eval-demeter-final --full.data.file processed/reformatted.Rdata --dest.dir processed/final --learning.rate=0.005 --randseed=1 --G.S=4e-5 --alpha.beta.gamma=0.8 --max.num.iter=110 > processed/final/stdout.txt 2> processed/final/stderr.txt

