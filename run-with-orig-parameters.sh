mkdir -p processed/orig-hp
./eval-demeter-final --full.data.file processed/reformatted.Rdata --dest.dir processed/orig-hp --learning.rate=0.005 --randseed=1 --G.S=4e-5 --alpha.beta.gamma=0.9 --max.num.iter=250 > processed/orig-hp/stdout.txt 2> processed/orig-hp/stderr.txt

