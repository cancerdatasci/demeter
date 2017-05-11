a.samples <- paste0("SA", seq(10))
b.samples <- paste0("SB", seq(10))

rand.dna <- function(len, count) {
  sapply(seq(count), function(x) { paste0(sample(c("A", "T", "G", "C"), size=len, replace = T), collapse="") })
}

seed.seqs <- rand.dna(7,4)
hps <- paste0(rand.dna(10, 20), sample(seed.seqs, 10, replace=T), rand.dna(4, 20))

gene.hps <- paste0(hps, "_G", as.character(sample(seq(4), replace=T, length(hps))), sep="")
# make a single gene family.  Some code breaks if we don't have at least one

a.guides <- gene.hps
b.guides <- gene.hps

a <- matrix(runif(length(a.samples) * length(a.guides)), nrow=length(a.guides), ncol=length(a.samples), dimnames=list(a.guides, a.samples))
b <- matrix(runif(length(b.samples) * length(b.guides)), nrow=length(b.guides), ncol=length(b.samples), dimnames=list(b.guides, b.samples))

# make a single gene family.  Some code breaks if we don't have at least one
rows.to.dup <- a[grep(".*_G1$", rownames(a)), , drop=F]
stopifnot(nrow(rows.to.dup) > 0)
rownames(rows.to.dup) <- sub("_G1", "_GDUP", rownames(rows.to.dup))

a <- rbind(a, rows.to.dup)

print(a)

batches <- setNames( c(rep(1, length(a.samples)), rep(2, length(b.samples))), c(a.samples, b.samples) )
learning.rate <- 0.005
lambda <- c(4e-5, 4e-5, 0.9, 0.9, 0.0, 0.9)
max.num.iter <- 5

#library(futile.logger)
#library(foreach)

library(assertthat)
source("dependencies.R")
source("demeter.R")
source("post_fit.R")

m <- merge.data.for.DEMETER(list(a, b), batches)
dir.create("test/results")
DEMETER(m, learning.rate, lambda, max.num.iter, "test/results")

