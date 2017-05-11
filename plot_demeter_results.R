load("results/m.solids.100iter.2seed.umap.Rdata")
m.2seed <- m
load("results/m.solids.100iter.1seed.umap.Rdata")
m.1seed <- m

# plot comparsion between 1 seed and 2 seed models
lonely.hairpins <- intersect(unlist(m.2seed$MapSH_1[m.2seed$lonely.solutions$S]),
                             unlist(m.2seed$MapSH_2[m.2seed$lonely.solutions$S]))

# number of hairpins which did not share seeds with anyone in m1
sum(m.1seed$lonely.solutions$S)
length(unlist(m.1seed$MapSH_2[m.1seed$lonely.solutions$S]))
# number of hairpins which did not share seeds with any other hairpins in m2
length(lonely.hairpins)

lonely.hairpins.mask<-rep(FALSE, length(m.2seed$alpha))
lonely.hairpins.mask[lonely.hairpins] <- TRUE

df <- data.frame(alpha=m.2seed$alpha[!lonely.hairpins.mask], gamma=m.2seed$gamma[!lonely.hairpins.mask])

require(ggplot2)
require(gridExtra)

hist_top <- ggplot(df)+geom_histogram(aes(x=alpha))

empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
       panel.background=element_blank(), 
       axis.text.x=element_blank(), axis.text.y=element_blank(),           
       axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- ggplot(df) + aes(x=alpha, y=gamma) + geom_point(alpha=0.05) + geom_smooth(method="lm", colour="red") 
hist_right <- ggplot(df)+geom_histogram(aes(x=gamma))+coord_flip()

grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

# plot comparision between alpha in old and new 
df <- data.frame(m.1seed.alpha=m.1seed$alpha, m.2seed.alpha=m.2seed$alpha)
ggplot(df) + aes(x=m.2seed.alpha, y=m.1seed.alpha) + geom_point() + geom_density2d() 

# Look at R^2 for hairpins
source("post_fit.R")
source("LGSmodel.R")
m.2seed.gdf<-generateHairpinReport(m.2seed)
m.1seed.gdf<-generateHairpinReport_old(m.1seed)

require(ggplot2)
# Generate R^2 hairpin comparison
ggplot(data.frame(total.R2.1seed=m.1seed.gdf$total.R2, total.R2.2seed=m.2seed.gdf$total.R2)) + aes(x=total.R2.1seed, y=total.R2.2seed) + geom_point() + geom_density2d() + geom_smooth(colour="red",method="lm")
ggplot(data.frame(total.R2.1seed=m.1seed.gdf$total.R2, total.R2.2seed=m.2seed.gdf$total.R2)) + aes(x=total.R2.1seed, y=total.R2.2seed) + geom_point(alpha=0.1) + geom_smooth(colour="red",method="lm")

require(reshape2)
ggplot(melt(m.2seed.gdf, id.vars=c("shRNAID"), measure.vars=c("gene.sol.R2", "seed1.sol.R2", "seed2.sol.R2", "total.R2"))) + aes(x=value) + geom_histogram() + facet_wrap(~variable, ncol=4)
ggplot(melt(m.1seed.gdf, id.vars=c("shRNAID"), measure.vars=c("gene.sol.R2", "seed.sol.R2", "total.R2"))) + aes(x=value) + geom_histogram() + facet_wrap(~variable, ncol=3)
help(melt)
names(m.2seed.gdf)

head(m.2seed.gdf$gene.seeds1)

hist(m.1seed.gdf$total.R2[m.1seed.gdf$total.R2>0.9])
hist(m.1seed.gdf$seed.sol.R2[m.1seed.gdf$total.R2>0.9])
hist(m.1seed.gdf$gene.sol.R2[m.1seed.gdf$total.R2>0.9])

# find points that got worse
worse_points <- m.1seed.gdf$total.R2>0.95 & m.2seed.gdf$total.R2 < 0.75
m1t<-table(m.1seed.gdf$gene.seed)
m2t<-table(union(m.2seed.gdf$gene.seed1, m.2seed.gdf$gene.seed2))
hist(m2t[m.1seed.gdf$gene.seed[worse_points]]  )

# write out gene solutions
dimnames(m.1seed$G)[[1]] <- m.1seed$names.genes
head(m.1seed$G)
write.table(m.1seed$G,"one-seed-gene-solutions.txt")
dimnames(m.2seed$G)[[1]] <- m.2seed$names.genes
head(m.2seed$G)
write.table(m.2seed$G,"two-seed-gene-solutions.txt")

plotGeneSolutionRsquared("SOX10", m.2seed)
write.csv(GenerateHairpinReportForGene("SOX10", m.2seed), file="~/dev/genedeps/sox10-shrna-quality.csv")

# plot a variety of comparisions
all(m.1seed$names.hps == m.2seed$names.hps)
df <- data.frame(hairpin = m.1seed$names.hps,
                 alpha.1seed = m.1seed$alpha,
                 alpha.2seed = m.2seed$alpha, 
                 beta.1seed = m.1seed$beta,
                 beta.2seed = m.2seed$beta,
                 gene.R2.1seed = m.1seed.gdf$gene.sol.R2,
                 gene.R2.2seed = m.2seed.gdf$gene.sol.R2,
                 seed1.R2.1seed = m.1seed.gdf$seed.sol.R2,
                 seed1.R2.2seed = m.2seed.gdf$seed1.sol.R2,
                 seed2.R2.2seed = m.2seed.gdf$seed2.sol.R2)
require(reshape2)
#dfm <- melt(df, "hairpin")
ggplot(df) + aes(x=alpha.2seed, y=alpha.1seed) + geom_point(alpha=0.1)
ggplot(df) + aes(x=beta.2seed, y=beta.1seed) + geom_point(alpha=0.1)
ggplot(df) + aes(x=gene.R2.2seed, y=gene.R2.1seed) + geom_point(alpha=0.1)
ggplot(df) + aes(x=seed1.R2.2seed, y=seed1.R2.1seed) + geom_point(alpha=0.1)
ggplot(df) + aes(x=seed1.R2.2seed, y=seed2.R2.2seed) + geom_point(alpha=0.1) + geom_density2d(colour="red")
hist(atan(df$seed1.R2.2seed/df$seed2.R2.2seed))

#add the distribution of Pearson correlation between 1seed and 2seed gene solutions
hist(sapply(seq_along(m.1seed$names.genes), function (i) { cor(m.1seed$G[i,], m.2seed$G[i,])  }), xlab='pearson correlation')
cor.per.gene <- sapply(seq_along(m.1seed$names.genes), function (i) { cor(m.1seed$G[i,], m.2seed$G[i,])  })
find.second.best <- function(gdf) {
  sapply(seq_along(m.1seed$names.genes), function (i) { 
    hairpins <- unlist(m.1seed$MapGH[i])
    geneR2 <- gdf$gene.sol.R2[hairpins]
    geneR2[order(-geneR2)[2]]
  })
}
second.best.R2.per.gene.1seed <- find.second.best(m.1seed.gdf)
second.best.R2.per.gene.2seed <- find.second.best(m.2seed.gdf)

par(mfrow=c(1,4))
hist(cor.per.gene)
hist(cor.per.gene[second.best.R2.per.gene.2seed>0.25])
hist(cor.per.gene[second.best.R2.per.gene.2seed>0.50])
hist(cor.per.gene[second.best.R2.per.gene.2seed>0.75])

plot(ecdf(second.best.R2.per.gene.1seed), xlab="R^2 of second best hairpin for gene", ylab="fraction of genes", main="ecdf")
plot(ecdf(second.best.R2.per.gene.2seed), col="red", add=T)
legend(0.6,0.3, # places a legend at the appropriate place 
       c("1-seed model","2-seed model"), # puts text in the legend 
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),col=c("black","red")) # gives the legend lines the correct color and width

# investigate top left corner of R^2s
df$hairpin.count.with.seed.1seed <- sapply(m.1seed$MapSH[m.1seed$MapHS_2[df$hairpin]], length)
df$seed1 <- m.2seed$MapHS_1[df$hairpin]
df$seed2 <- m.2seed$MapHS_2[df$hairpin]
df$gene <- m.2seed$MapHG[df$hairpin]
df$hairpin.count.with.gene.1seed <- sapply(m.1seed$MapGH[m.1seed$MapHG[df$hairpin]], length)
ggplot(df) + aes(x=gene.R2.2seed, y=gene.R2.1seed, colour=seed1.R2.1seed) + geom_point()

ggplot(df[df$gene.R2.1seed > 0.9 & df$gene.R2.2seed < 0.1,]) + aes(x=gene.R2.2seed, y=gene.R2.1seed) + geom_point()

funny <- df[df$gene.R2.1seed > 0.975 & df$gene.R2.2seed < 0.025,]
funny$gene <- m.2seed$names.genes[m.2seed$MapHG[funny$hairpin]]
funny$s1 <- m.2seed$names.seed[m.2seed$MapHS_1[funny$hairpin]]
funny$s2 <- m.2seed$names.seed[m.2seed$MapHS_2[funny$hairpin]]


summary(funny$gene)
sum(table(funny$s1) > 1)
sum(table(funny$s2) > 1)
t<-table(funny$s2)
t[t>1]
head(funny)

summary(df[df$gene.R2.1seed > 0.9 & df$gene.R2.2seed < 0.025, ])

what frac


#in how many genes there are hairpins sharing the same seed (seed 1, seed 2, both)

g<-m.2seed$MapHG
s1<-m.2seed$MapHS_1
s2<-m.2seed$MapHS_2
ms1<-min(s1, s2)
ms2<-max(s1, s2)
t<-table(c(paste(g, s1), paste(g, s2)))
head(paste(g, s1, s2))
sum(t[t > 1])

generateOutputsFromModel(m.2seed, num.cores=1)

funny[1,]
plotGeneSolutionRsquared("RNASEH2A", m.1seed, F)
plotGeneSolutionRsquared("RNASEH2A", m.2seed, F)

plotSeedSolutionRsquared("CCTCAAT", m.1seed, F)
plotSeedSolutionRsquared("CCTCAAT", m.2seed, F)


source("post_fit.R")
source("LGSmodel.R")
plotSeedSolutionRsquared("CCTCAAT", m.1seed)
GenerateHairpinReportForSeed("CGTTCTC", m.1seed)
output.path<-"./"
out.file <- file.path(output.path, 'seed.shRNA.explained.pdf')
pdf(out.file, h=6, w=6)
plotSeedSolutionRsquareForAllSeeds(m.2seed)
dev.off()


# create data frame with both sets of scores
m.2seed.gdf<-generateHairpinReport(m.2seed)[m.1seed$names.hps,]
m.1seed.gdf<-generateHairpinReport(m.1seed)[m.2seed$names.hps,]
all(m.1seed$names.hps == m.2seed$names.hps)
df <- data.frame(hairpin = m.1seed$names.hps,
                 alpha.1seed = m.1seed$alpha,
                 alpha.2seed = m.2seed$alpha, 
                 beta.1seed = m.1seed$beta,
                 beta.2seed = m.2seed$beta,
                 gene.R2.1seed = m.1seed.gdf$gene.sol.R2,
                 gene.R2.2seed = m.2seed.gdf$gene.sol.R2,
                 seed1.R2.1seed = m.1seed.gdf$seed2.sol.R2,
                 seed1.R2.2seed = m.2seed.gdf$seed1.sol.R2,
                 seed2.R2.2seed = m.2seed.gdf$seed2.sol.R2,
                 seed1.lonely.1seed = m.1seed$lonely.solutions$S[m.1seed$MapHS_2],
                 seed1.lonely.2seed = m.2seed$lonely.solutions$S[m.2seed$MapHS_1],
                 seed2.lonely.2seed = m.2seed$lonely.solutions$S[m.2seed$MapHS_2]
)
ggplot(df) + aes(y=gene.R2.1seed, x=gene.R2.2seed, colour=seed1.lonely.1seed) + geom_point(alpha=0.5)

ggplot(df[df$seed1.lonely.2seed & df$seed2.lonely.2seed,]) + aes(y=gene.R2.1seed, x=gene.R2.2seed) + geom_point(alpha=0.5)

all(df$seed1.lonely.2seed == df$seed1.lonely.2seed)

funny <- df[df$gene.R2.1seed > 0.975 & df$gene.R2.2seed < 0.025,]
head(funny)

plotGeneSolutionRsquared("COX17", m.1seed, FALSE)
plotGeneSolutionRsquared("COX17", m.2seed, FALSE)

plotSeedSolutionRsquared()

m.1seed.gdf["ACATCCTACTTCCTCAATGAA_RNASEH2A",]

# cor of genes which have one lonely seed hairpin
# first divide genes into two classes based on whether they have a lonely hairpin or not
genes.with.lonely.seed <- unique(m.1seed$names.genes[m.1seed$MapHG[m.1seed$lonely.solutions$S[m.1seed$MapHS_2]]])
gene.with.seed.class <- data.frame(gene=m.1seed$names.genes, has.lonely=m.1seed$names.genes %in% genes.with.lonely.seed)
row.names(gene.with.seed.class) <- gene.with.seed.class$gene
all(m.1seed.gdf$gene.symbol == m.2seed.gdf$gene.symbol)

second.best.R2.per.gene.1seed <- sapply(seq_along(m.1seed$names.genes), function (i) {
  hairpins <- unlist(m.1seed$MapGH[i])
  geneR2 <- m.1seed.gdf$gene.sol.R2[hairpins]
  geneR2[order(-geneR2)[2]]
})

hist(sapply(m.1seed$names.genes, function(g) { 
  h<-unlist(m.1seed$MapGH[g]);
  cor(m.1seed.gdf$beta[h], m.2seed.gdf$beta[h])
} 
))


plotBetaCorForGenes <- function(genes) {
  hist(sapply(genes, function(g) { 
    h<-unlist(m.1seed$MapGH[g]);
    cor(m.1seed$beta[h], m.2seed.gdf$beta[h])
  } 
  ))
}

match(c(4,5), 4)

plotBetaCorForGenes <- function(genes, main) {
  genes.index <- match(genes, m.1seed$names.genes)
  hist(sapply(genes.index, function(g) { 
    cor(m.1seed$G[g,], m.2seed$G[g,])
  } 
  ), main=main, xlab="correlation of gene solutions between models")
}


par(mfcol=c(2,1))
plotBetaCorForGenes(gene.with.seed.class$gene[gene.with.seed.class$has.lonely], 'For genes with a lonely seed')
plotBetaCorForGenes(gene.with.seed.class$gene[!gene.with.seed.class$has.lonely], 'For genes without any lonely seed')

head(second.best.R2.per.gene.1seed)

