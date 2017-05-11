# functions for analysis of ATARiS2 models (post fit)
require(Matrix)
stopifnot(require(plyr))

coerce.na.to.zero <- function(x) {
    ifelse(!is.na(x), x, 0)
}

GenerateHairpinReportForHairpinSelection <- function(hairpins, m, target.gene = NULL, 
    r2.formula = "traditional") { 
    stopifnot(length(hairpins) > 0)

    obs.hairpin <- m$H[hairpins, , drop = F]
    gene.seeds2 <- m$MapHS_2[hairpins]  # the seed IDs for shRNAs in the gene solution of 'gene'
    seed2.solutions <- m$S[gene.seeds2, , drop = F]
    gene.seeds1 <- m$MapHS_1[hairpins]  # the seed IDs for shRNAs in the gene solution of 'gene'
    seed1.solutions <- m$S[gene.seeds1, , drop = F]
    
    gene.symbols <- foreach(hairpin = hairpins, .combine = c) %do% {
        gene.ids <- setdiff(m$MapHG[[hairpin]], target.gene)
        if (length(gene.ids) == 0) {
            ""
        } else {
            do.call(paste, as.list(m$names.genes[gene.ids]))
        }
    }
    
    mu <- foreach(hairpin = hairpins, .combine = rbind) %do% {
        m$mu[hairpin, m$MapSampleBatch]
    }
    obs.hat <- obs.hairpin - mu
    
    packed.gene.Hhat <- foreach(hairpin = hairpins, .combine = rbind) %do% {
        genes <- m$MapHG[[hairpin]]
        gene.solutions.sum <- rep(0, ncol(m$G))
        other.gene.solutions.sum <- rep(0, ncol(m$G))
        
        for (gene.index in seq(length(genes))) {
            gene.solutions <- coerce.na.to.zero(m$G[genes[gene.index], , drop = F])
            beta <- m$beta[[hairpin]][gene.index]
            b = coerce.na.to.zero(beta) * gene.solutions
            if (!is.null(target.gene) && genes[gene.index] == target.gene) {
                gene.solutions.sum <- gene.solutions.sum + b
            } else {
                other.gene.solutions.sum <- other.gene.solutions.sum + b
            }
            stopifnot(sum(is.na(gene.solutions.sum)) == 0)
            stopifnot(sum(is.na(other.gene.solutions.sum)) == 0)
        }
        rbind(c(1, gene.solutions.sum), c(0, other.gene.solutions.sum))
    }
    
    gene.Hhat <- packed.gene.Hhat[packed.gene.Hhat[, 1] == 1, -1, drop = F]
    other.gene.Hhat <- packed.gene.Hhat[packed.gene.Hhat[, 1] != 1, -1, drop = F]
    
    seed1.Hhat <- coerce.na.to.zero(m$alpha[hairpins] * seed1.solutions)
    seed2.Hhat <- coerce.na.to.zero(m$gamma[hairpins] * seed2.solutions)
    total.Hhat <- seed1.Hhat + seed2.Hhat + gene.Hhat + other.gene.Hhat
    
    stopifnot(sum(is.na(seed1.Hhat)) == 0)
    stopifnot(sum(is.na(seed2.Hhat)) == 0)
    stopifnot(sum(is.na(gene.Hhat)) == 0)
    stopifnot(sum(is.na(other.gene.Hhat)) == 0)
    
    names.gene.seeds1 <- m$names.seeds[gene.seeds1]
    
    r.squared.per.hp <- function(h.hat) {
        r <- sapply(seq(nrow(h.hat)), function(row) {
            pred <- h.hat[row, ]
            obs <- obs.hat[row, ]
            mask <- !(is.na(obs) | is.na(pred))
            caret::R2(pred[mask], obs[mask], formula = r2.formula, na.rm = T)
        })
        r
    }
    
    seeds.Rsquared <- r.squared.per.hp(seed1.Hhat + seed2.Hhat)
    gene.Rsquared <- r.squared.per.hp(gene.Hhat)
    other.gene.Rsquared <- r.squared.per.hp(other.gene.Hhat)
    other.gene.Rsquared[gene.symbols == ""] <- NA
    seed1.Rsquared <- r.squared.per.hp(seed1.Hhat)
    seed2.Rsquared <- r.squared.per.hp(seed2.Hhat)
    total.Rsquared <- r.squared.per.hp(total.Hhat)
    
    # make sure if we have an 'other.gene' we hava a non-na other.gene.Rsquared
    stopifnot(sum(is.na(other.gene.Rsquared) & gene.symbols != "") == 0)
    
    # make sure we always have a total.Rsquared
    stopifnot(sum(is.na(total.Rsquared)) == 0)
    
    # print(m$names.hps[hairpins]) print(gene.symbols) print(names.gene.seeds1)
    # print(m$names.hps[hairpins]) print('other') print(other.gene.Rsquared)
    
    df <- data.frame(shRNAID = m$names.hps[hairpins], gene.symbol = gene.symbols, 
        seed1.name = names.gene.seeds1, seed2.name = m$names.seeds[gene.seeds2], 
        other.gene.sol.R2 = other.gene.Rsquared, gene.sol.R2 = gene.Rsquared, seed1.sol.R2 = seed1.Rsquared, 
        seed2.sol.R2 = seed2.Rsquared, total.R2 = total.Rsquared, seeds.sol.R2 = seeds.Rsquared, 
        alpha = m$alpha[hairpins], gamma = m$gamma[hairpins], row.names = m$names.hps[hairpins], 
        stringsAsFactors = F)
    
    if (!is.null(target.gene)) {
        hps <- m$MapGH[[target.gene]]
        df$beta <- sapply(hps, function(hp) {
            i <- match(target.gene, m$MapHG[[hp]])
            m$beta[[hp]][i]
        })
        
        df$other.gene.symbol <- gene.symbols
        df$gene.symbol <- m$names.genes[target.gene]
    }
    
    df
}


# generates a data frame with a row for each shRNA of gene 'gene'.  The DF
# contains for each shRNA the % of variance explained by the gene solution, by
# the seed solution, and by both of them (alpha*seed + beta*gene).
GenerateHairpinReportForGene <- function(gene.name, m) {
    gene <- which(gene.name == m$names.genes)
    GenerateHairpinReportForHairpinSelection(m$MapGH[[gene]], m, gene)
}


# generate table with info for each shRNA
generateHairpinReport <- function(m, num.cores = 6, genes = NULL) {
    require(doMC)
    require(foreach)
    registerDoMC(num.cores)
    
    # use specified genes or all of them if not specified
    if (is.null(genes)) 
        genes <- m$names.genes
    
    gdf.rows <- lapply(genes, function(gene) {
        GenerateHairpinReportForGene(gene, m)
    })
    
    names(gdf.rows) <- NULL
    gdf <- do.call(rbind, gdf.rows)
    
    gdf
}

# generates a data frame with a row for each shRNA of seed 'seed.seq'.  The DF
# contains for each shRNA the % of variance explained by the gene solution, by
# the seed solution, and by both of them (alpha*seed + beta*gene).
GenerateHairpinReportForSeed <- function(seed.seq, m) {
    seed <- which(seed.seq == m$names.seeds)
    if (is.null(m$MapSH_1)) {
        hairpins <- unlist(m$MapSH_2[seed])
    } else {
        hairpins <- union(unlist(m$MapSH_2[seed]), unlist(m$MapSH_1[seed]))
    }
    GenerateHairpinReportForHairpinSelection(hairpins, m)
}

plotGeneSolutionRsquaredFromGdf <- function(label, gdf, show.legend = T) {
    par(mar = c(5, 13, 4, 2))
    mp <- barplot(t(as.matrix(gdf[, c("gene.sol.R2", "other.gene.sol.R2", "seed1.sol.R2", 
        "seed2.sol.R2")])), xlim = c(0, 1), beside = T, col = c("#E69F00", "#000000", 
        "#009E73", "#F0E442"), axes = TRUE, axisnames = FALSE, xlab = "% variance explained (R^2)", 
        main = paste("Amount of", label, "shRNA effects explained\n by gene/seed effects"), 
        horiz = TRUE)
    
    ylabels <- strtrim(rownames(gdf), 21)
    
    text(par("usr")[1], mp[2, ] + 1, labels = ylabels, srt = 0, adj = c(1.1, 1.1), 
        xpd = TRUE, cex = 0.9)
    
    if (show.legend) {
        legend("topright", legend = c("gene", "other genes", "seed1", "seed2"), fill = c("#E69F00", 
            "#000000", "#009E73", "#F0E442"))
    }
}

plotGeneSolutionRsquared <- function(gene, m, show.legend = T) {
    gdf <- GenerateHairpinReportForGene(gene, m)
    plotGeneSolutionRsquaredFromGdf(gene, gdf, show.legend)
}

plotSeedSolutionRsquaredFromGdf <- function(seed, gdf, show.legend) {
    par(mar = c(5, 13, 4, 2))
    mp <- barplot(t(as.matrix(gdf[, c("other.gene.sol.R2", "seed1.sol.R2", "seed2.sol.R2")])), 
        xlim = c(0, 1), beside = T, col = c("palegreen3", "palevioletred", "darkred"), 
        axes = TRUE, axisnames = FALSE, xlab = "% variance explained (R^2)", main = paste("Fraction of seed", 
            seed, "shRNA effects\n", "explained by gene/seed effects"), horiz = TRUE)
    ylabels <- paste(strtrim(rownames(gdf), 21), "\n", substring(rownames(gdf), 23))
    
    text(par("usr")[1], mp[2, ] + 1, labels = ylabels, srt = 0, adj = c(1.1, 1.1), 
        xpd = TRUE, cex = 0.9)
    
    if (show.legend) {
        legend("topright", legend = c("genes", "seed1", "seed2"), fill = c("palegreen3", 
            "palevioletred", "darkred"))
    }
}

plotSeedSolutionRsquared <- function(seed, m, show.legend = TRUE) {
    gdf <- GenerateHairpinReportForSeed(seed, m)
    plotSeedSolutionRsquaredFromGdf(seed, gdf, show.legend)
}


plotSeedSolutionRsquareForAllSeeds <- function(m) {
    for (s in m$names.seeds) {
        try(plotSeedSolutionRsquared(s, m))
    }
}


# plots x-y plot of two versions of a gene solution - ATARiS1 vs. ATARiS2 ATARiS1
# data is given with rownames being gene solution IDs (KRAS_1_ 101001) ATARiS2
# data has gene symbol as row names
plotGeneSolutions <- function(gene.symbol, GS1, GS2) {
    stopifnot(all(colnames(GS1) == colnames(GS2)))
    gs1.idx <- grep(paste("^", gene.symbol, "_1_", sep = ""), rownames(GS1))
    gs1 <- GS1[gs1.idx, ]
    gs2 <- GS2[gene.symbol, ]
    data <- data.frame(ATARiS1_solution = gs1, ATARiS2_solution = gs2)
    require(ggplot2)
    p <- qplot(ATARiS1_solution, ATARiS2_solution, data = data, main = paste("ATARiS 1/2 gene solutions for", 
        gene.symbol)) + theme(aspect.ratio = 1) + coord_cartesian(xlim = c(-4, 3), 
        ylim = c(-4, 3))
    p + geom_abline(slope = 1, intercept = 0)
}

gen_shRNA_info_table <- function(m) {
    tab <- data.frame(shRNA_ID = m$names.hps, target_gene = m$names.genes[m$MapHG], 
        seed_seq = m$names.seeds[m$MapHS_1], stringsAsFactors = F)
}


# scales the gene solutions (and betas) so that for each gene there's one shRNA
# with a beta equals 1.  this puts all the gene solutions on the same scale
# (original shRNA scale), with a gene solutions score representing the maximum
# gene effect achievable by any shRNA of the target gene.
scaleGeneSolutions <- function(m) {
    require(foreach)
    
    # construct mapping of gene -> hairpin, geneindex
    gene.mapping <- foreach(g = seq_along(m$MapGH)) %do% {
        shRNAs <- m$MapGH[[g]]
        # for each hairpin, what is the index of this gene?
        gene.index <- sapply(m$MapHG[shRNAs], function(genes) {
            match(g, genes)
        })
        
        data.frame(shRNA = shRNAs, gene.index = gene.index)
    }
    
    # for each gene solution
    for (g in seq_along(gene.mapping)) {
        ## get shRNAs targeting it if # hairpins > 1
        df <- gene.mapping[[g]]
        if (nrow(df) > 1) {
            ## copute max beta
            max.beta <- foreach(shRNA = df$shRNA, gene.index = df$gene.index, .combine = max) %do% 
                {
                  m$beta[[shRNA]][gene.index]
                }
            ## divide betas by max stopifnot(max.beta != 0 && !is.nan(max.beta))
            if (max.beta == 0 || is.nan(max.beta)) {
                cat("Gene", g, "had 0 or na max.beta\n")
                foreach(shRNA = df$shRNA, gene.index = df$gene.index, .combine = max) %do% 
                  {
                    m$beta[[shRNA]][gene.index] <- NA
                  }
            } else {
                foreach(shRNA = df$shRNA, gene.index = df$gene.index, .combine = max) %do% 
                  {
                    m$beta[[shRNA]][gene.index] <- m$beta[[shRNA]][gene.index]/max.beta
                  }
            }
            ## multiply G by max
            m$G[g, ] <- m$G[g, ] * max.beta
        }
    }
    m
    
}

exportScaledModel <- function(m, output.file) {
    colnames(m$G) <- m$names.samples
    colnames(m$S) <- m$names.samples
    rownames(m$G) <- m$names.genes
    rownames(m$S) <- m$names.seeds
    m <- scaleGeneSolutions(m)
    save(m, file = output.file)
    m
}

exportSolutions <- function(m, output.path = "./") {
    out.file <- file.path(output.path, "GeneSols.csv")
    write.csv(signif(m$G, 4), file = out.file, quote = F)
    
    out.file <- file.path(output.path, "SeedSols.csv")
    write.csv(signif(m$S, 4), file = out.file, quote = F)
}

generateHairpinReports <- function(m, output.path = "./", num.cores = 2, write.pdfs = T) {
    gdf <- generateHairpinReport(m, num.cores = num.cores)
    out.file <- file.path(output.path, "shRNA.explained.csv")
    write.csv(gdf, file = out.file, quote = F, row.names = T)
    
    if (write.pdfs) {
        out.file <- file.path(output.path, "gene.shRNA.explained.pdf")
        pdf(out.file, h = 6, w = 6)
        genes <- sort(m$names.genes)
        tmp <- lapply(genes, plotGeneSolutionRsquared, m)
        dev.off()
        
        out.file <- file.path(output.path, "seed.shRNA.explained.pdf")
        pdf(out.file, h = 6, w = 6)
        plotSeedSolutionRsquareForAllSeeds(m)
        dev.off()
    }
}

generateFitParamHistograms <- function(m, out.file) {
    pdf(out.file, h = 6, w = 6)
    hist(m$G, breaks = 100, main = sprintf("mean(G)=%f\nsd(G)=%f", mean(m$G, na.rm = T), 
        sd(m$G, na.rm = T)))
    hist(m$S, breaks = 100, main = sprintf("mean(S)=%f\nsd(S)=%f", mean(m$S, na.rm = T), 
        sd(m$S, na.rm = T)))
    hist(m$alpha, breaks = 100, main = sprintf("mean(alpha)=%f\nsd(alpha)=%f", mean(m$alpha, 
        na.rm = T), sd(m$alpha, na.rm = T)))
    hist(m$gamma, breaks = 100, main = sprintf("mean(gamma)=%f\nsd(gamma)=%f", mean(m$gamma, 
        na.rm = T), sd(m$gamma, na.rm = T)))
    beta <- unlist(m$beta)
    hist(beta, breaks = 100, main = sprintf("mean(beta)=%f\nsd(beta)=%f", mean(beta, 
        na.rm = T), sd(beta, na.rm = T)))
    dev.off()
}

write.expanded.gene.families <- function(a, expanded.file, families.file) {
    expanded.names <- strsplit(rownames(a), "&")
    expanded.rows <- lapply(seq_along(expanded.names), function(i) {
        rep(i, length(expanded.names[[i]]))
    })
    unlist(expanded.names)
    
    result <- a[unlist(expanded.rows), ]
    rownames(result) <- unlist(expanded.names)
    write.csv(result, expanded.file)
    
    families <- which(sapply(expanded.names, length) > 1)
    families.csv <- unlist(sapply(expanded.names[families], function(x) {
        paste(x, collapse = ",")
    }))
    writeLines(families.csv, families.file)
}

clean.gene.solutions <- function(s) {
    # drop all genes named 'NO_CURRENT_...'
    s <- as.matrix(s)
    bad <- grep("NO_CURRENT_", rownames(s))
    if (length(bad) > 0) {
        s <- s[-bad, ]
    }
    
    # sanity check gene names, specifically looking for gene names that have been
    # converted to dates
    cat("Check the following for bad gene names:\n")
    print(row.names(s)[grep("-", row.names(s))])
    print(row.names(s)[grep(" ", row.names(s))])
    print(row.names(s)[grep("\\.", row.names(s))])
    print(rownames(s)[grep("^[0-9].+", rownames(s))])
    
    # sanity check that 'X' has not been prepended to column names that start with a
    # number
    stopifnot(length(grep("^X[0-9].+", colnames(s))) == 0)
    
    # drop rows which are all NA
    s <- s[which(rowSums(is.na(s)) < ncol(s)), ]
    
    s
}

generateSecondBestDistribution <- function(x, out.file) {
    # x <-
    # read.csv('/Users/pmontgom/data/ach-2.4.3-demeter/null.shRNA.explained.csv',row.names=1,as.is=T)
    gene.sol.R2.by.gene <- split(x$gene.sol.R2, x$gene.symbol)
    second.best <- sapply(gene.sol.R2.by.gene, function(r2s) {
        c(sort(r2s, decr = T), NA, NA)[[2]]
    })
    print(quantile(second.best, 0.95, na.rm = T))
    write.table(second.best, file = out.file)
}

write.shRNA.without.seed.eff <- function(m, filename) {
    s.eff <- (m$S[m$MapHS_1, ] * m$alpha) + (m$S[m$MapHS_2, ] * m$gamma)
    residuals <- m$H - s.eff
    write.csv(residuals, file = filename)
}


