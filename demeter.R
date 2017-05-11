library(foreach)

#######
# Linear Gene-Seed model
#######
# We model each shRNA data across multiple samples as a linear combination of
# gene effects and seed effects.
# I.e., hp_i ~= mu_i + alpha_i * Seed_(i) + beta_i * Gene_(i),
# with Seed_(i) being the seed effects (across all samples) that correspond to
# the seed sequence of hp_i. Same for Gene_(i).
# mu_i is a hp-specific constant effect
# alpha_i and beta_i define the mixture of seed and gene effects respectively.
# 
# To represent in matrix notation, first define N_samples, N_genes, N_seeds,
# N_hps to represent the number of samples, genes, seeds and shRNAs in the
# model, repsectively. Then we define:
# H (N_hps x N_samples) - (observed) shRNA readouts,
# S (N_seeds x N_samples) - (parameter) seed effects,
# G (N_genes x N_samples) - (parameter) gene effects,
# mu (N_hps x 1) - (parameter) constant shRNA effects,
# beta (N_hps) - (parameter) gene effects co-efficients, a vector per hairpin.  length(beta[h]) == length(MapHG[h])
# alpha (N_hps x 1) - (parameter) 11-17nt seed effects co-efficients,
# gamma (N_hps x 1) - (parameter) 12-18nt seed effects co-efficients,
# MapHG (N_hps) - hairpin to vector of index of target gene solution for each shRNA.
# MapGH (named list) N_genes entries - vector of shRNA indices in each gene entry.
# MapHS_1 (N_hps x 1) - index of target seed solution (seed region 11-17 nt) for each shRNA.
# MapSH_1 (named list) N_seeds entries - vector of shRNA indices in each seed (seed region 11-17 nt) entry.
# MapHS_2 (N_hps x 1) - index of target seed solution (seed region 12-18 nt) for each shRNA.
# MapSH_2 (named list) N_seeds entries - vector of shRNA indices in each seed (seed region 12-18 nt) entry.
# H_hat (N_hps x N_samples) - (inferred) inferred shRNAs readouts (H matrix)
# lonely.solutions (list) - logical vectors indicating a solution being "lonely" (= associated
#     with only one shRNA)
#  * G (N_genes, vector)
#  * S (N_seeds, vector)
# hp - hyper parameters (list)
#   lambda[1], ..., lambda[6] - lambda hyper parameters
#
# The model can then be expressed as:
#   H_hat = mu*1 + diag(alpha) * S[MapHS_1, ] + diag(beta) * G[MapHG, ] + diag(gamma) * S[MapHS_2, ]
# with diag(a) being a matrix (N_hps x N_hps) with elements of a (N_hps x 1) in
# its diagonal and all others zeros, and 1 being a vector (1 x N_hps) of 1's.

# Loss function:
# 
# L(S,G,mu,alpha,beta,gamma) = sum_ij(|H_ij - H_hat_ij|^2) 
#     + lambda_1*S^2
#     + lambda_2*G^2
#     + lambda_3*alpha^2 
#     + lambda_4*beta^2
#     + lambda_5*mu^2 
#     + lambda_6*gamma^2

# library(Rcpp)
# sourceCpp("SGD.cpp")
library(futile.logger)

# Constructs an empty model object (the structure of the model)
# Needs to be initialized before used (InitModelParams)
# no.gene.solutions is a logical matrix, where TRUE means we have prior knowledge that 
# there is no gene effect in this cell line for the gene.
ConstructModel <- function(N_hps, N_samples, N_seeds, N_genes, MapHS_1, MapSH_1, MapHS_2, MapSH_2, MapHG, MapGH, MapSampleBatch,
                           .names.hps, .names.samples, .names.genes, .names.seeds, no.gene.solutions) {
  model <- list()

#  model$lambda <- lambda # c(0.01, 0.01, 0.01, 0.01, 0) # regularization hyper parameters
  model$N_hps <- N_hps
  model$N_samples <- N_samples
  model$N_seeds <- N_seeds
  model$N_genes <- N_genes
  
  model$MapHS_1 <- MapHS_1
  model$MapSH_1 <- MapSH_1
  model$MapHS_2 <- MapHS_2
  model$MapSH_2 <- MapSH_2
  model$MapHG <- MapHG
  model$MapGH <- MapGH
  model$MapSampleBatch <- MapSampleBatch
  model$lonely.solutions <- list()
  
  model <- within(model, {
    names.hps <- .names.hps
    names.samples <- .names.samples
    names.genes <- .names.genes
    names.seeds <- .names.seeds
  })

  if(is.null(no.gene.solutions)) {
    empty.mask <- matrix(FALSE, ncol=length(.names.samples), nrow=length(.names.genes))
    model$shouldBeZeroInG <- empty.mask
  } else {
    model$shouldBeZeroInG <- create.fullsized.mask(.names.genes, .names.samples, t(no.gene.solutions))
  }  
  model
}


# initializes the model params
# shRNA.data.wide is the observed matrix of scrren readouts as data frame of gct
# available methods: "zeros", "random", "means"
InitModelParams <- function(model, shRNA.data.wide,
                            method="random") {
  flog.debug("InitModelParams: Called with method='%s'", method)

  shRNA.data <- as.matrix(shRNA.data.wide)
  
  stopifnot(class(shRNA.data) == "matrix")
  library(foreach)
  
  model <- within(model, {
    batches = length(unique(MapSampleBatch))
    stopifnot(ncol(shRNA.data) == N_samples && nrow(shRNA.data) == N_hps)
    H <- shRNA.data
    N_datapoints <- sum(!is.na(H))
    
    if (method=="zeros") {
      S <- matrix(0, nrow=N_seeds, ncol=N_samples)
      G <- matrix(0, nrow=N_genes, ncol=N_samples)
      
      mu <- matrix(0, nrow=N_hps, ncol=batches)
      alpha <- matrix(0, nrow=N_hps, ncol=1)
      beta <- matrix(0, nrow=N_hps, ncol=1)
      gamma <- matrix(0, nrow=N_hps, ncol=1)
    }
    if (method=="random") {
      S <- matrix(rnorm(N_seeds*N_samples), nrow=N_seeds, ncol=N_samples)
      G <- matrix(rnorm(N_genes*N_samples), nrow=N_genes, ncol=N_samples)
      
      mu <- matrix(rnorm(N_hps), nrow=N_hps, ncol=batches)
      alpha <- matrix(rnorm(N_hps), nrow=N_hps, ncol=1)
      beta <- matrix(rnorm(N_hps), nrow=N_hps, ncol=1)      
      gamma <- matrix(rnorm(N_hps), nrow=N_hps, ncol=1)
    }
    if (method=="means") {
      # below line sums for each seed seq all the hps associated with it and 
      # divides by their number.
      # TODO: use scale() instead

      flog.debug("Calculating mu")
      mu <- (foreach(batch=seq(batches), .combine=cbind) %do% {
      	matrix(rowMeans(H[,MapSampleBatch == batch], na.rm=T), nrow=N_hps, ncol=1) } )
      mu[is.nan(mu)] <- 0

      flog.debug("Calculating Hstar")
      Hstar = H - (mu[,MapSampleBatch])
      flog.debug("Calculating S")
      S <- matrix(foreach(seed = names.seeds, .combine=rbind) %do% {
        seed.hps <- union(MapSH_1[[seed]], MapSH_2[[seed]])
        colMeans(Hstar[seed.hps, , drop=F], na.rm=T)
        }, nrow=length(names.seeds))
      flog.debug("Calculating G")
      G <- foreach(gene.hps = MapGH, .combine=rbind) %do% {
        colMeans(Hstar[gene.hps, , drop=F], na.rm=T)
      }
      rownames(G) <- names.genes
      flog.debug("Done with G and S")
      rm(Hstar)

      # zero out elements of G and S that are NAs due to NAs in H
      S[is.na(S)] <- 0
      G[is.na(G)] <- 0
      
      flog.debug("Creating alpha")
      alpha <- matrix(0.5, nrow=N_hps, ncol=1)
      flog.debug("Creating Beta")
      beta <- lapply(MapHG, function(genes) { rep(0.5, length(genes)) })
      flog.debug("Creating gamma")
      gamma <- matrix(0.5, nrow=N_hps, ncol=1)
      
      stopifnot(sum(is.na(mu)) == 0)
    }
    
  })
    
  flog.debug("InitModelParams: Done.")
  
  model
}

mask.G.and.S <- function(model) {
  model <- InitLonelySolutions(model)

  # NA out elements in G for those elements in no.gene.solutions
  flog.info("NAs before masking %d", sum(is.na(model$G)))
  if(!is.null(model$no.gene.solutions)) {
    model$G <- masked.update.mat(model$G, t(model$no.gene.solutions), NA)
  }
  flog.info("NAs after masking %d", sum(is.na(model$G)))
  flog.info("zeros after masking %d", sum(model$G==0, na.rm=T))

  # initialize elements to zero when we have no expression
  # zero out elements in G that we expect to be zero
  model$G[model$shouldBeZeroInG] <- 0

  return(model)
}

# computes and returns H_hat - the estimated H matrix, based on all parameters/hidden vars.
# H hat shouldn't contain any NAs
# H_hat = mu*1 + diag(alpha) * S[MapHS_1, ] + diag(beta) * G[MapHG, ]
# with diag(a) being a matrix (N_hps x N_hps) with elements of a (N_hps x 1) in
# its diagonal and all others zeros, and 1 being a vector (1 x N_samples) of 1's.
ComputeHHat <- function(model) {
  
  hhat <- with(model, {
    S.with.zero.contrib <- S
    S.with.zero.contrib[is.na(S.with.zero.contrib)] <- 0
    seed_1 <- as.vector(alpha) * S.with.zero.contrib[MapHS_1, , drop=F]
    seed_2 <- as.vector(gamma) * S.with.zero.contrib[MapHS_2, , drop=F]
    gene <- CalculateGeneEffectC(MapHG, beta, G)
    tmu <- mu[,MapSampleBatch]
    seed_1 + seed_2 + gene + tmu
  }) 
  
  dimnames(hhat) <- dimnames(model$H)
  hhat
}


# computes (H-Hhat)^2
TotalSquaredError <- function(model, hhat=NULL, training.points) {
  if (is.null(hhat)) {
    hhat <- ComputeHHat(model)
  }
  
  # make sure we have a prediction for everywhere we have an observation
  #stopifnot(all( (is.na(hhat) | is.na(model$H) ) == is.na(model$H) ))

  sum((model$H-hhat)^2, na.rm=T)
}

# Finds the best alpha and beta values for each shRNA, given the current
# associated gene and seed solutions. Fitting is done using simple linear
# regression. This function is used to initialize the alpha's and beta's before
# fitting of all parameters.
FitAlphaBeta <- function(m) {
    
    coeffs <- foreach(h = 1:m$N_hps) %do% {
        if ((h%%10000) == 1 || h < 10) {
            flog.debug("FitAlphaBeta h=%d", h)
        }
        
        seed1_eff = t(m$S[m$MapHS_1[h], , drop = F])
        seed2_eff = t(m$S[m$MapHS_2[h], , drop = F])
        gene_eff = t(m$G[m$MapHG[[h]], , drop = F])
        
        seed1_eff[is.na(seed1_eff)] <- 0
        seed2_eff[is.na(seed2_eff)] <- 0
        gene_eff[is.na(gene_eff)] <- 0
        
        data = as.data.frame(cbind(m$H[h, ], seed1_eff, seed2_eff, gene_eff))
        names(data) <- c("h", "seed1_eff", "seed2_eff", paste("g", seq(ncol(gene_eff)), 
            sep = ""))
        
        if (all(is.na(data[, 1]))) {
            # if we have NAs for all H, then that's fishy...
            flog.info("hairpin %s had no measurements, so skipping lm fit, and setting alpha, beta, gamma, mu to 0", h)
            list(mu = 0, alpha = 0, gamma = 0, beta = rep(0, ncol(gene_eff)))
        } else if (all(is.na(data[, -1]))) {
            # if we have NAs for all G and S for this hairpin
            cat("hairpin %s had all NAs, so skipping lm fit, and setting alpha, beta, gamma, mu to 0", h)
            list(mu = 0, alpha = 0, gamma = 0, beta = rep(0, ncol(gene_eff)))
        } else {
            res <- tryCatch({
                lm(h ~ ., data = data, na.action = na.exclude, model = F, qr = F)
            }, error = function(e) {
                flog.error("error: %s", e)
            })
            
            # if seed1=seed2, then seed2_eff is a free variable and populated with NA
            res$coeff[is.na(res$coeff)] <- 0
            
            list(mu = res$coeff[[1]], alpha = res$coeff[[2]], gamma = res$coeff[[3]], 
                beta = as.numeric(res$coeff[4:length(res$coeff)]))
        }
    }
    
    m$mu[, 1] <- sapply(coeffs, function(x) {
        x$mu
    })
    m$mu[is.nan(m$mu)] <- 0
    m$alpha <- matrix(sapply(coeffs, function(x) {
        x$alpha
    }), ncol = 1)
    m$gamma <- matrix(sapply(coeffs, function(x) {
        x$gamma
    }), ncol = 1)
    m$beta <- lapply(coeffs, function(x) {
        x$beta
    })
    
    # constrain alpha/gamma/beta by > 0 so that we start our search inside the
    # feasible region
    m$alpha[, 1] <- ifelse(m$alpha[, 1] < 0, 0, m$alpha[, 1])
    m$gamma[, 1] <- ifelse(m$gamma[, 1] < 0, 0, m$gamma[, 1])
    m$beta <- lapply(m$beta, function(x) {
        ifelse(x < 0, 0, x)
    })
    
    m
}

# Mask out elements in G and S which only have a <=1 non-NA observation.
InitLonelySolutions <- function(m) {
    m <- within(m, {
        flog.debug("Setting lonely solution G")
        
        # browser() create a matrix of genes x sample, where each element is true if
        # there are > 1 non-NA elements in H for that gene in that cell line
        lonely.gene.samples <- foreach(hps = MapGH, .combine = rbind) %do% {
            apply(!is.na(H[hps, , drop = F]), 2, sum) <= 1
        }
        
        G[lonely.gene.samples] <- NA
        
        flog.debug("Setting lonely solution S")
        lonely.seed.samples <- foreach(seed = names.seeds, .combine = rbind) %do% 
            {
                hps <- union(MapSH_1[[seed]], MapSH_2[[seed]])
                apply(!is.na(H[hps, , drop = F]), 2, sum) <= 1
            }
        
        S[lonely.seed.samples] <- NA
    })
    
    flog.debug("finished lonely solutions")
    m
}


# returns the r-squared between each column of m and v
computeRSquaredMatVec <- function(m, v) {
    stopifnot(nrow(m) == length(v))
    
    sapply(1:ncol(m), function(i) {
        if (all(is.na(m[, i])) || all(is.na(v))) {
            NA
        } else {
            summary(lm(m[, i] ~ v, na.action = na.exclude))$r.squared
        }
    })
}


# returns a vector of R^2 for column-column pairs of two input matrics.  I.e.,
# cor(a[,1], b[,1])^2, cor(a[,2], b[,2])^2, ...
computeRSquareMatMat <- function(m1, m2) {
    stopifnot(all(dim(m1) == dim(m2)))
    
    sapply(1:ncol(m1), function(i) {
        if (all(is.na(m1[, i])) || all(is.na(m2[, i]))) {
            NA
        } else {
            summary(lm(m1[, i] ~ m2[, i], na.action = na.exclude))$r.squared
        }
    })
}


debug.msg <- function(s) {
    # require(pryr) bytes <- as.numeric(mem_used())
    bytes <- 0
    flog.debug("[mem %.2f GB] %s", bytes/(1024 * 1024 * 1024), s)
}

create.fullsized.mask <- function(.rownames, .colnames, mask) {
    full.mask <- matrix(FALSE, ncol = length(.colnames), nrow = length(.rownames))
    rownames(full.mask) <- .rownames
    colnames(full.mask) <- .colnames
    shared.rows <- intersect(rownames(mask), .rownames)
    shared.cols <- intersect(colnames(mask), .colnames)
    full.mask[shared.rows, shared.cols] <- mask[shared.rows, shared.cols]
    return(full.mask)
}

# shRNA.data - in long format: columns: shRNA_ID, Description, cell_line, batch,
# value lambda - vector of 5 regularization term coefficients
prepareLGSModel <- function(shRNA.data.wide, MapHairpinGene, MapSampleBatch, init.method = "means", 
    num.cores = 1, shuffle = FALSE, no.gene.solutions = matrix(FALSE, 0, 0)) {
    library(futile.logger)
    stopifnot(!is.null(MapHairpinGene))
    stopifnot(!is.null(MapSampleBatch))
    stopifnot(sum(is.na(MapSampleBatch)) == 0)
    stopifnot(is.matrix(shRNA.data.wide))
    stopifnot(is.list(MapHairpinGene))
    
    debug.msg("Entered fitLGSModel()")
    
    # only shuffle after we've collapsed our data down to one row per hairpin
    if (shuffle) {
        debug.msg("Shuffling gene/shRNA labels")
        rownames(shRNA.data.wide) <- sample(rownames(shRNA.data.wide))
    }
    
    debug.msg("Creating model object")
    m <- CreateLGSModelFromGCTDataFrame(shRNA.data.wide, MapSampleBatch, MapHairpinGene, 
        no.gene.solutions)
    
    debug.msg("Initializing model...")
    m <- InitModelParams(m, shRNA.data.wide, method = init.method)
    debug.msg("Done. (Initializing model)")
    
    return(m)
}


is.tail.strictly.increasing <- function(series, samples = 4) {
    # looking for a consistent increase over n samples
    tail.samples <- tail(series, samples)
    stopifnot(!is.null(tail.samples))
    return(all(!is.nan(tail.samples)) && length(tail.samples) == samples && all(order(tail.samples) == 
        seq(samples)))
}

optimizeLGSModel <- function(m, lambda, learning.rate, max.num.iter, training.datapoints, 
    performance.log.file = NULL) {
    flog.info("optimizeLGSModel with lambda: %s", paste(lambda, collapse = " "))
    m$lambda <- lambda
    
    require(foreach)
    library(Rcpp)
    sourceCpp("SGD.cpp")
    
    performance <- printTotalLoss(m)
    performance$iteration <- NA
    performance.log <- list(performance)
    tryCatch({
        
        tic()
        performance <- printTotalLoss(m, training.datapoints)
        performance$iteration <- 0
        performance.log[[length(performance.log) + 1]] <- performance
        last_loss <- performance$total
        
        flog.info("Completed loss calculation (%.1fs)", toc(print = F))
        
        ### datapoints <- which(!is.na(m$H)) validation.set.size <- length(datapoints) *
        ### 0.05 validation.set.mask <- sample(c(rep(TRUE, validation.set.size), rep(FALSE,
        ### length(datapoints)-validation.set.size))) training.datapoints <-
        ### datapoints[!validation.set.mask] validation.datapoints <-
        ### datapoints[validation.set.mask]
        
        last.checkpoint.time <- Sys.time()
        
        dynamic.learning.rate <- learning.rate
        set.seed(1)
        randomized.datapoints <- sample(training.datapoints)
        
        tic()
        # learning.rate <- find.best.learning.rate(learning.rate, m,
        # randomized.datapoints[seq(length(randomized.datapoints)*0.5)])
        flog.info("Picked learning rate: %f (%.1fs)", learning.rate, toc(print = F))
        # set decay such that learning rate is halved after 50 iterations
        learning.rate.decay <- (learning.rate * 2 - 1)/(50 * learning.rate)
        
        for (iter in 1:max.num.iter) {
            flog.info("Starting iteration #%d/%d", iter, max.num.iter)
            
            # use each datapoint to update parameters in a random order
            if (iter > 1) {
                randomized.datapoints <- sample(training.datapoints)
            }
            
            # if(iter > 10) { dynamic.learning.rate <- learning.rate / (1 + learning.rate *
            # learning.rate.decay * iter) } else {
            dynamic.learning.rate <- learning.rate
            # }
            flog.info("New learning rate: %f (%.1fs)", dynamic.learning.rate, toc(print = F))
            
            tic()
            updated <- LearnStochasticIterationC(m, randomized.datapoints, list(learning.rate = rep(dynamic.learning.rate, 
                6)))
            m$mu <- updated$mu
            m$alpha <- updated$alpha
            m$beta <- updated$beta
            m$gamma <- updated$gamma
            m$S <- updated$S
            m$G <- updated$G
            
            flog.info("Completed iteration #%d (%.1fs)", iter, toc(print = F))
            
            tic()
            performance <- printTotalLoss(m, training.datapoints)
            performance$iteration <- iter
            cur_loss <- performance$total
            performance.log[[length(performance.log) + 1]] <- performance
            flog.info("Completed loss calculation (%.1fs)", toc(print = F))
            
            if (is.tail.strictly.increasing(sapply(performance.log, function(x) {
                x$VSE
            }), 6)) {
                flog.info("Last 6 iterations were increasing VSE.  Stopping.")
                break
            }
            
            # save a checkpoint every two minutes now <- Sys.time() if(now -
            # last.checkpoint.time > 2*60) { save(m, file='m.checkpoint.Rdata')
            # last.checkpoint.time <- Sys.time() }
        }
        
        colnames(m$G) <- m$names.samples
        rownames(m$G) <- m$names.genes
        
        m
    }, finally = {
        performance.log <- do.call(rbind, performance.log)
        save(performance.log, file = performance.log.file)
    })
}

printTotalLoss <- function(m, training.points = NULL, skip.print = F) {
    library(futile.logger)
    
    if (is.null(training.points)) {
        training.points <- which(!is.na(m$H))
        validation.points <- NULL
    } else {
        training.mask <- matrix(FALSE, nrow = nrow(m$H), ncol = ncol(m$H))
        training.mask[training.points] <- TRUE
        validation.points <- which(!training.mask)
        validation.pt.count <- length(validation.points)
    }
    
    train.error <- calculateTotalErrorC(m, training.points)
    # print(train.error)
    
    TSE <- round(train.error[2]/train.error[1], 5)
    RET <- round(train.error[3]/train.error[1], 5)
    CPE <- round(train.error[4]/train.error[1], 5)
    VSE <- NA
    
    stopifnot(!is.nan(TSE))
    stopifnot(!is.nan(RET))
    if (is.null(validation.points)) {
        flog.info("Loss = %f (MSE) + %f (MRT) + %f (CPE) = %f (Total)", TSE, RET, 
            CPE, TSE + RET + CPE)
    } else {
        validation.error <- calculateTotalErrorC(m, validation.points)
        # print('val') print(validation.error)
        VSE <- round(validation.error[2]/validation.error[1], 5)
        if (!skip.print) {
            flog.info("Loss = %f (MSE) + %f (MRT) + %f (CPE) = %f (Total), VSE: %f", 
                TSE, RET, CPE, TSE + RET + CPE, VSE)
        }
    }
    return(data.frame(TSE = TSE, RET = RET, CPE = CPE, VSE = VSE, total = TSE + RET + 
        CPE))
}

CreateMapGH <- function(MapHG, names.genes) {
    flog.debug("Populating HairpinGeneRelationship")
    
    hgr.count <- sum(sapply(MapHG, length))
    
    hp <- rep(0, hgr.count)
    gene <- rep(0, hgr.count)
    
    next.pos <- 1
    for (i in seq_along(MapHG)) {
        genes <- MapHG[[i]]
        indices <- next.pos:(length(genes) + next.pos - 1)
        gene[indices] <- names.genes[genes]
        hp[indices] <- i
        next.pos <- next.pos + length(genes)
    }
    
    flog.debug("Populating HairpinGeneRelationshipGrouped")
    HairpinGeneRelationshipGrouped <- split(hp, gene)
    MapGH <- HairpinGeneRelationshipGrouped[names.genes]
}

# processes a DF generated from a gct file of shRNA data and creates an empty LGS
# model from it df should have the following format: * rownames: shRNA ID -
# sequence * additional columns: data.  * colnames: cell line names a named list
# of hairpin name -> gene name
CreateLGSModelFromGCTDataFrame <- function(df, MapSampleBatch, MapHairpinGenes, no.gene.solutions) {
    stopifnot(!is.null(rownames(df)))
    stopifnot(!is.null(colnames(df)))
    stopifnot(is.list(MapHairpinGenes))
    
    flog.debug("CreateLGSModelFromGCTDataFrame")
    names.hps <- rownames(df)
    N_hps <- nrow(df)
    names.samples <- colnames(df)
    N_samples <- ncol(df)
    names.genes <- foreach(genes = MapHairpinGenes, .combine = union) %do% {
        genes
    }
    N_genes <- length(names.genes)  # unique genes
    
    gene.to.index <- setNames(seq(length(names.genes)), names.genes)
    
    seed_seqs_2_8 <- substring(names.hps, 11, 17)
    seed_seqs_1_7 <- substring(names.hps, 12, 18)
    names.seeds <- unique(union(seed_seqs_2_8, seed_seqs_1_7))
    N_seeds <- length(names.seeds)
    
    flog.debug("Populating MapHG")
    MapHG <- lapply(names.hps, function(hp) {
        gene.to.index[MapHairpinGenes[[hp]]]
    })
    names(MapHG) <- names.hps
    
    MapGH <- CreateMapGH(MapHG, names.genes)
    
    flog.debug("Populating MapHS_1")
    MapHS_1 <- match(seed_seqs_2_8, names.seeds)
    names(MapHS_1) <- rownames(df)
    
    MapHS_2 <- match(seed_seqs_1_7, names.seeds)
    names(MapHS_2) <- rownames(df)
    
    # MapSH <- lapply(names.seeds, function(s) which(names.seeds[MapHS] == s))
    # equivalent to above, much much faster:
    grouped.MapSH_1 <- split(seq_along(MapHS_1), names.seeds[MapHS_1])
    grouped.MapSH_2 <- split(seq_along(MapHS_2), names.seeds[MapHS_2])
    
    MapSH_1 <- sapply(names.seeds, function(x) {
        c()
    })
    MapSH_1[match(names(grouped.MapSH_1), names.seeds)] = grouped.MapSH_1
    MapSH_2 <- sapply(names.seeds, function(x) {
        c()
    })
    MapSH_2[match(names(grouped.MapSH_2), names.seeds)] = grouped.MapSH_2
    
    stopifnot(length(MapSH_1) == length(MapSH_2))
    
    
    flog.debug("calling ConstructModel")
    # build a model
    m <- ConstructModel(N_hps, N_samples, N_seeds, N_genes, MapHS_1, MapSH_1, MapHS_2, 
        MapSH_2, MapHG, MapGH, MapSampleBatch, names.hps, names.samples, names.genes, 
        names.seeds, no.gene.solutions)
    
    return(m)
}


tic <- function(gcFirst = TRUE, type = c("elapsed", "user.self", "sys.self")) {
    type <- match.arg(type)
    assign(".type", type, envir = baseenv())
    if (gcFirst) 
        gc(FALSE)
    tic <- proc.time()[type]
    assign(".tic", tic, envir = baseenv())
    invisible(tic)
}

# print - should elapsed time be printed
toc <- function(print = TRUE) {
    type <- get(".type", envir = baseenv())
    toc <- proc.time()[type]
    tic <- get(".tic", envir = baseenv())
    if (print) {
        print(toc - tic)
    }
    invisible(toc - tic)
}


# print numbers formatted as percentage
percent <- function(x, digits = 2, format = "f", ...) {
    paste(formatC(100 * x, format = format, digits = digits, ...), "%", sep = "")
}

replace.bad.gene.names <- function(orig.genes) {
  genes <- orig.genes

  bad <- c("Mar", "Dec", "Sep")
  corrected <- c("MARCH", "DEC", "SEPT")
  for(i in seq(length(bad))) {
    genes <- sub(paste0("0?([0-9]+)-", bad[[i]]), paste0(corrected[[i]], "\\1"), genes)
  }

  genes <- sub("(^|_)MAR([0-9]+)$", paste0("\\1","MARCH", "\\2"), genes)

  mismatches <- which(orig.genes != genes)
  if(length(mismatches) > 0) {
    cat("Remapping names which look mangled by excel:\n")
    for(i in mismatches) {
      cat("\t", orig.genes[[i]], "->", genes[[i]], "\n")
    }
  }

  genes
}

assert.no.bad.gene.names <- function(genes) {
    bad.genes <- grep("[0-9]+-(Mar|Sep|Dec)", genes, v = T)
    bad.genes <- union(grep("^MAR[0-9]$", genes, v = T), bad.genes)
    if(length(bad.genes) > 0) {
      print("bad gene names")
      print(bad.genes)
      stopifnot(length(bad.genes) == 0)
    }
}

assert.no.bad.sample.names <- function(samples) {
    bad.samples <- grep("char", samples, v = T)
#    print(bad.samples)
    stopifnot(length(bad.samples) == 0)
}

# Given a GCT as a dataframe, return a matrix which has samples for column names,
# hairpin as the row names as well as a list mapping hairpins -> gene
collapse.to.unique.hairpins <- function(gct) {
    flog.info("starting collapse.to.unique.hairpins")
    
    hairpins <- substr(rownames(gct), 1, regexpr("_", rownames(gct)) - 1)
    genes <- substr(rownames(gct), regexpr("_", rownames(gct))+1, nchar(rownames(gct)))
    genes <- replace.bad.gene.names(genes)
    assert.no.bad.gene.names(genes)
    
    unique.hairpins <- unique(hairpins)
    rows.to.keep <- match(unique.hairpins, hairpins)

    # pull out one row per hairpin, drop the 'Description' column and convert it to a
    # matrix
    data <- gct[rows.to.keep, ,drop=F]
    rownames(data) <- unique.hairpins
    
    assert.no.bad.sample.names(colnames(data))
    
    flog.info("done (starting collapse.to.unique.hairpins)")
    list(data = data, hg = data.frame(hairpins = hairpins, genes = genes, stringsAsFactors = F))
}

combine.long.and.gene.map <- function(x, y) {
    m <- merge(x$data, y$data, all = T, by = 0)
    data <- as.matrix(m[, -1])
    row.names(data) <- m[, 1]
    stopifnot(!is.null(colnames(data)))
    stopifnot(!is.null(rownames(data)))
    
    list(data = data, hg = unique(rbind(x$hg, y$hg)), MapSampleBatch = c(x$MapSampleBatch, 
        y$MapSampleBatch))
}

drop.promiscuous.hairpins <- function(collapsed, max.targeted) {
    flog.info("starting drop.promiscuous.hairpins")
    grouped <- split(collapsed$hg$genes, collapsed$hg$hairpins)
    grouped <- sapply(grouped, length)
    bad.hairpins <- names(grouped[grouped > max.targeted])
    flog.info("done (drop.promiscuous.hairpins) dropped %d hairpins targeting > %d genes", 
        length(bad.hairpins), max.targeted)
    
    list(data = collapsed$data[!(rownames(collapsed$data) %in% bad.hairpins), , drop = F], 
        hg = collapsed$hg[!(collapsed$hg$hairpins %in% bad.hairpins), ])
}

cbind.by.rownames <- function(a, b) {
    m1 <- cbind(a, b[match(rownames(a), rownames(b)), ,drop=F])
    missing.in.a <- setdiff(rownames(b), rownames(a))
    m2 <- cbind(matrix(NA, ncol=ncol(a), nrow=length(missing.in.a), dimnames=list(missing.in.a, colnames(a))), b[missing.in.a,,drop=F])
    rbind(m1, m2)
}

preprocess.fc.matrices <- function(fc.matrices, max.genes.targeted = 10) {
    # need to merge gene -> hairpin mappings
    
    combine.long.and.gene.map <- function(x, y) {
        data <- cbind.by.rownames(x$data, y$data)
        stopifnot(!is.null(colnames(data)))
        stopifnot(!is.null(rownames(data)))
        
        list(data = data, hg = unique(rbind(x$hg, y$hg)), MapSampleBatch = c(x$MapSampleBatch, 
            y$MapSampleBatch))
    }
    
    all.batches <- (foreach(i = seq_along(fc.matrices), .combine = combine.long.and.gene.map) %do% 
        {
            gct <- fc.matrices[[i]]
            collapsed <- collapse.to.unique.hairpins(gct)
            
            flog.info("rows before dropping promiscuous hairpins: %s", nrow(collapsed$hg$hairpins))
            collapsed <- drop.promiscuous.hairpins(collapsed, max.genes.targeted)
            flog.info("rows after dropping promiscuous hairpins: %s", nrow(collapsed$hg$hairpins))

            collapsed$MapSampleBatch <- setNames(rep(i, ncol(collapsed$data)), colnames(collapsed$data))
         
            stopifnot(!is.null(colnames(collapsed$data)))
            stopifnot(!is.null(rownames(collapsed$data)))
            
            collapsed
        })

    assert_that(!is.null(all.batches))
    
    all.batches$hg <- collapse.identical.genes(all.batches$hg$genes, all.batches$hg$hairpins)
    
    all.batches
}

merge.data.for.DEMETER <- function(fc.matrices, batch.per.sample) {
    assert_that(is.list(fc.matrices))
    # check fc.matrices
    for(fc.matrix in fc.matrices) {
        assert_that(is.matrix(fc.matrix))
        assert_that(!is.null(rownames(fc.matrix)))
        assert_that(!is.null(colnames(fc.matrix)))
    }
    
    assert_that(is.numeric(batch.per.sample))
    assert_that(!is.null(names(batch.per.sample)))

    processed <- preprocess.fc.matrices(fc.matrices)
    processed$MapSampleBatch <- batch.per.sample[names(processed$MapSampleBatch)]

    MapHairpinGene <- split(processed$hg$gene, processed$hg$hairpin)

    m <- prepareLGSModel(processed$data, MapHairpinGene, processed$MapSampleBatch, no.gene.solutions = NULL)
    m
}

DEMETER <- function(m, learning.rate, lambda, max.num.iter, dest.dir) {
    dir.create(dest.dir, recursive=T)
    m <- mask.G.and.S(m)
    flog.info("Fitting alphas and betas...")
    tic()
    m <- FitAlphaBeta(m)
    flog.info(sprintf("Done. (%.1fs)", toc(print = F)))
    
    flog.info("Optimizing model\n")
    non.na.points <- which(!is.na(m$H))
    training.datapoints <- non.na.points
    m <- optimizeLGSModel(m, lambda = lambda, learning.rate = learning.rate, max.num.iter = max.num.iter, 
        training.datapoints = training.datapoints, performance.log.file = paste0(dest.dir, "/perf.log"))

    write.shRNA.without.seed.eff(m, paste0(dest.dir, "/shRNA-without-seed.csv"))
    m <- exportScaledModel(m, paste0(dest.dir, "/m.betascaled.Rdata"))
    generateFitParamHistograms(m, paste0(dest.dir, "/param_histograms.pdf"))
    exportSolutions(m, dest.dir)
    generateHairpinReports(m, dest.dir)
    
    g <- read.csv(paste0(dest.dir, "/GeneSols.csv"), row.names = 1, check.names = F)
    g <- clean.gene.solutions(g)
    write.csv(g, paste0(dest.dir, "/GeneSolsCleaned.csv"))
    
    # clean up seed solutions
    s <- read.csv(paste0(dest.dir, "/SeedSols.csv"), row.names = 1, check.names = F)
    write.csv(s, paste0(dest.dir, "/SeedSolsCleaned.csv"))

    # z-score gene solutions
    zg <- (g - mean(g, na.rm=T)) / sd(g, na.rm=T)

    write.expanded.gene.families(zg, paste0(dest.dir, "/ExpandedGeneZSolsCleaned.csv"), 
        paste0(dest.dir, "/GeneFamilies.csv"))
}

