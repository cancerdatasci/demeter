
fetch.and.preprocess <- function(urls.to.fetch) {
    require(foreach)
    
    # need to merge gene -> hairpin mappings
    
    all.batches <- (foreach(i = seq_along(urls.to.fetch), .combine = combine.long.and.gene.map) %do% 
        {
            flog.info("reading %s", urls.to.fetch[i])
            # if(file.exists(urls.to.fetch[i])) { read the file as a gct file
            data <- read.table(urls.to.fetch[i], skip = 2, head = T, row.names = 1, 
                check.names = F, quote = NULL)
            # load(urls.to.fetch[i]); } else { load(url(urls.to.fetch[i])); }
            gct <- as.data.frame(data)
            flog.info("done reading")
            
            collapsed <- collapse.to.unique.hairpins(gct)
            
            flog.info("rows before dropping promiscuous hairpins: %s", nrow(collapsed$hg$hairpins))
            save(collapsed, file = paste("collapsed-", i, ".Rdata", sep = ""))
            collapsed <- drop.promiscuous.hairpins(collapsed, 10)
            flog.info("rows after dropping promiscuous hairpins: %s", nrow(collapsed$hg$hairpins))
            
            collapsed$MapSampleBatch <- setNames(rep(i, ncol(collapsed$data)), colnames(collapsed$data))
            
            # list(data=collapsed$data, MapGeneToHairpins=collapsed$MapGeneToHairpins,
            # MapSampleBatch=MapSampleBatch)
            stopifnot(!is.null(colnames(collapsed$data)))
            stopifnot(!is.null(rownames(collapsed$data)))
            
            collapsed
        })
    
    all.batches$hg <- collapse.identical.genes(all.batches$hg$genes, all.batches$hg$hairpins)
    
    all.batches
}

collapse.identical.genes <- function(genes, hairpins) {
    flog.info("start collapse.identical.genes")
    
    require(data.table)
    require(foreach)
    dt <- data.table(hairpin = hairpins, gene = genes)
    per.gene <- dt[, list(hairpins = do.call(paste, as.list(sort(hairpin)))), by = gene]
    gene.to.row <- data.table(gene = dt$gene, row.index = seq(nrow(dt)))
    setkey(gene.to.row, gene)
    mkgenelist <- function(genes) {
        genes <- as.list(genes)
        genes[["sep"]] = "&"
        do.call(paste, genes)
    }
    identical.genes <- per.gene[, list(first.gene = gene[[1]], genes = mkgenelist(gene)), 
        by = hairpins]
    per.gene.list <- split(hairpins, genes)
    new.table <- foreach(first.gene = identical.genes$first.gene, gene = identical.genes$genes, 
        .combine = rbind) %do% {
        data.frame(gene = gene, hairpin = per.gene.list[[first.gene]], stringsAsFactors = F)
    }
    
    flog.info("finished collapse.identical.genes")
    
    new.table
}

