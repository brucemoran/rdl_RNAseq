#! R

##load libraries
libs <- c("sleuth","tidyverse","ggrepel","biomaRt","pheatmap","RColorBrewer","PoiClaClu","DESeq2","rhdf5","ggplot2","devtools","cowplot", "reshape2", "limma", "edgeR", "DESeq2")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

##input parameters
OUTDIR <- "kallisto-sleuth_full"
RDATADIR <- "kallisto-sleuth_full/RData"
dir.create(RDATADIR, showWarnings=FALSE, recursive=TRUE)

CLIN <- "rdl_RNAseq.metadata.tsv"
metadata <- read_tsv(CLIN) %>%
                   dplyr::filter(! sampleID %in% c("S2","S12")) %>%
                   dplyr::mutate(path = file.path(sampleID, "kallisto", "abundance.h5")) %>%
                   dplyr::rename("sample" = sampleID)

##source functions
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/plot/poissonHeatmap.func.R")
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/sleuth/plot_expression_sleuth.func.R")
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/plot/pcaPlot.func.R")

##input parameters for loading sleuth data
memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))
ramPerCoreGB <- 15
threadAlloc <- round((memfree/1000000)/ramPerCoreGB, 0)-1

##rat annotation
OUTPUT <- paste0(RDATADIR, "/tx2gene.gene2ext.human_rat_intersect_anno.RData")
if(!file.exists(OUTPUT)){
    RGENOME <- "rnorvegicus_gene_ensembl"
    HGENOME <- "hsapiens_gene_ensembl"
    rmart <- biomaRt::useMart(biomart = "ensembl", dataset = RGENOME)
    hmart <- biomaRt::useMart(biomart = "ensembl", dataset = HGENOME)

    tx2gene <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = hmart)
    colnames(tx2gene)[1] <- "target_id"
    gene2ext <- biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = hmart)

    rtx2gene <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id",    "external_gene_name"), mart = rmart)) %>%
                dplyr::arrange(external_gene_name) %>%
                dplyr::mutate(human_external_gene_name = toupper(external_gene_name))
    rtx2gene_tid <- as_tibble(biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = rmart))
    colnames(rtx2gene_tid)[1] <- "target_id"

    htx2gene <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = hmart)) %>%
                dplyr::arrange(external_gene_name) %>%
                dplyr::rename(human_external_gene_name = "external_gene_name")
    rhtx2gene <- inner_join(rtx2gene, htx2gene, by="human_external_gene_name")
    rhtx2gene <- rhtx2gene %>% dplyr::rename(ensembl_gene_id = "ensembl_gene_id.x", human_ensembl_gene_id = "ensembl_gene_id.y")
    save(rtx2gene_tid, rtx2gene, htx2gene, rhtx2gene, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##read sleuth inputs
OUTPUT <- paste0(RDATADIR, "/rdl_RNAseq.sleuthObject.RData")
if(!file.exists(OUTPUT)){
    h5closeAll()
    so <- sleuth_prep(metadata,
                      full_model =~ Tissue + Diet + Batch,
                      target_mapping = rtx2gene_tid,
                      aggregation_column = "ensembl_gene_id",
                      gene_mode = TRUE,
                      num_cores = threadAlloc,
                      max_bootstrap = 50)
    save(so, metadata, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##counts from gene_mode (non-integer) sleuth object
OUTPUT <- paste0(RDATADIR, "/obs_norm_filt_df.conds.all.RData")
if(!file.exists(OUTPUT)){

    ##create DF
    obs_norm_filt_df <- dcast(so$obs_norm_filt,
                              target_id ~ sample,
                              value.var="scaled_reads_per_base") %>%
                        dplyr::rename("ensembl_gene_id" = target_id) %>%
                        dplyr::arrange(ensembl_gene_id) %>%
                        dplyr::mutate_if(is.numeric, round, 0) %>%
                        as.data.frame() %>%
                        column_to_rownames(., var="ensembl_gene_id")

    ##input conditions
    conds <- metadata %>% dplyr::select(sample, Individual, Tissue, Intervention, Diet, Batch) %>% as.data.frame()
    conds[] <- lapply(conds, factor)

    save(obs_norm_filt_df, conds, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##distribution of log(scaled_reads_per_base) to filter low expression
OUTPUT <- paste0(RDATADIR, "/obs_norm_filt2_df.conds.RData")
if(!file.exists(OUTPUT)){

    ##v. similar distribution for Tissue, Batch and other factors
    ##indciates log(3) as cutoff generally, =~20 scaled_reads_per_base
    samplet <- conds$sample[conds$Tissue=="Tumour"]
    sot <- so$obs_norm_filt %>% dplyr::filter(sample %in% samplet)
    samplen <- conds$sample[conds$Tissue=="Normal"]
    son <- so$obs_norm_filt %>% dplyr::filter(sample %in% samplen)

    ggp <- ggplot() +
           geom_density(data=son, aes(log(scaled_reads_per_base), colour="Tumour")) +
           geom_density(data=sot, aes(log(scaled_reads_per_base), colour="Normal")) +
           xlim(c(-0.5,11))

    ##filter obs_norm_filt_df
    print(paste0("Sleuth filtered: ", dim(obs_norm_filt_df), " genes"))
    obs_norm_filt2_df <- as_tibble(obs_norm_filt_df,
                                    rownames="ensembl_gene_id") %>% dplyr::mutate(mean = rowMeans(dplyr::select_if(., is.numeric))) %>%
                                                                    dplyr::filter(mean > 20) %>%
                                                                    dplyr::select(-mean) %>%
                                                                    column_to_rownames("ensembl_gene_id")
    print(paste0("Density filtered: ", dim(obs_norm_filt2_df), " genes"))

    save(sot, son, obs_norm_filt2_df, conds, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##limma-voom including duplicateCorrelation as multiple Individuals have several replicates
OUTPUT <- paste0(RDATADIR, "/limma-voom.DGEList-corfit.RData")
if(!file.exists(OUTPUT)){

    ##make input model.matrix design
    mmdf <- data.frame(Contrasts=paste(conds$Tissue, conds$Diet, conds$Intervention, sep="."),
                       Batch=conds$Batch)
    all.design <- model.matrix(~0 + Contrasts + Batch, data = mmdf)
    colnames(all.design) <- c(levels(mmdf$Contrasts), "Batch")

    ##make DGEList from filtered counts
    all.dge <- DGEList(counts=obs_norm_filt2_df)
    all.keep <- filterByExpr(all.dge,
                             all.design)
    all.dge <- all.dge[all.keep,
                       keep.lib.sizes=FALSE]
    all.dge <- calcNormFactors(all.dge,
                               method="TMM")

    ##voom, with 2 rounds of dupcor
    all.voom <- voom(all.dge,
                     all.design,
                     plot = FALSE)
    all.cor <- duplicateCorrelation(all.voom,
                                    block = conds$Individual)
    print(paste0("Duplicate correlation 1: ", all.cor$consensus.correlation))
    all.voom1 <- voom(all.dge,
                      all.design,
                      plot = FALSE,
                      correlation = all.cor$consensus.correlation, block=conds$Individual)
    all.cor1 <- duplicateCorrelation(all.voom1,
                                     block = conds$Individual)
    print(paste0("Duplicate correlation 2: ", all.cor1$consensus.correlation))

    ##fit
    all.fit <- lmFit(all.voom1,
                     all.design,
                     block = conds$Individual,
                     correlation = all.cor1$consensus.correlation)

    save(mmdf, all.design, all.voom1, all.cor1, all.fit, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##contrasts applied to tumour, normal across desired comparisons
OUTPUT <- paste0(RDATADIR, "/limma-voom.contrastList.RData")
if(!file.exists(OUTPUT)){

    cont.matrix <- makeContrasts(
        t_Diet_HFD_vs_NC=Tumour.HFD.NO-Tumour.NC.NO,
        t_Int_LIR_vs_NO=Tumour.HFD.LIR-Tumour.HFD.NO,
        t_Int_RES_vs_NO=Tumour.HFD.RES-Tumour.HFD.NO,
        t_Int_LIR_vs_RES=Tumour.HFD.LIR-Tumour.HFD.RES,
        n_Diet_HFD_vs_NC=Normal.HFD.NO-Normal.NC.NO,
        n_Int_LIR_vs_NO=Normal.HFD.LIR-Normal.HFD.NO,
        n_Int_RES_vs_NO=Normal.HFD.RES-Normal.HFD.NO,
        n_Int_LIR_vs_RES=Normal.HFD.LIR-Normal.HFD.RES,
        levels=all.design)
    all.fit2 <- contrasts.fit(all.fit, cont.matrix)
    all.fite <- eBayes(all.fit2, robust=TRUE)
    contrastList <- as.list(1:8)
    contrastList[[1]] <- topTable(all.fite, number=Inf, coef="t_Diet_HFD_vs_NC")
    contrastList[[2]] <- topTable(all.fite, number=Inf, coef="t_Int_LIR_vs_NO")
    contrastList[[3]] <- topTable(all.fite, number=Inf, coef="t_Int_RES_vs_NO")
    contrastList[[4]] <- topTable(all.fite, number=Inf, coef="t_Int_LIR_vs_RES")
    contrastList[[5]] <- topTable(all.fite, number=Inf, coef="n_Diet_HFD_vs_NC")
    contrastList[[6]] <- topTable(all.fite, number=Inf, coef="n_Int_LIR_vs_NO")
    contrastList[[7]] <- topTable(all.fite, number=Inf, coef="n_Int_RES_vs_NO")
    contrastList[[8]] <- topTable(all.fite, number=Inf, coef="n_Int_LIR_vs_RES")
    names(contrastList) <- c("t_Diet_HFD_vs_NC",
                             "t_Int_LIR_vs_NO",
                             "t_Int_RES_vs_NO",
                             "t_Int_LIR_vs_RES",
                             "n_Diet_HFD_vs_NC",
                             "n_Int_LIR_vs_NO",
                             "n_Int_RES_vs_NO",
                             "n_Int_LIR_vs_RES")
    save(cont.matrix, all.fit2, all.fite, all.cor1, contrastList, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##iterate over contrastList, find those with DE
OUTPUT <- paste0(RDATADIR, "/limma-voom.DEList.RData")
if(!file.exists(OUTPUT)){
    DEList <- lapply(seq_along(contrastList),function(f){
        dimcl <- as_tibble(contrastList[[f]], rownames="ensembl_gene_id") %>%
                    inner_join(rtx2gene,.) %>%
                    dplyr::filter(adj.P.Val < 0.1) %>%
                    dplyr::select(ensembl_gene_id, external_gene_name, logFC, adj.P.Val) %>%
                    dplyr::mutate_if(is.numeric, round, 3) %>%
                    dplyr::arrange(desc(logFC))
        if(dim(dimcl)[1] > 0){
            return(dimcl)
        }
    })
    names(DEList) <- names(contrastList)
    DEList <- DEList[!sapply(DEList, is.null)]
    lapply(seq_along(DEList), function(f){
        write_excel_csv(x=DEList[[f]],
                        path=paste0(OUTDIR, "/", names(DEList)[f],".DE.csv"),
                        col_names=TRUE,,
                        append=FALSE)
    })
    save(DEList, file=OUTPUT)

}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##native plots
all.fite$names <- as_tibble(all.fite$coefficients, rownames="ensembl_gene_id") %>%
                  left_join(., rhtx2gene) %>%
                  dplyr::select(-human_ensembl_gene_id) %>%
                  unique() %>%
                  dplyr::select(external_gene_name) %>%
                  unlist()

pv1 <- function(f){
    if(f<0.1){
        return("Y")
    }
    else{
        return("N")
    }
}

cindVec <- c(1, 4, 6, 8)
contVec <- c("Diet_HFD_vs_NC", "Int_LIR_vs_RES", "nInt_LIR_vs_NO", "nInt_LIR_vs_RES")
deVec <- c(1:4)

setwd(OUTDIR)
for (x in 1:length(cindVec)){
    cl1 <- as_tibble(contrastList[[cindVec[x]]], rownames="rn") %>% dplyr::arrange(rn) %>% dplyr::mutate(HL = unlist(lapply(adj.P.Val, pv1))) %>% dplyr::select(rn, HL) %>% column_to_rownames("rn") %>% unlist()

    pdf(paste0("volcanoMAplot.", contVec[x], ".pdf"))
        volcanoplot(all.fite, coef = cindVec[x], style = "p-value", highlight = dim(DEList[[deVec[x]]])[1], hl.col="red",
                xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, names=all.fite$names)
    plotMA(all.fite, coef = cindVec[x], hl.col="red",
           xlab = "Log2 Fold Change", pch=16, cex=0.35,
           ylab = "Expression log-ratio", status=cl1)
    dev.off()
}

##expression plots
dindVec <- c(1, 4, 6, 8)
dentList <- list(c("Diet","HFD","NC", "Tumour tissue; Diet: High Fat vs. Normal Chow", ""),
                 c("Intervention", "LIR", "RES", "Tumour tissue; Intervention: Liraglutide vs. Restriction", ""),
                 c("Intervention", "LIR", "NO", "Normal tissue; Intervention: Liraglutide vs. None", "n"),
                 c("Intervention", "LIR", "RES", "Normal tissue; Intervention: Liraglutide vs. Restriction", "n"))
deVec <- c(1:4)
lapply(seq_along(deVec), function(x){
    tfive <- left_join(DEList[[x]], rhtx2gene) %>% na.omit() %>% head(5) %>% dplyr::select(ensembl_gene_id, external_gene_name)
    bfive <- left_join(DEList[[x]], rhtx2gene) %>% na.omit() %>% tail(5) %>% dplyr::select(ensembl_gene_id, external_gene_name)

    ggtfive <- plot_expression_so_tpm(so, genes=tfive[,1], mapgenes=tfive, conds=metadata, condcol=dentList[[x]][1], factors=c(dentList[[x]][2], dentList[[x]][3]), tag=paste0("Top 5; ", dentList[[x]][4]))
    ggsave(ggtfive, filename=paste0("top5.", dentList[[x]][5], dentList[[x]][1],".",dentList[[x]][2], "_vs_", dentList[[x]][3], ".pdf"))

    ggbfive <- plot_expression_so_tpm(so, genes=bfive[,1], mapgenes=bfive, conds=metadata, condcol=dentList[[x]][1], factors=c(dentList[[x]][2], dentList[[x]][3]), tag=paste0("Bottom 5; ", dentList[[x]][4]))
    ggsave(ggbfive, filename=paste0("bottom5.", dentList[[x]][5], dentList[[x]][1],".",dentList[[x]][2], "_vs_", dentList[[x]][3], ".pdf"))
})
