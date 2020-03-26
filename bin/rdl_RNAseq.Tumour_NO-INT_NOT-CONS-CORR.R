#! R

##load libraries
libs <- c("sleuth","tidyverse","biomaRt","pheatmap","RColorBrewer","PoiClaClu","DESeq2","rhdf5","ggplot2","devtools","cowplot", "reshape2", "limma", "edgeR", "DESeq2")
libsLoaded <- lapply(libs,function(l){suppressWarnings(suppressMessages(library(l, character.only = TRUE)))})

##input parameters
OUTDIR <- "kallisto-sleuth_any-no_corfit"
RDATADIR <- "kallisto-sleuth_any-no_corfit/RData"
dir.create(RDATADIR, showWarnings=FALSE, recursive=TRUE)

CLIN <- "rdl_RNAseq.metadata.tsv"
metadata <- read_tsv(CLIN) %>%
            dplyr::filter(! sampleID %in% c("S2","S12")) %>%
            dplyr::filter(! Tissue %in% c("Normal","Normal_Tumour")) %>%
            dplyr::mutate(path = file.path(sampleID, "kallisto", "abundance.h5")) %>%
            dplyr::rename("sample" = sampleID) %>%
            dplyr::mutate(Individual = paste0("R", Individual)) %>%
            dplyr::filter(Diet %in% "HFD") %>%
            dplyr::mutate(Int_Any = unlist(lapply(Intervention, function(f){
                                            if(f=="NO"){return("NO")}
                                            else{return("YES")}
                                    })))

##source plotting functions
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/plot/poissonHeatmap.func.R")
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/plot/pcaPlot.func.R")
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/sleuth/plot_expression_sleuth.func.R")

##input parameters for loading sleuth data
memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE))
ramPerCoreGB <- 15
threadAlloc <- round((memfree/1000000)/ramPerCoreGB, 0)-1

##txp2gene annotation
##rat-human annotation
OUTPUT <- paste0(RDATADIR, "/tx2gene.gene2ext.human_rat_intersect_anno.RData")
if(!file.exists(OUTPUT)){
    RGENOME <- "rnorvegicus_gene_ensembl"
    HGENOME <- "hsapiens_gene_ensembl"
    rmart <- biomaRt::useMart(biomart = "ensembl", dataset = RGENOME)
    hmart <- biomaRt::useMart(biomart = "ensembl", dataset = HGENOME)

    tx2gene <- biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = hmart)
    colnames(tx2gene)[1] <- "target_id"
    gene2ext <- biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = hmart)

    rtx2gene <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = rmart)) %>%
                dplyr::arrange(external_gene_name) %>%
                dplyr::mutate(human_external_gene_name = toupper(external_gene_name))
    rtx2gene_tid <- as_tibble(biomaRt::getBM(attributes=c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = rmart))
    colnames(rtx2gene_tid)[1] <- "target_id"

    htx2gene <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = hmart)) %>%
                dplyr::arrange(external_gene_name) %>%
                dplyr::rename(human_external_gene_name = "external_gene_name")
    rhtx2gene <- inner_join(rtx2gene, htx2gene, by="human_external_gene_name")
    rhtx2gene <- rhtx2gene %>% dplyr::rename(ensembl_gene_id = "ensembl_gene_id.x", human_ensembl_gene_id = "ensembl_gene_id.y")
    save(rtx2gene_tid, htx2gene, rhtx2gene, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##read sleuth inputs
OUTPUT <- paste0(RDATADIR, "/sleuthObject.no_corfit.RData")
if(!file.exists(OUTPUT)){
    h5closeAll()
    so <- sleuth_prep(metadata,
                      full_model =~ Batch,
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
OUTPUT <- paste0(RDATADIR, "/obs_norm_filt_df.conds.no_corfit.RData")
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
    conds <- metadata %>% dplyr::select(sample, Individual, Intervention, Int_Any, Batch) %>% as.data.frame()
    conds[] <- lapply(conds, factor)

    save(obs_norm_filt_df, conds, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##distribution of log(scaled_reads_per_base) to filter low expression
OUTPUT <- paste0(RDATADIR, "/obs_norm_filt2_df.conds.no_corfit.RData")
if(!file.exists(OUTPUT)){

    ##v. similar distribution for Tissue, Batch and other factors
    ##indciates log(3) as cutoff generally, =~20 scaled_reads_per_base
    sot <- so$obs_norm_filt
    ggp <- ggplot() +
           geom_density(data=sot, aes(log(scaled_reads_per_base), colour="Int_Any")) +
           xlim(c(-0.5,11))

    ##filter obs_norm_filt_df
    print(paste0("Sleuth filtered: ", dim(obs_norm_filt_df), " genes"))
    obs_norm_filt2_df <- as_tibble(obs_norm_filt_df,
                                    rownames="ensembl_gene_id") %>%
                         dplyr::mutate(mean = rowMeans(dplyr::select_if(., is.numeric))) %>%
                         dplyr::filter(mean > 20) %>%
                         dplyr::select(-mean) %>%
                         column_to_rownames("ensembl_gene_id")
    print(paste0("Density filtered: ", dim(obs_norm_filt2_df), " genes"))

    save(ggp, obs_norm_filt2_df, conds, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##limma-voom including duplicateCorrelation as multiple Individuals have several replicates
OUTPUT <- paste0(RDATADIR, "/limma-voom.DGEList.no_corfit.RData")
if(!file.exists(OUTPUT)){

    ##make input model.matrix design
    mmdf <- data.frame(Contrasts=conds$Int_Any)
    all.design <- model.matrix(~0 + Contrasts, data = mmdf)
    colnames(all.design) <- c(levels(mmdf$Contrasts))

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
                     all.design)

    ##fit
    all.fit <- lmFit(all.voom,
                     all.design)

    save(mmdf, all.design, all.voom, all.fit, file=OUTPUT)
}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

##contrasts applied to tumour, normal across desired comparisons
OUTPUT <- paste0(RDATADIR, "/limma-voom.contrastOut.no_corfit.RData")
if(!file.exists(OUTPUT)){

    cont.matrix <- makeContrasts(
        HFD_vs_ANY=NO-YES,
        levels=all.design)
    all.fit2 <- contrasts.fit(all.fit, cont.matrix)
    all.fite <- eBayes(all.fit2, robust=TRUE)
    contrastOut <- topTable(all.fite, number=Inf, coef="HFD_vs_ANY")

    contrastSig <- as_tibble(contrastOut, rownames="ensembl_gene_id") %>%
                     inner_join(rhtx2gene,.) %>%
                     dplyr::filter(adj.P.Val < 0.1) %>%
                     dplyr::select(ensembl_gene_id, external_gene_name, human_external_gene_name, human_ensembl_gene_id, logFC, adj.P.Val) %>%
                     dplyr::mutate_if(is.numeric, round, 3) %>%
                     dplyr::arrange(desc(logFC))

    save(contrastOut, contrastSig, all.fit2, all.fite, file=OUTPUT)

}
if(file.exists(OUTPUT)){
    load(OUTPUT)
}

write_excel_csv(x=contrastSig,
                path=paste0(OUTDIR, "/rat_diet_liraglutide.Tumour_HFD_NO-ANY.DE.NOTcorfit.csv"),
                col_names=TRUE,,
                append=FALSE)

all.fite$names <- rhtx2gene %>%
                  dplyr::filter(ensembl_gene_id %in% rownames(all.fite$coefficients)) %>%
                  dplyr::arrange(ensembl_gene_id) %>%
                  dplyr::select(external_gene_name) %>%
                  unlist()
setwd(OUTDIR)

pdf("rat-diet-liraglutide_RNAseq.Tumour_NO-INT_NOT-CONS-CORR.volcano.pdf")
    volcanoplot(all.fite, coef = "HFD_vs_ANY", style = "p-value", highlight = dim(contrastSig)[1], hl.col="red",
                xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, names=all.fite$names)
dev.off()

top.tibble.sig <- contrastSig %>% dplyr::filter(logFC %in% c(head(logFC, 10), tail(logFC, 10)))
NO_vs_ANYup <- unlist(top.tibble.sig[c(1:10), "ensembl_gene_id"])
NO_vs_ANYdown <- unlist(top.tibble.sig[c(11:20), "ensembl_gene_id"])


plot_expression_so_tpm(so, genes=NO_vs_ANYup, mapgenes=rhtx2gene, conds, condcol="Int_Any", factors=c("NO", "YES"), tag="Int_Any, Up in NO")
plot_expression_so_tpm(so, genes=NO_vs_ANYdown, mapgenes=rhtx2gene, conds, condcol="Int_Any", factors=c("NO", "YES"), tag="Int_Any, Down in NO")
