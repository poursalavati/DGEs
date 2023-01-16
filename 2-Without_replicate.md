### 1-load packages
    library(edgeR)
    library(biomaRt)
    library(RColorBrewer)
    library(gplots)
    library(Glimma)
    library(pheatmap)

### 2-Change working directory to where you saved Rmatrix file:
  	File>Change dir... 

### 3-Load Rmatrix
    counts <- read.delim("Rmatrix", row.names = 1)
    dim(counts)
    head(counts)

### 4-Select our targets
    counts_S2_FA <- counts[,5:6]
    head(counts_S2_FA)

### 5-Import to edgeR and create groups
    x <- DGEList(counts_S2_FA)
    samplenames <- substring(colnames(x), 0, nchar(colnames(x)))
    colnames(x) <- samplenames
    group <- as.factor(c("T1", "T2"))
    x$samples$group <- group

### 6-Add annotation, metadata using BiomaRt
    listMarts(host="plants.ensembl.org")
    ensembl=useMart(host="plants.ensembl.org", "plants_mart")
    listDatasets(ensembl)
    ensembl = useDataset("vvinifera_eg_gene",mart=ensembl)
    filters = listFilters(ensembl)
    attributes = listAttributes(ensembl)
    attributes[1:10,]
    filters[1:10,]

    q <- rownames(counts_S2_FA)
    length(q)
    dim(x)

    annot <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'percentage_gene_gc_content', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position'), filters = 'ensembl_gene_id', values = q, mart = ensembl)

    rownames(annot) <- annot$ensembl_gene_id
    x$genes <- annot [rownames(x),]
    x


### 7-Transformations from the raw-scale
    cpm <- cpm(x)
    lcpm <- cpm(x, log=TRUE)
    L <- mean(x$samples$lib.size) * 1e-6
    M <- median(x$samples$lib.size) * 1e-6
    c(L, M)
    summary(lcpm)

### 8-Removing genes that are lowly expressed
    keep.exprs <- filterByExpr(x, group=group)
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    lcpm.cutoff <- log2(10/M + 2/L)

### 9-Visualizing raw data and filtered data density (PDF)
    nsamples <- ncol(x)
    col <- brewer.pal(nsamples, "Paired")
    par(mfrow=c(1,2))
    pdf(file ="Fig1_log-CPM-density.pdf")
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(main="A. Raw data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
    }
    
    legend("topright", samplenames, text.col=col, bty="n")
    lcpm <- cpm(x, log=TRUE)
    plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    title(main="B. Filtered data", xlab="Log-cpm")
    abline(v=lcpm.cutoff, lty=3)
    for (i in 2:nsamples){
    den <- density(lcpm[,i])
    lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", samplenames, text.col=col, bty="n")
   
    dev.off()
    dev.off()

### 10-Normalising gene expression distributions
    x <- calcNormFactors(x, method = "TMM")
    x$samples$norm.factors
    lcpm <- cpm(x, log=TRUE)

### 11-Create design and contrast matrix 
    design <- model.matrix(~0+group)
    colnames(design) <- gsub("group", "", colnames(design))
    design

    contr.matrix <- makeContrasts(T1vsT2 = T1 - T2,levels = colnames(design))
    contr.matrix

### 12-Dispersion value
    Simply pick a reasonable dispersion value, based on your experience with similar data, and use that for exactTest or glmFit. Typical values for the common BCV (square-root dispersion) for datasets arising from well-controlled experiments are 0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates. 

    bcv <- 0.2

### 13-Differential expression analysis

#### 13.1-Exact Test
    et <- exactTest(x, dispersion=bcv^2)
    summary(decideTests(et))
    topTags(et)

#### 13.2-likelihood ratio tests (glmFit):
    fit <- glmFit(x,dispersion=bcv^2,design)
    lrt <- glmLRT(fit,contrast=contr.matrix)
    summary(decideTests(lrt))


#### 13.3-Using threshold
    EdgeR offers a rigorous statistical test for thresholded hypotheses under the GLM framework. It is analogous to TREAT but much more powerful than the original TREAT method. Given a fold-change (or log-fold-change) threshold, the thresholded testing can be done by calling the function glmTreat() on a DGEGLM object produced by either glmFit() or glmQLFit().

    tr <- glmTreat(fit, contrast=contr.matrix, lfc=1)
    summary(decideTests(tr))

### Depending on the type of experiment and our targets, we can continue with either of these two tests:
### Here I continue with tr:

    tab <- topTags(tr, n=Inf, p=0.05)$table
    dim(tab)
    UpGenes <- tab[ tab$logFC > 0, ]
    DownGenes <- tab[ tab$logFC < 0,]

    write.table(topTags(et,n=4000), file='All_diff_et.txt', sep = "\t", row.names = F)
    write.table(UpGenes, file='UP_et_tr.txt', sep = "\t", row.names = F)
    write.table(DownGenes, file='DOWN_et_tr.txt', sep = "\t", row.names = F)

### 14-Heatmap for top DGEs
    top_tr_50 <- topTags(tr,n=50)
    lcpm_top50 <- lcpm[rownames(top_tr_50),]

    pheatmap(scale(lcpm_top50),cellwidth = 49, cellheight = 4, fontsize = 3,scale = "column",show_colnames = T, show_rownames = T,legend = TRUE,cluster_rows=T, cluster_cols=F,labels_col = c("S2.FarmA.T1", "S2.FarmA.T2"),border_color = "NA",main = "Top S2 DGE Genes in Farm A: T1 vs T2",, angle_col = "0",fontsize_col = 6,filename = "Fig2_Heatmap_V_dsj.pdf")

### 15-Venn diagram for samples
    pdf(file ="Fig3_venn.pdf")
    vennDiagram(counts_S2_FA, circle.col = c("red", "green3"),names=c("S2.FarmA.T1", "S2.FarmA.T2"),counts.col=c("orange"))
    dev.off()
    dev.off()

### 16-Gene ontology using BiomaRt
    tr_forGO <- topTags(tr,n=1085)
    tr_id_GO <- as.character(tr_forGO$table$ensembl_gene_id)

    all_tr_GO <- getBM(attributes=c('ensembl_gene_id', 'external_transcript_name', 'source', 'external_synonym', 'go_id', 'name_1006', 'definition_1006', 'namespace_1003', 'entrezgene_id', 'transcript_length', 'description'), filters = 'ensembl_gene_id', values = tr_id_GO, mart = ensembl)

    mergedtr_GO <- data.frame(all_tr_GO, tr_forGO[match(all_tr_GO$ensembl_gene_id, tr_forGO$table$ensembl_gene_id),])

    write.table(mergedtr_GO, file = "Final_GO_tr.txt", sep = "\t", quote = FALSE, row.names=F)

### 17-Import to the Excel
    Import UP, Down, All_diff_tr.txt and Final_GO_tr.txt to seprate sheet and rename:
    No_Annotation and GO
    With_Annotation and GO (sort by Pvalue or FDR. Also remove endemble_id.01 column)


### 18-MD plot for DGEs
    pdf(file ="Fig4_plotMD.pdf")
    plotMD(tr)
    abline(h=c(-1, 1), col="blue")
    dev.off()
    dev.off()

### 19-Create and customize interactive graph using Glimma
    dt_tr <- decideTests(tr)
    summary(dt_tr)

    dim(mergedtr_GO)
    all_trid_GO <- as.character(tr$genes$ensembl_gene_id)
    length(all_trid_GO)

    dim(tr)

    all_tr_GO <- getBM(attributes=c('ensembl_gene_id', 'external_transcript_name', 'source', 'external_synonym', 'go_id', 'name_1006', 'definition_1006', 'namespace_1003', 'entrezgene_id', 'transcript_length', 'description'), filters = 'ensembl_gene_id', values = all_trid_GO, mart = ensembl)

    tr2 <- tr
    merged_all_tr_GO <- data.frame(all_tr_GO, tr2$genes[match(all_tr_GO$ensembl_gene_id, tr2$genes$ensembl_gene_id),])

    tr2$genes <- merged_all_tr_GO[rownames(tr2$genes),]
    colnames(tr2$genes)
    colnames(tr2$genes)[1] <- "Ensemb ID"
    colnames(tr2$genes)[6] <- "GO name"
    colnames(tr2$genes)[7] <- "GO definition"
    colnames(tr2$genes)[8] <- "GO domain"
    colnames(tr2$genes)[9] <- "NCBI ID"
    colnames(tr2$genes)[15] <- "Biotype"
    colnames(tr2$genes)[14] <- "GC"
    colnames(tr2$genes)

    glMDPlot(tr2, coef=1, status=dt_tr, main=colnames(tr2)[1], side.main="Ensembl.ID", counts=lcpm, groups=group, launch=FALSE, display.columns=c("Ensemb ID", "NCBI ID", "GC", "GO domain", "GO name", "GO definition", "Biotype"), folder="Intractive-MD_tr")

    write.table(summary(dt_tr), file='DGE_Summ_tr.txt', sep = "\t", row.names = T,quote = FALSE)

### Save workspace and history for each experiment

################################################
