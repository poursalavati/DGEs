### 1-If you have not R installed:
    https://cloud.r-project.org/bin/windows/base/R-4.1.0-win.exe

### 2-Insall and load packages:
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

    BiocManager::install(c("edgeR", "Glimma", "gplots", "RColorBrewer", "biomaRt"))

    install.packages("pheatmap")


#### Load packages:


    library(edgeR) 
    library(biomaRt)
    library(RColorBrewer)
    library(gplots)
    library(Glimma)
    library(pheatmap)


### 3-Change working directory to where you saved Rmatrix file:
    File>Change dir... 
	
### 4-Load Rmatrix
    counts <- read.delim("Rmatrix", row.names = 1)

    dim(counts)
    
    head(counts)

### 5-Select our targets
    counts_S3 <- counts[,c(9,11,10,12)]
    
    dim(counts_S3)  
    
    head(counts_S3)

### 6-Import to edgeR and create groups
    x <- DGEList(counts_S3)
    samplenames <- substring(colnames(x), 0, nchar(colnames(x)))

    x

    colnames(x) <- samplenames

    colnames(x)

    group <- as.factor(c("T1", "T1", "T2", "T2"))

    group

    x$samples$group <- group

### 7-Add annotation, metadata using BiomaRt   

    listMarts(host="plants.ensembl.org")   
.  

    bensembl=useMart(host="plants.ensembl.org", "plants_mart")

    listDatasets(ensembl)

    ensembl = useDataset("vvinifera_eg_gene",mart=ensembl)
    
    filters = listFilters(ensembl)
    
    attributes = listAttributes(ensembl)
    
    attributes[1:10,]
    
    filters[1:10,]
.  

    q <- rownames(counts)
    
    length(q)
    
    dim(x)
.

    annot <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'percentage_gene_gc_content', 'gene_biotype', 'chromosome_name', 'start_position', 'end_position'), filters = 'ensembl_gene_id', values = q, mart = ensembl)
.  

    rownames(annot) <- annot$ensembl_gene_id
    
    x$genes <- annot [rownames(x),]
    
    head(x)

### 8-Transformations from the raw-scale
    cpm <- cpm(x)
    
    lcpm <- cpm(x, log=TRUE)    
    
    L <- mean(x$samples$lib.size) * 1e-6
    
    M <- median(x$samples$lib.size) * 1e-6
    
    c(L, M)
    
    summary(lcpm)

### 9-Removing genes that are lowly expressed
    keep.exprs <- filterByExpr(x, group=group)
    
    x <- x[keep.exprs,, keep.lib.sizes=FALSE]

    lcpm.cutoff <- log2(10/M + 2/L)


### 10-Visualizing raw data and filtered data density 

####  (you can skip this step and go directly to PDF creation)
    lcpm2 <- lcpm
    
    x2 <- x
    
    nsamples <- ncol(x2)
    
    col <- brewer.pal(nsamples, "Paired")
    
    par(mfrow=c(1,2))
    
    plot(density(lcpm2[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    
    title(main="A. Raw data", xlab="Log-cpm")
    
    abline(v=lcpm.cutoff, lty=3)
    
    for (i in 2:nsamples){
        den <- density(lcpm2[,i]) 
        lines(den$x, den$y, col=col[i], lwd=2)
        }

    legend("topright", samplenames, text.col=col, bty="n")
    
    lcpm2 <- cpm(x2, log=TRUE)
    
    plot(density(lcpm2[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    
    title(main="B. Filtered data", xlab="Log-cpm")
    
    abline(v=lcpm.cutoff, lty=3)
    
    for (i in 2:nsamples){
        den <- density(lcpm2[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
        }
    
    legend("topright", samplenames, text.col=col, bty="n")

### 11-PDF output of the graph:
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

### 12-Normalising gene expression distributions
    x <- calcNormFactors(x, method = "TMM")
    
    x$samples$norm.factors

### 13-Create design and contrast matrix 
    design <- model.matrix(~0+group)
    
    colnames(design) <- gsub("group", "", colnames(design))
    
    design

    contr.matrix <- makeContrasts(T1vsT2 = T1 - T2,levels = colnames(design))

    contr.matrix

### 14-Estimate dispersion value
    x <- estimateDisp(x,design)

    x

### 15-MDS plot for samples
    lcpm <- cpm(x, log=TRUE)
    
    pdf(file ="Fig2_PlotMDS.pdf")
    
    plotMDS(x, col = as.numeric(group))
    
    title(main="Sample groups")

    dev.off()
    dev.off()

### 16-Differential expression analysis
    fit <- glmQLFit(x, design)
    
    qlf <- glmQLFTest(fit, contrast=contr.matrix)
    
    summary(decideTests(qlf))

    tabqlf <- topTags(qlf, n=Inf, p=0.05)$table

    dim(tabqlf)

    UpGenesqlf <- tabqlf[ tabqlf$logFC > 0, ]

    DownGenesqlf <- tabqlf[ tabqlf$logFC < 0,]

    write.table(topTags(qlf,n=180), file='All_diff_qlf.txt', sep = "\t", row.names = F)

    write.table(UpGenesqlf, file='UP_qlf.txt', sep = "\t", row.names = F)
    
    write.table(DownGenesqlf, file='DOWN_qlf.txt', sep = "\t", row.names = F)

### 17-Heatmap for top DGEs
    top_qlf_50 <- topTags(qlf,n=50)
    
    lcpm_qlf_top50 <- lcpm[rownames(top_qlf_50),]

    pheatmap(scale(lcpm_qlf_top50),cellwidth = 49, cellheight = 4, fontsize = 3,scale = "column",show_colnames = T, show_rownames = T,legend = TRUE,cluster_rows=T, cluster_cols=F,labels_col = c("S3.FarmA.T1","S3.FarmA.T2","S3.FarmB.T1","S3.FarmB.T2"),border_color = "NA",main = "Top S3 DGEs: T1 vs T2", angle_col = "0",fontsize_col = 6,filename = "Fig3_Heatmap.pdf") 

### 18-Gene ontology using BiomaRt
    qlf_forGO <- topTags(qlf,n=180)

    qlf_id_GO <- as.character(qlf_forGO$table$ensembl_gene_id)

    all_qlf_GO <- getBM(attributes=c('ensembl_gene_id', 'external_transcript_name', 'source', 'external_synonym', 'go_id', 'name_1006', 'definition_1006', 'namespace_1003', 'entrezgene_id', 'transcript_length', 'description'), filters = 'ensembl_gene_id', values = qlf_id_GO, mart = ensembl)

    mergedQLF_GO <- data.frame(all_qlf_GO, qlf_forGO[match(all_qlf_GO$ensembl_gene_id, qlf_forGO$table$ensembl_gene_id),])

    write.table(mergedQLF_GO, file = "Final_GO_QLF.txt", sep = "\t", quote = FALSE, row.names=F)

### 19-Import to the Excel
#### Import UP, Down, All_diff_qlf.txt and Final_GO_QLF.txt to seprate sheet and rename:
- No_Annotation and GO  
- With_Annotation and GO (sort by Pvalue or FDR. Also remove endemble_id.01 column)

### 20-MD plot for DGEs
    pdf(file ="Fig4_plotMD.pdf")
    
    plotMD(qlf)
    
    abline(h=c(-1, 1), col="blue")

    dev.off()
    dev.off()

### 21-Create and customize interactive graph using Glimma
    dt_qlf <- decideTests(qlf)
    
    summary(dt_qlf)

    dim(mergedQLF_GO)

    all_qlfid_GO <- as.character(qlf$genes$ensembl_gene_id)

    length(all_qlfid_GO)

    dim(qlf)

    all_qlf_GO <- getBM(attributes=c('ensembl_gene_id', 'external_transcript_name', 'source', 'external_synonym', 'go_id', 'name_1006', 'definition_1006', 'namespace_1003', 'entrezgene_id', 'transcript_length', 'description'), filters = 'ensembl_gene_id', values = all_qlfid_GO, mart = ensembl)

    qlf2 <- qlf

    merged_all_qlf_GO <- data.frame(all_qlf_GO, qlf2$genes[match(all_qlf_GO$ensembl_gene_id, qlf2$genes$ensembl_gene_id),])

    qlf2$genes <- merged_all_qlf_GO[rownames(qlf2$genes),]
    colnames(qlf2$genes)
    colnames(qlf2$genes)[1] <- "Ensemb ID"
    colnames(qlf2$genes)[6] <- "GO name"
    colnames(qlf2$genes)[7] <- "GO definition"
    colnames(qlf2$genes)[8] <- "GO domain"
    colnames(qlf2$genes)[9] <- "NCBI ID"
    colnames(qlf2$genes)[15] <- "Biotype"
    colnames(qlf2$genes)[14] <- "GC"
    colnames(qlf2$genes)
.  

    glMDPlot(qlf2, coef=1, status=dt_qlf, main=colnames(qlf2)[1], side.main="Ensembl.ID", counts=lcpm, groups=group, launch=FALSE, display.columns=c("Ensemb ID", "NCBI ID", "GC", "GO domain", "GO name", "GO definition", "Biotype"), folder="Intractive-MD_QLf")

###   22-Venn diagram for samples
    pdf(file ="Fig5_Venn.pdf")

    vennDiagram(cpm, circle.col = c("red", "green3"),counts.col=c("orange"))

    dev.off()
    dev.off()

    write.table(summary(dt_qlf), file='DGE_Summ_qlf.txt', sep = "\t", row.names = T,quote = FALSE)

### Save workspace and history for each experiment

################################################
