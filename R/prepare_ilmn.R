
# the tabfile should be generated using
# bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ALLELE_A\t%INFO/ALLELE_B\t[ %GT]\t[ %BAF]\t[ %LRR]\n' milsch.bcf > taboutput.tsv
# for single sample vcfs.
# TODO: for the corner case where the normal and tumor sample have been assayed on different
# arrays, we need to add a check so that we only keep SNPs present on both chips.
prepare.impute2.ilmn <- function(tabfile, tumorbaf, outfile.start, chrom, sex, imputeinfo.file) {
    snp.cols <- c("chrom", "pos", "id",
        "ref", "alt", "allele_a", "allele_b",
        "genotype", "baf", "lrr")

    snp.data <- read.table(tabfile, sep="\t", header=F, 
        stringsAsFactors=FALSE, strip.white=TRUE)
    colnames(snp.data) <- snp.cols
    
    tbaf <- read.table(tumorbaf, sep="\t", header=F, stringsAsFactors=FALSE)
    colnames(tbaf) <- c("chrom", "tbaf_pos", "tbaf")
    tbaf <- tbaf[tbaf$chrom==chrom,]    

    impute.info <- parse.imputeinfofile(imputeinfo.file, sex, chrom=chrom)
    known_SNPs <- read.table(impute.info$impute_legend[1], sep=" ", header=T, stringsAsFactors=FALSE)

    snp.data <- snp.data[snp.data[,1]==chrom,]
    snp.data$a.allele <- ifelse(snp.data$allele_a==1,
        snp.data$alt, snp.data$ref)
    snp.data$b.allele <- ifelse(snp.data$allele_b==1,
        snp.data$alt, snp.data$ref)

    snp.data <- merge(x=snp.data, y=tbaf, by.x="pos", by.y="tbaf_pos")   
    # fetch those in 1kG
    snp.data <- snp.data[snp.data$pos %in% known_SNPs$position,]

    # then subset 1kG to keep only info on those covered on array.
    known_SNPs <- known_SNPs[known_SNPs$position %in% snp.data$pos,]
    known_SNPs <- known_SNPs[order(snp.data$pos),]

    # is eur.maf the correct allele freq?
    snp.data.annot <- merge(x=snp.data, y=known_SNPs, by.x="pos", by.y="position")
   
    is.het <- snp.data.annot$genotype == "0/1"
    is.hom.ref <- snp.data.annot$genotype=="0/0"
    is.hom.alt <- snp.data.annot$genotype=="1/1"
    het.snps <- snp.data.annot[is.het, c("pos", "a.allele", "b.allele", "tbaf", "id.x", "pos", "a0", "a1")]
    colnames(het.snps) <- c("Physical.Position", "Allele.A", "Allele.B", "allele.frequency",
        "id", "position", "a0", "a1")
    
    # flip the tbaf if necessary (if a0 != Allele.a)
    new_freq = ifelse(het.snps[,2]==het.snps[,7], as.numeric(het.snps[,4]), 1.0-as.numeric(het.snps[,4]))
    het.snps[,4] = new_freq

    write.csv(het.snps, file=paste0(outfile.start, chrom, "_withAlleleFreq.csv", sep=""), quote=F, row.names=F)

    genotypes <- array(0, c(nrow(snp.data.annot), 3))
    genotypes[is.hom.ref, 1] = 1
    genotypes[is.het, 2] = 1
    genotypes[is.hom.alt, 3] = 1
    is.genotyped = (is.het | is.hom.ref | is.hom.alt)
    snp.names <- paste0("snp",1:sum(is.genotyped), sep="")
    out.data <- cbind(snp.names,
        snp.data.annot[is.genotyped, c("id.x","pos","a0","a1")],     
        genotypes[is.genotyped,])
    write.table(out.data, file=paste0(outfile.start,chrom,".txt"), row.names=F, col.names=F, quote=F)

    if (chrom == 23) {
        sample.g.file = paste0(outfile.start,"sample_g.txt")
        sample_g_data = data.frame(ID_1=c(0,"INDIVI1"),ID_2=c(0,"INDIVI1"),missing=c(0,0),sex=c("D",2))
        write.table(sample_g_data, file=sample.g.file, row.names=F, col.names=T, quote=F)
    }
}

prepare_ilmn <- function(tumor_tabsep_file, normal_tabsep_file, tname,
                             chrom_names) {
    cat("Loading tumor data...\n")
    tumordata <- read.table(
        tumor_tabsep_file,
        sep="\t",
        header=F,
        stringsAsFactors=FALSE,
        strip.white=TRUE)

    cat("Processing tumor data...\n")
    tbaf <- tumordata[,c(1, 2, 9)]
    tlrr <- tumordata[,c(1, 2, 10)]

    cat("Writing tumor data...\n")
    # do we need further error checks?
    tumor.baf.out <- paste0(tname, "_mutantBAF.tab", sep="")
    write.table(tbaf, tumor.baf.out, row.names=F, quote=F, sep="\t")
    tumor.lrr.out <- paste0(tname, "_mutantLogR.tab", sep="")
    write.table(tlrr, tumor.lrr.out, row.names=F, quote=F, sep="\t")

    cat("Loading germline data...\n")
    germlinedata <- read.table(
        normal_tabsep_file,
        sep="\t",
        header=F,
        stringsAsFactors=FALSE,
        strip.white=TRUE)

    cat("Processing germline data...\n")

    nbaf <- germlinedata[,c(1, 2, 9)]
    nlrr <- germlinedata[,c(1, 2, 10)]

    cat("Writing germline data...\n")
    normal.baf.out <- paste0(tname, "_germlineLogR.tab", sep="")
    write.table(nbaf, normal.baf.out, row.names=F, quote=F, sep="\t")

    normal.lrr.out <- paste0(tname, "_germlineBAF.tab", sep="")
    write.table(nlrr, normal.lrr.out, row.names=F, quote=F, sep="\t")
    # TODO
}


