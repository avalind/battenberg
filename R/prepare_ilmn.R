
# the tabfile should be generated using
# bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ALLELE_A\t%INFO/ALLELE_B\t[ %GT]\t[ %BAF]\t[ %LRR]\n' milsch.bcf > taboutput.tsv
# for single sample vcfs.
# TODO: for the corner case where the normal and tumor sample have been assayed on different
# arrays, we need to add a check so that we only keep SNPs present on both chips.
prepare.impute2.ilmn <- function(tabfile, outfile.start, chrom, sex, imputeinfo.file) {
    snp.cols <- c("chrom", "pos", "id",
        "ref", "alt", "allele_a", "allele_b",
        "genotype", "baf", "lrr")

    snp.data <- read.table(tabfile, sep="\t", header=F, 
        stringsAsFactors=FALSE, strip.white=TRUE)
    colnames(snp.data) <- snp.cols

    impute.info <- parse.imputeinfofile(imputeinfo.file, sex, chrom=chrom)
    known_SNPs <- read.table(impute.info$impute_legend[1], sep=" ", header=T)

    snp.data <- snp.data[snp.data[,1]==chrom,]
    snp.data$a.allele <- ifelse(snp.data$allele_a==1,
        snp.data$alt, snp.data$ref)
    snp.data$b.allele <- ifelse(snp.data$allele_b==1,
        snp.data$alt, snp.data$ref)
    
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
    het.snps <- snp.data.annot[is.het, c("pos", "a.allele", "b.allele", "eur.maf", "id.x", "pos", "a0", "a1")]
    colnames(het.snps) <- c("Physical.Position", "Allele.A", "Allele.B", "allele.frequency",
        "id", "position", "a0", "a1")
    write.csv(het.snps, file=paste0(outfile.start, chrom, "withAlleleFreq.csv", sep=""), quote=F, row.names=F)

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

    # add chromosome 23 special case
    return(het.snps)
}

# from battenberg: copied and pasted for devel.
parse.imputeinfofile = function(imputeinfofile, is.male, chrom=NA) {
  impute.info = read.table(imputeinfofile, stringsAsFactors=F)
  colnames(impute.info) = c("chrom", "impute_legend", "genetic_map", "impute_hap", "start", "end", "is_par")
  # Remove the non-pseudo autosomal region (i.e. where not both men and woman are diploid)
  if(is.male){ impute.info = impute.info[impute.info$is_par==1,] }
  chr_names=unique(impute.info$chrom)
  # Subset for a particular chromosome
  if (!is.na(chrom)) {
    impute.info = impute.info[impute.info$chrom==chr_names[chrom],]
  }
  return(impute.info)
}

prepare_illumina <- function(tumor_tabsep_file, normal_tabsep_file, tname,
                             chrom_names) {
    cat("Loading tumor data...")
    tumordata <- read.table(
        tumor_tabsep_file,
        sep="\t",
        header=F,
        stringsAsFactors=FALSE,
        strip.white=TRUE)

    cat("Processing tumor data...")
    tbaf <- tumordata[,c(1, 2, 9)]
    tlrr <- tumordata[,c(1, 2, 10)]

    cat("Writing tumor data...")
    # do we need further error checks?
    tumor.baf.out <- paste0(tname, "_mutantBAF.tab", sep="")
    write.table(tbaf, tumor.baf.out, row.names=F, quote=F, sep="\t")
    tumor.lrr.out <- paste0(tname, "_mutantLogR.tab", sep="")
    write.table(tlrr, tumor.lrr.out, row.names=F, quote=F, sep="\t")

    cat("Loading germline data...")
    germlinedata <- read.table(
        normal_tabsep_file,
        sep="\t",
        header=F,
        stringsAsFactors=FALSE,
        strip.white=TRUE)

    cat("Processing germline data...")

    nbaf <- germlinedata[,c(1, 2, 9)]
    nlrr <- germlinedata[,c(1, 2, 10)]

    normal.baf.out <- paste0(tname, "_germlineLogR.tab", sep="")
    write.table(nbaf, normal.baf.out, row.names=F, quote=F, sep="\t")

    normal.lrr.out <- paste(tname, "_germlineBAF.tab", sep="")
    write.table(nlrr, normal.lrr.out, row.names=F, quote=F, sep="\t")
    # TODO
}


