
# Import
import os
import sys

# Fastq List
counter = 1
fl = "./input_alignment/"
il = "/projects/ag-papan/eduardo/Papantonis_Stag/Data/"
fastqList = ["control/fasta/HCT116_ChIP-seq_Control_BROAD_SE_1_ENCFF413RQG", "control/fasta/HCT116_ChIP-seq_Control_BROAD_SE_2.1_ENCFF957PGR", "control/fasta/HCT116_ChIP-seq_Control_BROAD_SE_2.2_ENCFF796ERZ", "control/fasta/HCT116_ChIP-seq_Control_HAIB_SE_1_ENCFF000PBW", "control/fasta/HCT116_ChIP-seq_Control_HAIB_SE_2_ENCFF000PBY", "control/fasta/HCT116_ChIP-seq_Control_USC_SE_1_ENCFF000VCY", "control/fasta/HCT116_ChIP-seq_Control_USC_SE_2_ENCFF000VCW", "control/fasta/HCT116_ChIP-seq_Control_UW_SE_1_ENCFF001HME", "dnase/fasta/HCT116_DNase-seq_UW_SE_1_ENCFF001DCL", "dnase/fasta/HCT116_DNase-seq_UW_SE_2_ENCFF001DCK", "histones/fasta/HCT116_ChIP-seq_H2AFZ_BROAD_SE_1_ENCFF023NCH", "histones/fasta/HCT116_ChIP-seq_H3K20me1_BROAD_SE_1_ENCFF478VDD", "histones/fasta/HCT116_ChIP-seq_H3K27ac_USC_SE_1_ENCFF000VCQ", "histones/fasta/HCT116_ChIP-seq_H3K27ac_USC_SE_2_ENCFF000VCS", "histones/fasta/HCT116_ChIP-seq_H3K27me3_BROAD_SE_1_ENCFF851YMW", "histones/fasta/HCT116_ChIP-seq_H3K36me3_USC_SE_1_ENCFF002AAO", "histones/fasta/HCT116_ChIP-seq_H3K36me3_USC_SE_2_ENCFF002AAN", "histones/fasta/HCT116_ChIP-seq_H3K4me1_USC_SE_1_ENCFF000VCI", "histones/fasta/HCT116_ChIP-seq_H3K4me1_USC_SE_2_ENCFF000VCK", "histones/fasta/HCT116_ChIP-seq_H3K4me2_BROAD_SE_1_ENCFF943AHJ", "histones/fasta/HCT116_ChIP-seq_H3K4me2_BROAD_SE_2_1_ENCFF936MMN", "histones/fasta/HCT116_ChIP-seq_H3K4me2_BROAD_SE_2_2_ENCFF937OOL", "histones/fasta/HCT116_ChIP-seq_H3K4me3_UW_SE_1_ENCFF001FIS", "histones/fasta/HCT116_ChIP-seq_H3K4me3_UW_SE_2_ENCFF001FIZ", "histones/fasta/HCT116_ChIP-seq_H3K79me2_BROAD_SE_1_ENCFF517MPZ", "histones/fasta/HCT116_ChIP-seq_H3K9ac_USC_SE_1_ENCFF398KLJ", "histones/fasta/HCT116_ChIP-seq_H3K9ac_USC_SE_2_ENCFF408RRT", "histones/fasta/HCT116_ChIP-seq_H3K9me2_BROAD_SE_1_ENCFF190MTF", "histones/fasta/HCT116_ChIP-seq_H3K9me3_USC_SE_1_ENCFF002AAK", "histones/fasta/HCT116_ChIP-seq_H3K9me3_USC_SE_2_ENCFF002AAM", "tf/fasta/HCT116_ChIP-seq_ATF3_HAIB_SE_1_ENCFF000OYD", "tf/fasta/HCT116_ChIP-seq_ATF3_HAIB_SE_2_ENCFF000OYI", "tf/fasta/HCT116_ChIP-seq_CBX3_HAIB_SE_1_ENCFF000OYG", "tf/fasta/HCT116_ChIP-seq_CBX3_HAIB_SE_2_ENCFF000OYN", "tf/fasta/HCT116_ChIP-seq_CEBPB_HAIB_SE_1_ENCFF000OYS", "tf/fasta/HCT116_ChIP-seq_CEBPB_HAIB_SE_2_ENCFF000OYU", "tf/fasta/HCT116_ChIP-seq_CTCF_BROAD_SE_1_ENCFF345KQX", "tf/fasta/HCT116_ChIP-seq_CTCF_HAIB_SE_1_ENCFF000OYZ", "tf/fasta/HCT116_ChIP-seq_CTCF_HAIB_SE_2_ENCFF000OZC", "tf/fasta/HCT116_ChIP-seq_CTCF_UW_SE_1_ENCFF001HLV", "tf/fasta/HCT116_ChIP-seq_CTCF_UW_SE_2_ENCFF001HLW", "tf/fasta/HCT116_ChIP-seq_EGR1_HAIB_SE_1_ENCFF000OZI", "tf/fasta/HCT116_ChIP-seq_EGR1_HAIB_SE_2_ENCFF000OZL", "tf/fasta/HCT116_ChIP-seq_ELF1_HAIB_SE_1_ENCFF000OZU", "tf/fasta/HCT116_ChIP-seq_ELF1_HAIB_SE_2_ENCFF000OZS", "tf/fasta/HCT116_ChIP-seq_EZH2_BROAD_SE_2.1_ENCFF758IZH", "tf/fasta/HCT116_ChIP-seq_EZH2_BROAD_SE_2.2_ENCFF254OVT", "tf/fasta/HCT116_ChIP-seq_EZH2phosphoT487_BROAD_SE_1_ENCFF961WCQ", "tf/fasta/HCT116_ChIP-seq_EZH2phosphoT487_BROAD_SE_2_ENCFF005ONO", "tf/fasta/HCT116_ChIP-seq_FOSL1_HAIB_SE_1_ENCFF000OZV", "tf/fasta/HCT116_ChIP-seq_FOSL1_HAIB_SE_2_ENCFF000OZZ", "tf/fasta/HCT116_ChIP-seq_MAX_HAIB_SE_1_ENCFF000PAP", "tf/fasta/HCT116_ChIP-seq_MAX_HAIB_SE_2_ENCFF000PAQ", "tf/fasta/HCT116_ChIP-seq_POL2RA_USC_SE_1_ENCFF000WXI", "tf/fasta/HCT116_ChIP-seq_POL2RA_USC_SE_2_ENCFF000WXG", "tf/fasta/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB_SE_1_ENCFF000PBA", "tf/fasta/HCT116_ChIP-seq_POLR2AphosphoS5_HAIB_SE_2_ENCFF000PBG", "tf/fasta/HCT116_ChIP-seq_REST_HAIB_SE_1_ENCFF000PAW", "tf/fasta/HCT116_ChIP-seq_REST_HAIB_SE_2_ENCFF000PAX", "tf/fasta/HCT116_ChIP-seq_SIN3A_HAIB_SE_1_ENCFF000PCD", "tf/fasta/HCT116_ChIP-seq_SIN3A_HAIB_SE_2_ENCFF000PCE", "tf/fasta/HCT116_ChIP-seq_SP1_HAIB_SE_1_ENCFF000PCM", "tf/fasta/HCT116_ChIP-seq_SP1_HAIB_SE_2_ENCFF000PCT", "tf/fasta/HCT116_ChIP-seq_TCF7L2_USC_SE_1_ENCFF000WXH", "tf/fasta/HCT116_ChIP-seq_TCF7L2_USC_SE_2_ENCFF000WXJ", "tf/fasta/HCT116_ChIP-seq_TEAD4_HAIB_SE_1_ENCFF000PDC", "tf/fasta/HCT116_ChIP-seq_TEAD4_HAIB_SE_2_ENCFF000PDE", "tf/fasta/HCT116_ChIP-seq_USF1_HAIB_SE_1_ENCFF000PBK", "tf/fasta/HCT116_ChIP-seq_USF1_HAIB_SE_1_ENCFF000PDL", "tf/fasta/HCT116_ChIP-seq_USF1_HAIB_SE_2_ENCFF000PBP", "tf/fasta/HCT116_ChIP-seq_USF1_HAIB_SE_2_ENCFF000PDJ", "tf/fasta/HCT116_ChIP-seq_YY1_HAIB_SE_1_ENCFF000PDQ", "tf/fasta/HCT116_ChIP-seq_YY1_HAIB_SE_2_ENCFF000PDT", "tf/fasta/HCT116_ChIP-seq_ZBTB33_HAIB_SE_1_ENCFF000PEC", "tf/fasta/HCT116_ChIP-seq_ZBTB33_HAIB_SE_2_ENCFF000PDZ", "tf/fasta/HCT116_ChIP-seq_ZNF274_USC_SE_1_ENCFF002AAP", "tf/fasta/HCT116_ChIP-seq_ZNF274_USC_SE_2_ENCFF002AAU", "tf/fasta/hg19_HCT116_ChIP-seq_SRF_HAIB_SE_1_ENCFF000PCS", "tf/fasta/hg19_HCT116_ChIP-seq_SRF_HAIB_SE_2_ENCFF000PCV", "tf/fasta/HCT116_ChIP-seq_JUND_HAIB_SE_2_ENCFF000PAH", "tf/fasta/HCT116_ChIP-seq_JUND_HAIB_SE_1_ENCFF000PAI", "tf/fasta/HCT116_ChIP-seq_RAD21_HAIB_SE_1_ENCFF000PBK", "tf/fasta/HCT116_ChIP-seq_RAD21_HAIB_SE_2_ENCFF000PBP", ["tf/fasta/HCT116_ChIP-seq_ZFX_SYDH_PE_1_1_ENCFF811IJJ", "tf/fasta/HCT116_ChIP-seq_ZFX_SYDH_PE_1_2_ENCFF489XKU"], ["tf/fasta/HCT116_ChIP-seq_ZFX_SYDH_PE_2_1_ENCFF223VZW", "tf/fasta/HCT116_ChIP-seq_ZFX_SYDH_PE_2_2_ENCFF792PVA"]]

# Fastq Loop
for fastq in fastqList:

  # SE or PE
  if(isinstance(fastq, basestring)):

    # Auxiliary
    if("tf/fasta" in fastq): ol = "tf/bam/"
    elif("histones/fasta" in fastq): ol = "histones/bam/"
    elif("dnase/fasta" in fastq): ol = "dnase/bam/"
    elif("control/fasta" in fastq): ol = "control/bam/"

    # Parameters
    alignType = "SE"
    minQuality = "20"
    ncores = "10"
    fastqFileName = il+fastq+".fastq.gz"
    indexFileName = "/projects/ag-papan/genomes/BowtieIndexes/hg19.zip"
    tempLocation = "/scratch/eduardo/stag_align_se/"
    outputLocation = il+ol

    # Creating files
    inFileName = fl+str(counter)+".txt"
    inFile = open(inFileName,"w")
    inFile.write("\n".join([alignType, minQuality, ncores, fastqFileName, indexFileName, tempLocation, outputLocation]))
    inFile.close()
    counter += 1

  elif(isinstance(fastq, list)):

    # Auxiliary
    if("tf/fasta" in fastq[0]): ol = "tf/bam/"
    elif("histones/fasta" in fastq[0]): ol = "histones/bam/"
    elif("dnase/fasta" in fastq[0]): ol = "dnase/bam/"
    elif("control/fasta" in fastq[0]): ol = "control/bam/"

    # Parameters
    alignType = "PE"
    minQuality = "20"
    ncores = "10"
    fastqFileName = il+fastq[0]+".fastq.gz,"+il+fastq[1]+".fastq.gz"
    indexFileName = "/projects/ag-papan/genomes/BowtieIndexes/hg19.zip"
    tempLocation = "/scratch/eduardo/stag_align_pe/"
    outputLocation = il+ol

    # Creating files
    inFileName = fl+str(counter)+".txt"
    inFile = open(inFileName,"w")
    inFile.write("\n".join([alignType, minQuality, ncores, fastqFileName, indexFileName, tempLocation, outputLocation]))
    inFile.close()
    counter += 1


