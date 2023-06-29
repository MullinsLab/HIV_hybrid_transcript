SHELL := /bin/bash
export SHELLOPTS:=errexit:pipefail
.DELETE_ON_ERROR:

# pipeline identifying IS from 5'/3' ends

LTR ?= 3
ERRORRATE ?= 0.1
MINMAPLEN ?= 20
MINIDENTITY ?= 0.99
HIVMINIDENTITY ?= 0.9
THREADS ?= 10
MINQUALITY ?= 20
MINREADLEN ?= 75
BIN := script
INPUT := data
OUTPUT := output_$(LTR)LTR
GFF := human_gff/GRCh38.p2_gene.gff
HGENOME := human_genome/human_genome_GRCh38.p2/GCF_000001405.28_GRCh38.p2_genomic.fna
HXB2 := HXB2/HXB2-K03455.fasta
LTRFILE := seqs/$(LTR)ltr.fas
ADPTFILE := seqs/adapter.fas
LINKERFILE := seqs/linker.fas

TEXT := $(shell script/parse_seqs.pl $(LTRFILE) $(ADPTFILE) $(LINKERFILE))
LTRSEQ := $(word 1, $(TEXT))
ADPTSEQ := $(word 2, $(TEXT))
LINKERSEQ := $(word 3, $(TEXT))
LTROVLP := $(word 4, $(TEXT))
ADPTOVLP := $(word 5, $(TEXT))
LINKEROVLP := $(word 6, $(TEXT))

all : $(OUTPUT)/$(sample)/$(sample)_consensus_R1_R2_mapping2HIV.csv
# sickle quality trimming
$(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq : $(INPUT)/$(sample)_R1.fastq.gz $(INPUT)/$(sample)_R2.fastq.gz | mkdir-output-$(sample)
	sickle pe -q $(MINQUALITY) -l $(MINREADLEN) -n -t sanger \
	-f <(gunzip -c $(INPUT)/$(sample)_R1.fastq.gz) \
	-r <(gunzip -c $(INPUT)/$(sample)_R2.fastq.gz) \
	-o $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq \
	-p $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq \
	-s $(OUTPUT)/$(sample)/$(sample)_sickle_single.fastq \
	> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim LTR in filtered R1 reads
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR.fastq : $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq
	cutadapt --trimmed-only -g $(LTRSEQ) -O $(LTROVLP) -e $(ERRORRATE) -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# retrieve paired R2 reads
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR.fastq : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR.fastq
	$(BIN)/retrievePairReadsFromFastq.pl $^ $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq $@  >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim illumina reverse adapter in paired R2 reads
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR.fastq
	cutadapt --trimmed-only -g $(ADPTSEQ) -O $(ADPTOVLP) -e $(ERRORRATE) -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim linker in paired R2 reads
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA.fastq
	cutadapt --trimmed-only -g $(LINKERSEQ) -O $(LINKEROVLP) -e $(ERRORRATE) -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# retrieve paired R1 reads
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK.fastq
	$(BIN)/retrievePairReadsFromFastq.pl $^ $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR.fastq $@  >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# retrieve R1 trimmed sequences
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK_trimmed.csv : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK.fastq
	$(BIN)/retrieveTrimmedSeq.pl $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $^ $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK_trimmed.fasta $@  >> $(OUTPUT)/$(sample)/$(sample)_log.txt 
# retrieve R2 trimmed sequences
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed.fasta : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK_trimmed.csv
	$(BIN)/retrieveTrimmedSeq.pl $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK.fastq $@ $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed.csv  >> $(OUTPUT)/$(sample)/$(sample)_log.txt 
# trim illumina reverse adapter in _R2_sickle_$(LTR)LTR_RA_LK_trimmed.fasta
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_UMI_LK.fasta : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed.fasta
	cutadapt --trimmed-only -g $(ADPTSEQ) -O $(ADPTOVLP) -e $(ERRORRATE) -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim linker in _R2_sickle_$(LTR)LTR_UMI_LK.fasta
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_UMI.fasta : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_UMI_LK.fasta
	cutadapt --trimmed-only -a $(LINKERSEQ) -O $(LINKEROVLP) -e $(ERRORRATE) -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# append UMI to $(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed.csv
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed_umi.csv : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_UMI.fasta
	$(BIN)/append_umi.pl $^ $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed.csv $@
# map to human genome
$(OUTPUT)/$(sample)/$(sample)_bwa_human.sam : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed_umi.csv
	bwa mem -t $(THREADS) -o $@ $(HGENOME) $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK.fastq
# parse human sam file
$(OUTPUT)/$(sample)/$(sample)_bwa_human_parsed.csv : $(OUTPUT)/$(sample)/$(sample)_bwa_human.sam
	$(BIN)/parseSamBreakpoint.pl $^ $@ $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK_trimmed.csv $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed_umi.csv $(LTR) $(MINMAPLEN)
# get consensus IS and breakpoint mapping to human at least 99%
$(OUTPUT)/$(sample)/$(sample)_consensus_IS_BP.csv : $(OUTPUT)/$(sample)/$(sample)_bwa_human_parsed.csv	
	$(BIN)/getConsensusISBPWithIdentityFusionRepetitiveMapping.pl $^ $(GFF) $(LTR) $@ $(MINIDENTITY)
# get final IS summary file
$(OUTPUT)/$(sample)/$(sample)_consensus_IS_BP_Final.csv : $(OUTPUT)/$(sample)/$(sample)_consensus_IS_BP.csv	
	$(BIN)/add_quality_control.pl $^ $@
# map to HIV (HXB2) genome
$(OUTPUT)/$(sample)/$(sample)_bwa_hxb2.sam : $(OUTPUT)/$(sample)/$(sample)_consensus_IS_BP_Final.csv
	bwa mem -t $(THREADS) -o $@ $(HXB2) $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK.fastq
# parse hiv sam file
$(OUTPUT)/$(sample)/$(sample)_bwa_hxb2_parsed.csv : $(OUTPUT)/$(sample)/$(sample)_bwa_hxb2.sam
	$(BIN)/parseSamBreakpoint.pl $^ $@ $(OUTPUT)/$(sample)/$(sample)_R1_sickle_$(LTR)LTR_RA_LK_trimmed.csv $(OUTPUT)/$(sample)/$(sample)_R2_sickle_$(LTR)LTR_RA_LK_trimmed_umi.csv $(LTR) $(MINMAPLEN)
# get consensus R1 and R2 (trimmed all primers, linker and adaptors) mapping to HXB2 at least 90%
$(OUTPUT)/$(sample)/$(sample)_consensus_R1_R2_mapping2HIV.csv : $(OUTPUT)/$(sample)/$(sample)_bwa_hxb2_parsed.csv
	$(BIN)/getHIVOnlyConsensus.pl $^ $@ $(HIVMINIDENTITY)


mkdir-output-%:
	mkdir -p $(OUTPUT)/$*

clean : 
	rm -rfv $(OUTPUT)/


