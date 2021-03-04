SHELL := /bin/bash
export SHELLOPTS:=errexit:pipefail
.DELETE_ON_ERROR:

BIN := script
INPUT := data
OUTPUT := output
GFF := human_gff/GRCh38.p2_gene.gff
HGENOME := ~/human_genome_GRCh38.p2/GCF_000001405.28_GRCh38.p2_genomic.fna


all : $(OUTPUT)/$(sample)/$(sample)_ISBP_R1.fastq
# sickle quality trimming
$(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq : $(INPUT)/$(sample)_R1.fastq.gz $(INPUT)/$(sample)_R2.fastq.gz | mkdir-output-$(sample)
	sickle pe -q 20 -w 10 -l 75 -n -t sanger \
	-f <(gunzip -c $(INPUT)/$(sample)_R1.fastq.gz) \
	-r <(gunzip -c $(INPUT)/$(sample)_R2.fastq.gz) \
	-o $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq \
	-p $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq \
	-s $(OUTPUT)/$(sample)/$(sample)_sickle_single.fastq \
	> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim illumina reverse adapter in filtered R2 reads
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq
	cutadapt --trimmed-only -g CCGCTCCGTCCGACGACTCACTATA -O 23 -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim linker in filtered R2 reads
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA.fastq
	cutadapt --trimmed-only -g GTTATGGTACTT -O 11 -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# extract 5' template + primer sequences
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_5tnp.fasta : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK.fastq
	$(BIN)/template_primer_fasta.pl $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA.fastq $^ 5 >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# extract 5' template id
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_5t.fasta : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_5tnp.fasta
	cutadapt --trimmed-only -a GTTATGGTACTT -O 11 -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# retrieve template ids from _5t files
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_templates.txt : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_5t.fasta
	$(BIN)/retrieve_templateid_from_tnpfile.pl $^ $@  >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# retrieve paired R1 reads
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_templates.txt
	$(BIN)/retrievePairReads.pl $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $^ $@  >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# trim 3' LTR
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR.fastq : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair.fastq
	cutadapt --trimmed-only -g GTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA -O 59 -o $@ $^ >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# make trimmed 3'LTR R1 template file
$(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR_templates.txt : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR.fastq
	$(BIN)/makeR1Templates.pl $^ $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_templates.txt $@ $(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR_templates.fastq >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# retrieve paired R2 reads
$(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_pair.fastq : $(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR_templates.txt
	$(BIN)/retrievePairReads.pl $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK.fastq $^ $@  >> $(OUTPUT)/$(sample)/$(sample)_log.txt
# map to human genome
$(OUTPUT)/$(sample)/$(sample)_bwa_human.sam : $(OUTPUT)/$(sample)/$(sample)_R2_sickle_RA_LK_pair.fastq
	bwa mem -t 10 $(HGENOME) $(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR_templates.fastq $^ > $@
# parse sam file
$(OUTPUT)/$(sample)/$(sample)_bwa_human_parsed.csv : $(OUTPUT)/$(sample)/$(sample)_bwa_human.sam
	$(BIN)/parseSamUMIBreakpoint.pl $^ $@
# get consensus IS and breakpoint
$(OUTPUT)/$(sample)/$(sample)_template_consensus_IS_breakpoint.csv : $(OUTPUT)/$(sample)/$(sample)_bwa_human_parsed.csv	
	$(BIN)/getConsensusISBreakpoint.pl $^ $(OUTPUT)/$(sample)/$(sample)_R1_sickle_pair_3LTR_templates.txt $(GFF) $@
# retrieve R1 and R2 reads with IS and breakpoint
$(OUTPUT)/$(sample)/$(sample)_ISBP_R1.fastq : $(OUTPUT)/$(sample)/$(sample)_template_consensus_IS_breakpoint.csv	
	$(BIN)/retrieveISBPR1R2Reads.pl $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq $(OUTPUT)/$(sample)/$(sample)_bwa_human_parsed.csv $@ $(OUTPUT)/$(sample)/$(sample)_ISBP_R2.fastq


	
# identify template ID from R2 reads
#$(OUTPUT)/$(sample)/$(sample)_R2_tid.txt : $(OUTPUT)/$(sample)/$(sample)_R2_sickle.fastq
#	$(BIN)/locateTemplateIdFromR2.pl $^ $@ >> $(OUTPUT)/$(sample)/$(sample).log	
# trim primer and LTR in R1 reads
#$(OUTPUT)/$(sample)/$(sample)_R1_trimLTR.fastq : $(OUTPUT)/$(sample)/$(sample)_R2_tid.txt
#	$(BIN)/trimR1PrimerLTR.pl $(OUTPUT)/$(sample)/$(sample)_R1_sickle.fastq $^ $@ $(OUTPUT)/$(sample)/$(sample)_R1_tid.txt >> $(OUTPUT)/$(sample)/$(sample).log
# map to human genome
#$(OUTPUT)/$(sample)/$(sample)_R1_trimLTR_bwa_human.sam : $(OUTPUT)/$(sample)/$(sample)_R1_trimLTR.fastq
#	bwa mem -t 10 $(HGENOME) $^ > $@
# parse sam file
#$(OUTPUT)/$(sample)/$(sample)_R1_trimLTR_bwa_human_parsed.csv : $(OUTPUT)/$(sample)/$(sample)_R1_trimLTR_bwa_human.sam
#	$(BIN)/parseSam.pl $^ $@
# get consensus IS
#$(OUTPUT)/$(sample)/$(sample)_template_consensus_IS.csv : $(OUTPUT)/$(sample)/$(sample)_R1_trimLTR_bwa_human_parsed.csv	
#	$(BIN)/getConsensusIS.pl $^ $(OUTPUT)/$(sample)/$(sample)_R1_tid.txt $(GFF) $@
	
#	rm $(OUTPUT)/$(sample)/*.fastq
#	rm $(OUTPUT)/$(sample)/$(sample)_sickle_*

mkdir-output-%:
	mkdir -p $(OUTPUT)/$*

clean : 
	rm -rfv $(OUTPUT)/
	

