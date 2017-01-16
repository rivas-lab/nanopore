########################################################
### Declare paths
########################################################
SAMTOOLS='/srv/gsfs0/projects/bustamante/progs/samtools-0.1.20/samtools'
GENOME_GRCH37='/srv/gsfs0/projects/bustamante/hcosta_projects/resources/Gencode/Release19_GRCh37/GRCh37.p13.genome.fa'
GENOME_BCRABL_DIR='/srv/gsfs0/projects/bustamante/hcosta_projects/resources/BCR-ABL_genome'
FASTQ='/srv/gsfs0/projects/bustamante/hcosta_projects/projects/minion/data/20161001_bcrabl/20161001_bcrabl.fq'
LAST='/srv/gsfs0/projects/bustamante/progs/last-759'
RESULTS_DIR='/srv/gsfs0/projects/bustamante/hcosta_projects/projects/minion/results'
PICARD='/srv/gsfs0/projects/bustamante/progs/picard-tools-1.138/picard.jar'


########################################################
### QLOGIN
########################################################

screen
qlogin -l h_vmem=40G -q extended


########################################################
### Install poretools
########################################################
# install anaconda2.7
# this is now default python, to load scg python do 'module load python'

git clone https://github.com/arq5x/poretools
cd poretools
python setup.py install


########################################################
### Get fastq 
########################################################

poretools \
	fastq \
	/srv/gsfs0/projects/bustamante/hcosta_projects/projects/minion/data/20161001_bcrabl/reads \
	> /srv/gsfs0/projects/bustamante/hcosta_projects/projects/minion/data/20161001_bcrabl/201601001_bcrabl.fq


########################################################
### Get run stats
########################################################

poretools \
	fastq \
	/srv/gsfs0/projects/bustamante/hcosta_projects/projects/minion/data/20161001_bcrabl/reads \
	> /srv/gsfs0/projects/bustamante/hcosta_projects/projects/minion/data/20161001_bcrabl/stats.txt


########################################################
### Align to BCR and ABL1 transcripts (follow LAST split DNA-genome alignment protocol as we're mapping to transcriptome)
########################################################

### 1) Concatenate BCR and ABL transcript fasta sequences
# ABL1 (https://www.ncbi.nlm.nih.gov/nuccore/NM_005157.5)
# BCR (https://www.ncbi.nlm.nih.gov/nuccore/NM_004327.3)
# Merge into one fa file
$GENOME_BCRABL_DIR/transcriptome.bcr.abl.loci.whole.fa


### 2) Index merged fasta file
$SAMTOOLS faidx \
	$GENOME_BCRABL_DIR/transcriptome.bcr.abl.loci.whole.fa


### 3) Create LAST index
$LAST/src/lastdb \
	-uNEAR \
	-R01 \
	$GENOME_BCRABL_DIR/transcriptome.bcr.abl.loci.whole \
	$GENOME_BCRABL_DIR/transcriptome.bcr.abl.loci.whole.fa


### 4) Align
$LAST/src/lastal \
	-s 2 \
	-T 0 \
	-Q 1 \
	-a 1 \
	-D100 \
	$GENOME_BCRABL_DIR/transcriptome.bcr.abl.loci.whole $FASTQ \
	| \
	$LAST/src/last-split \
	-c0 \
	> $RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.maf


########################################################
### Mask repetitive alignments
########################################################

$LAST/scripts/last-postmask \
	$RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.maf \
	> $RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.maf


########################################################
### Convert to bam, sort and index
########################################################

### 1) MAF to SAM
$LAST/scripts/maf-convert \
	sam \
	-r 'ID:BCR-ABL CN:BUSTAMANTE LB:BCR-ABL PL:NANOPORE PU:BCR-ABL SM:BCR-ABL' \
	$RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.maf \
	> $RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.sam


### 2) SAM to BAM
$SAMTOOLS \
	view \
	-bt \
	$GENOME_BCRABL_DIR/transcriptome.bcr.abl.loci.whole.fa.fai \
	$RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.sam \
	> $RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.bam


### 3) Sort and index bam
java -jar $PICARD SortSam \
	I=$RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.bam \
	O=$RESULTS_DIR/20161009_bcrabl_mapping_transcriptome_bcr_abl_loci_whole/transcriptome.bcr.abl.loci.whole.fastParam.masked.sorted.bam \
	SORT_ORDER=coordinate \
	CREATE_INDEX=TRUE

