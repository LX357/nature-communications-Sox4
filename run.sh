# Create directories for output files
mkdir -p filterData rRNA/smallRNA align/exon/index

# Build combined index for rRNA and small RNA, then build bowtie2 index
perl fetch_ncgb.pl -t rRNA -s Homo GenBank_gbanimal.fa > ref4riboseq/rRNA/rRNA.fa
bowtie2-build ref4riboseq/rRNA/rRNA.fa ref4riboseq/rRNA/rRNA
cat Rfam.fasta ref4riboseq/hsa.ncgb.fa > rRNA/smallRNA/indexxxx.fa
bowtie2-build rRNA/smallRNA/indexxxx.fa rRNA/smallRNA/smallrnas.fa

# Build bowtie2 index for exon sequences
gffread ref4riboseq/hsa.gtf -g hsa.fa -w hsa_exon.fa -W
bowtie2-build ref4riboseq/hsa_exon.fa align/exon/index/exon.fa

# Process each sample in batch
cat sample.list | while read i
do 
	# Quality control and trimming raw reads using fastp
	fastp -i ${i}.fq.gz -o filterData/${i}_1.fq.gz -a TGGAATTCT -l 10 -q 20 -u 50 -w 4 -j filterData/$i.json -h filterData/$i.html

	# Align clean reads to rRNA index, remove rRNA reads
	bowtie2 -S /dev/null -p 4 --mm -x ref4riboseq/rRNA/rRNA -U filterData/${i}_1.fq.gz --un-gz rRNA/$i.fq.gz 2> rRNA/$i.log

	# Align non-rRNA reads to small RNA index, remove small RNA reads
	bowtie2 -p 4 --mm --no-unal -x rRNA/smallRNA/smallrnas.fa -U rRNA/$i.fq.gz --un-gz rRNA/smallRNA/$i.fq.gz -S rRNA/smallRNA/$i.sam 2> rRNA/smallRNA/$i.sRNA.log

	# Filter reads by length (20-40 nt)
	perl filter_length.pl rRNA/smallRNA/$i.fq.gz 20 40 $i

	# Align filtered reads to reference genome using STAR (riboseq mapping)
	STAR --runThreadN 4 --runMode alignReads --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic --outSAMtype BAM Unsorted --outSAMunmapped None --genomeDir ref4riboseq --readFilesIn rRNA/smallRNA/${i}_ft.fq.gz --outFileNamePrefix align/${i}.

	# Process BAM file and annotate ribosome location
	perl deal_sam_ribo_v2.pl align/$i.Aligned.out.bam ref4riboseq/hsa_rib_loc.bed align/$i

	# Align reads to exon sequences, convert to sorted BAM
	bowtie2 -p 4 -k 10 --no-unal -x align/exon/index/exon.fa -U rRNA/smallRNA/${i}_ft.fq.gz 2> align/exon/$i.log | samtools view -Sb - | samtools sort - align/exon/$i

	# Predict ribosome pausing sites
	perl offline_pausepred.m.pl align/exon/$i.bam 1000 20 ref4riboseq/hsa_exon.fa align/exon/$i.xls

	# Generate final stop table
	perl maketable.pl align/exon/$i.annot.xls ref4riboseq/hsa_exon.fa > align/exon/$i.stop.xls

	# Generate 3-nucleotide periodicity profile 
	perl draw_ud_v3.pl align/exon/$i.bam ref4riboseq/hsa.gtf
done
