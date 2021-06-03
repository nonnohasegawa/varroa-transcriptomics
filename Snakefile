##############################################################################
##############################################################################
##############################################################################
################ TRANSCRIPTOMIC ANALYSIS PIPELINE FOR BIG DATA ###############
######## TAKEN SNAKEMAKE FILE FROM MAEVA TECHER AND EDITED BY NONNO H  #######
##############################################################################
##############################################################################
##############################################################################

localrules: getHaps, all

### SET DIRECTORY PATHS REQUIRED
READDir = "/bucket/MikheyevU/Nonno/VarroaWorld/renamed"
OUTDir = "/flash/MikheyevU/Nonno/temp/trimmed_reads"
REFDir = "/flash/MikheyevU/Nonno/ref2020"
SCRATCH = "/flash/MikheyevU/Nonno/scratch"
DATADir = "/flash/MikheyevU/Nonno/temp/data"


### PATHS FOR VARROA DESTRUCTOR GENOME AND REGIONS SPLIT
VIRUSESRef = REFDir + "/viruses/viruses2020.fasta"
VIRUSrenamed = REFDir + "/viruses/ENArenamed/viruses2020-ENA.fasta"
renamedIndex = REFDir + "/viruses/ENArenamed/viruses-ENA"
VIRUSESBowtieIndex = REFDir + "/viruses/viruses"


### SAMPLES LIST AND OTHER PARAMETERS
SAMPLES, = glob_wildcards(OUTDir + "/{sample}_R1.fastq") #this one is the exact raw input
CONTIGS = ["CEND01000001.1", "KR819915.1", "KX578272.1", "KY354234.1", "KY354240.1", "MG571081.1", "MG571088.1", "MK032464.1", "MK032465.1", "MK032466.1", "MK032467.1", "MK032468.1", "MK032469.1", "MK032470.1", "NC_002066.1", "NC_002548.1", "NC_003784.1", "NC_004807.1", "NC_004830.2", "NC_006494.1", "NC_009025.1", "NC_010711.1", "NC_014137.1", "NC_027619.1", "NC_027631.1", "NC_032433.1", "NC_035071.1", "NC_040601.1"]


rule all:
        input:
        		#expand(DATADir + "/bowtie2/virus/ENA-renamed/alignments/{sample}.bam", sample = SAMPLES),
        		expand(DATADir + "/bowtie2/virus/ENA-renamed/stats/flagstat/{sample}.txt", sample = SAMPLES),
        		expand(DATADir + "/bowtie2/virus/ENA-renamed/stats/idxstats/{sample}.txt", sample = SAMPLES),
        		#DATADir + "/bowtie2/virus/ENA-renamed/varscan/ENA-sorted.mpileup",
				#DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/ENA-sorted.vcf.gz",
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf.gz",
        		expand(DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/{contig}/{contig}-{sample}-aligned.fasta", sample = SAMPLES, contig = CONTIGS)


rule bowtie2virus:
		input:
			read1 = DATADir + "/bowtie2/varroa/reads_unmapped/{sample}.1",
			read2 = DATADir + "/bowtie2/varroa/reads_unmapped/{sample}.2",
		threads: 12
		output:
			alignment = temp(DATADir + "/bowtie2/virus/ENA-renamed/alignments/{sample}.bam"),
			bam = DATADir + "/bowtie2/virus/ENA-renamed/{sample}.bam"
		shell:
		 	"""
		 	bowtie2 -p {threads} --very-sensitive-local -x {renamedIndex} -1 {input.read1} -2 {input.read2} | samtools view -Su - | samtools sort - -m 20G -T {SCRATCH}/bowtie2/{wildcards.sample} -o - | samtools rmdup - {output.alignment}
		 	samtools index {output.alignment}
		 	"""


rule stats:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/{sample}.bam"
		output:
			flag = DATADir + "/bowtie2/virus/ENA-renamed/stats/flagstat/{sample}.txt",
			idx = DATADir + "/bowtie2/virus/ENA-renamed/stats/idxstats/{sample}.txt"
		threads: 2
		shell:
			"""
			samtools flagstat {input} > {output.flag}
			samtools idxstats {input} > {output.idx}
			"""


rule variant:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/alignments/{sample}.bam"
		threads: 8
		output:
			DATADir + "/bowtie2/virus/ENA-renamed/variantbam/{sample}.bam"
		shell:
			"""
			variant {input} -m 200 -b -o {output}
			samtools index {output}
			"""


rule mpileup:
		input:
			expand(DATADir + "/bowtie2/virus/ENA-renamed/variantbam/{sample}.bam", sample = SAMPLES)
		threads: 8
		output:
			pileup = DATADir + "/bowtie2/virus/ENA-renamed/varscan/ENA-sorted.mpileup",
			txt = DATADir + "/bowtie2/virus/ENA-renamed/varscan/feed.varscan.txt"
		shell:
			"""
			echo {input} | tr " " "\n" > /flash/MikheyevU/Nonno/temp/data/bowtie2/virus/ENA-renamed/varscan/bamlist.txt
			samtools mpileup -f {VIRUSrenamed} --bam-list /flash/MikheyevU/Nonno/temp/data/bowtie2/virus/ENA-renamed/varscan/bamlist.txt > {output.pileup}
			for i in {input}; do basename $i .bam; done > {output.txt}
			"""


rule varscan:
		input:
			mpileup = DATADir + "/bowtie2/virus/ENA-renamed/varscan/ENA-sorted.mpileup",
			txt = DATADir + "/bowtie2/virus/ENA-renamed/varscan/feed.varscan.txt"
		threads: 4
		output:
			vcf = DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/ENA-sorted.vcf",
			gz = DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/ENA-sorted.vcf.gz"
		shell:
			"""
			java -jar VarScan.jar mpileup2cns {input.mpileup} --vcf-sample-list {input.txt} --output-vcf 1 > {output.vcf} 
			bgzip {output.vcf}
			"""


rule tabix:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/ENA-sorted.vcf.gz"
		output:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/ENA-sorted.vcf.gz.tbi"
		shell:
			"""
			tabix {input}
			"""


rule indel:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/ENA-sorted.vcf.gz"
		threads: 6
		output:
			prefix = DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel",
			vcf = DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf"
		shell:
			"""
			vcftools --gzvcf {input} --remove-indels --recode --recode-INFO-all --out {output.prefix}
			"""

rule bgzip:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf"
		output:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf.gz"
		shell:
			"""
			bgzip {input}
			"""

rule tabix2:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf.gz"
		output:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf.gz.tbi"
		shell:
			"""
			tabix {input}
			"""

rule decompose:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel.recode.vcf.gz"
		threads: 6
		output:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel-decomp.vcf"
		shell:
			"""
			bcftools norm -a -O z --output {output} {input}
			"""


rule bgzip2:
		input:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel-decomp.vcf"
		output:
			DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel-decomp.vcf.gz"
		shell:
			"""
			bgzip {input}
			"""


rule consensus_mtDNA:
        input:
        	DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/indel-decomp.vcf.gz"
        output:
        	DATADir + "/bowtie2/virus/ENA-renamed/varscan/fasta/{contig}/{contig}-{sample}-aligned.fasta"
        threads: 6
        shell:
        	"""
        	samtools faidx {VIRUSrenamed} {wildcards.contig} | bcftools consensus -s
        	{wildcards.sample} -I {input} > {output}
        	if [ ! -f {output}]
        	then
        		echo "No output exists, generating dummy output".
        		samtools faidx {VIRUSrenamed} {wildcards.contig} > {output}
        	fi
        	ID={wildcards.sample}
        	VIRUS={wildcards.contig}
        	sed -i -e "1s/.*/>$ID-$VIRUS/" {output}
        	"""