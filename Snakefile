SAMPLES = ["A", "B"]


rule all:
    input:
        "plots/quals.svg"


rule bwa_map:
    input:
        genome="data/genome.fa",
        sample="data/samples/{sample}.fastq",
	index1="data/genome.fa.amb",
	index2="data/genome.fa.ann",
	index3="data/genome.fa.bwt",
	index4="data/genome.fa.fai",
	index5="data/genome.fa.pac",
	index6="data/genome.fa.sa"	
    output:
        "mapped_reads/{sample}.bam"
    conda: "environment.yaml"
    shell:
        "bwa mem {input.genome} {input.sample} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    conda: "environment.yaml"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda: "environment.yaml"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    conda: "environment.yaml"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    conda: "environment.yaml"
    script:
        "scripts/plot-quals.py"
