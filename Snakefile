configfile: "config.yaml"

SAMPLES = ["A001_91"]
fastqPath = config["fastqPath"]
#fastqPath = "/mnt/lustre/groups/herv_project/Brainbank/rnaseq/run1/"
outPath = config["outPath"]

rule all:
    input:
        #expand(outPath + "bwa_mem/{sample}.sam", sample=SAMPLES),
        #expand(outPath + "sorted/{sample}.bam", sample=SAMPLES),
        expand(outPath + "counts/{sample}.cntfile", sample=SAMPLES),
        expand(outPath  + "bwa_mem/{sample}.sam", sample=SAMPLES),
        expand(outPath  + "star/{sample}.Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        expand(outPath  + "htseq/{sample}.htseq.cnt", sample=SAMPLES)
rule interleave:
    input:
        fastqPath + "{sample}_R1_001.fastq.gz",
        fastqPath + "{sample}_R2_001.fastq.gz"
    output: 
        outPath  +  "interleaved/{sample}.fq.gz"
        #"/mnt/lustre/groups/herv_project/snakemake/ervmap/interleaved/{sample}.fastq"
    conda:
        "envs/bbmap.yaml"
    resources:
        mem_mb=16000
    shell:
        """
        reformat.sh in1={input[0]} in2={input[1]} out={output}
        """

rule clipAdapters:
     input:
         outPath  +  "interleaved/{sample}.fq.gz"
     output:
         outPath  +  "cleaned/{sample}.fq.gz"
     conda:
        "envs/bbmap.yaml"
     resources:
        mem_mb=16000
     shell:
     	 """
         bbduk.sh -Xmx4g in={input} out={output} ref={config[resources]}original.adapters threads=2 ktrim=r k=21 mink=11 minlen=25 hdist=1 tbo tpe 
         """

rule qualityFilter:
      input:
         outPath  +  "cleaned/{sample}.fq.gz"
      output:
         outPath  +  "qual/{sample}.fq.gz" 
      conda:
        "envs/bbmap.yaml"
      resources:
        mem_mb=16000
      shell:
         """
         bbduk.sh -Xmx4g in={input} out={output} threads=2 qtrim=rl trimq=20 minlen=25 
         """

rule bwa_mem_align:
      input:
        outPath  +  "qual/{sample}.fq.gz"
      output:
        outPath + "bwa_mem/{sample}.sam"
      conda:
        "envs/bwa.yaml"
      shell:
        """
        bwa mem -t 2 -p {config[ref_hg38]} {input} > {output}
        """

rule parse_bam_ERVmap:
      input:
        outPath + "bwa_mem/{sample}.sam"
      output:
        outPath + "ervMap/{sample}.bam"
      conda:
         "envs/samtools.yaml"
      shell:
        """
        samtools view -Sh -F4 {input} > {output}.tmp1
        perl {config[scripts]}/parse_bam.pl {output}.tmp1 > {output}.tmp2
        rm {output}.tmp1
        samtools view -bSh {output}.tmp2 > {output}.tmp3
        rm {output}.tmp2
        samtools sort -@ 2 {output}.tmp3 -o {output} > {output} 
        samtools index {output}
        rm {output}.tmp3
        """


rule ERVmap_counts:
      input:
        outPath + "ervMap/{sample}.bam"
      output:
        outPath + "counts/{sample}.cntfile"
      conda:
        "envs/bedtools.yaml"  
      shell:
        """
        bedtools multicov -bams {input} -bed {config[ERVmapBed]} > {output} 
        """

rule STAR:
      input:
         outPath  +  "qual/{sample}.fq.gz" 
      output:
         outPath + "star/{sample}.Aligned.sortedByCoord.out.bam",
         outPath + "star/{sample}.Aligned.sortedByCoord.out.bam.bai"
      conda:
        "envs/star.yaml"
      threads: 8
      resources:
        mem_mb=80000
      shell:
         """
         STAR --runThreadN 6 --genomeDir {config[genomeDir]} --sjdbGTFfile {config[gtf]} --sjdbOverhang 149 --readFilesIn {input} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix {wildcards.sample}. --outStd BAM_SortedByCoordinate > {output[0]}
         samtools index {output[0]}
         """


rule HTseqCount:
      input:
         outPath + "star/{sample}.Aligned.sortedByCoord.out.bam"
      output:
         outPath + "htseq/{sample}.htseq.cnt"
      conda:
         "envs/htseq.yaml"
      shell:
         """
         samtools view {input} | htseq-count -f sam -s no - {config[gtf]} > {output}
         """
