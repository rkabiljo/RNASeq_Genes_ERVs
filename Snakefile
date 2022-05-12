configfile: "config.yaml"

#SAMPLES = ["A001_91"]
fastqPath = config["fastqPath"]
#fastqPath = "/mnt/lustre/groups/herv_project/Brainbank/rnaseq/run1/"
outPath = config["outPath"]

rule all:
    input:
        #expand(outPath + "bwa_mem/{sample}.sam", sample=SAMPLES),
        #expand(outPath + "sorted/{sample}.bam", sample=SAMPLES),
        expand(outPath + "counts/{sample}/{sample_s}.cntfile", zip, sample=config["SAMPLES"], sample_s=config["SAMPLES_WITH_S_NUM"]), 
        expand(outPath  + "bwa_mem/{sample}/{sample_s}.sam", zip, sample=config["SAMPLES"], sample_s=config["SAMPLES_WITH_S_NUM"]),
        expand(outPath  + "star/{sample}/{sample_s}.Aligned.sortedByCoord.out.bam", zip, sample=config["SAMPLES"], sample_s=config["SAMPLES_WITH_S_NUM"]),
        expand(outPath  + "htseq/{sample}/{sample_s}.htseq.cnt", zip, sample=config["SAMPLES"], sample_s=config["SAMPLES_WITH_S_NUM"])
rule interleave:
    input:
        fastqPath + "Sample_{sample}/{sample_s}_R1_001.fastq.gz",
        fastqPath + "Sample_{sample}/{sample_s}_R2_001.fastq.gz"
    output: 
        outPath  +  "interleaved/{sample}/{sample_s}.fq.gz"
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
         outPath  +  "interleaved/{sample}/{sample_s}.fq.gz"
     output:
         outPath  +  "cleaned/{sample}/{sample_s}.fq.gz"
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
         outPath  +  "cleaned/{sample}/{sample_s}.fq.gz"
      output:
         outPath  +  "qual/{sample}/{sample_s}.fq.gz" 
      shell:
         """
         /scratch/users/k1929424/Miniconda3/pkgs/bbmap-38.18-0/bin/bbduk.sh -Xmx4g in={input} out={output} threads=2 qtrim=rl trimq=20 minlen=25 
         """

rule bwa_mem_align:
      input:
        outPath  +  "qual/{sample}/{sample_s}.fq.gz"
      output:
        outPath + "bwa_mem/{sample}/{sample_s}.sam"
      envmodules:
        "apps/bwa/0.7.17-singularity"
      shell:
        """
        bwa mem -t 2 -p {config[ref_hg38]} {input} > {output}
        """

rule parse_bam_ERVmap:
      input:
        outPath + "bwa_mem/{sample}/{sample_s}.sam"
      output:
        outPath + "bwa_mem/{sample}/{sample_s}.bam"
      envmodules:
        "apps/samtools/0.1.19-singularity"
      shell:
        """
        samtools view -Sh -F4 {input} > {output}.tmp1
        perl {config[scripts]}/parse_bam.pl {output}.tmp1 > {output}.tmp2
        rm {output}.tmp1
        samtools view -bSh {output}.tmp2 > {output}
        rm {output}.tmp2
        """

rule sort_index:
      input:
        outPath + "bwa_mem/{sample}/{sample_s}.bam"
      output:
        outPath + "sorted/{sample}/{sample_s}.bam"
      envmodules:
        "apps/samtools/0.1.19-singularity"
      shell:
        """
        samtools sort -@ 2 {input} -o {output} > {output} 
        samtools index {output}
        """

rule ERVmap_counts:
      input:
        outPath + "sorted/{sample}/{sample_s}.bam"
      output:
        outPath + "counts/{sample}/{sample_s}.cntfile"
      shell:
        """
        bedtools multicov -bams {input} -bed {config[ERVmapBed]} > {output} 
        """

rule STAR:
      input:
         outPath  +  "qual/{sample}/{sample_s}.fq.gz" 
      output:
         outPath + "star/{sample}/{sample_s}.Aligned.sortedByCoord.out.bam",
         outPath + "star/{sample}/{sample_s}.Aligned.sortedByCoord.out.bam.bai"
      envmodules:
         "apps/star/2.7.3a",
         "apps/samtools/0.1.19-singularity"
      shell:
         """
         STAR --runThreadN 6 --genomeDir {config[genomeDir]} --sjdbGTFfile {config[gtf]} --sjdbOverhang 149 --readFilesIn {input} --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix {wildcards.sample}. --outStd BAM_SortedByCoordinate > {output[0]}
         samtools index {output[0]}
         """


rule HTseqCount:
      input:
         outPath + "star/{sample}/{sample_s}.Aligned.sortedByCoord.out.bam"
      output:
         outPath + "htseq/{sample}/{sample_s}.htseq.cnt"
      envmodules:
        "apps/samtools/0.1.19-singularity"
      shell:
         """
         pip install --upgrade numpy
         #samtools view {input} | python -m HTSeq.scripts.count -f sam -s no - {config[gtf]} > {output}
         python -m HTSeq.scripts.count -f bam -s no {input} {config[gtf]} > {output}
         """
