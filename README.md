# RNASeqHERV
RNAseq pipeline adapted from Ashley Jones 

<div align="center">
    <img src="RNAseq.jpg" width="200px"</img> 
</div>


# Description
This pipeline aims to recreate the differential expression analysis described in https://www.nature.com/articles/s41598-021-93742-3.pdf

# Steps

1.Interleave fastqs <br>
2.Trim adapters  <br>
3. Quality Control  <br>
4.BWA MEM align  <br>
5.Parse the alignments using adapted ERVmap scripts  <br>
6.STAR align (output of Step 3)  <br>
7.HTseq count   <br>

# Dependencies
When ran on Rosalind, environment modules are used and will automatically be loaded (--use-envmodules should be used when invoking snakemake). All other dependencies (e.g. reference genomes, HERV locations annotation, etc.) are in hervk_group folders and you will need permissions to these folders.

# Execution on Rosalind

Provided that the targets are specified in rule all, run the pipeline with:

```
nohup snakemake --use-conda --use-envmodules --cores <HOWEVER_MANY_SAMPLES> --cluster "sbatch -p brc --mem-per-cpu=35G" 
```

## Example of targets requested in rule 'all'
 ```
rule all:
    input:
        expand(outPath + "counts/{sample}/{sample_s}.cntfile", zip, sample=config["SAMPLES"], sample_s=config["SAMPLES_WITH_S_NUM"]), 
        expand(outPath  + "htseq/{sample}/{sample_s}.htseq.cnt", zip, sample=config["SAMPLES"], sample_s=config["SAMPLES_WITH_S_NUM"])
  ```
  This will produce counts for ERVs (.cntfile) and counts for cellular data (.htseq.cnt) for each sample
  
  Samples and the subfolders for each sample are defined in config.yaml. 
  
  # After running the pipeline - merge counts for ERVs and for cellular
  
  Step1: <br>
  copy to one directory using python script, remove with sed these extra few lines in the cellular counts files. Change the paths to reflect your own: <br>
  
  ```
import glob
import re
import os

os.system("mkdir -p output/erv")
os.system("mkdir -p output/cellular")

for file in /mnt/lustre/groups/projectmine/RNAseqRunningDir/output/cellular/*.htseq.cnt; do sed -i '/__/d' $file; done

for filepath in glob.iglob('/mnt/lustre/groups/projectmine/RNAseqRun*hg18/htseq/*/*.htseq.cnt'):
      os.system("cp %s output/cellular/" %(filepath))

for filepath in glob.iglob('/mnt/lustre/groups/projectmine/RNAseqRun*hg18/counts/*/*.cntfile'):
      os.system("cp %s output/erv/" %(filepath))
```
 
 Step2: <br>
 Go to directory cellular. If the file has two endings, remove the last one.  E.g. my files are A158_14_S29.htseq.cnt, etc, before the merge script get rid of the last ending: .cnt
 
 ```
 cd output/cellular
 for file in *.cnt; do     mv "$file" "${file%%.cnt}"; done
 ```
 
  Step2:<br>
  Merge using ERVmap scripts <br>
  ```
  perl scripts/merge.pl 3 6 cntfile ./output/erv > ./output/erv/merged_erv.txt
  perl scripts/merge.pl 0 1 htseq.cnt ./output/cellular > ./output/cellular/merged_cellular.txt
```
