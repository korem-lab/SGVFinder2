# Tips for using ICRA on large datasets with a large reference database

## Overview

The ICRA algorithm can be broken down into four steps in the view of generating new files:

1. fastq -> bam
2. bam -> pmp
3. pmp -> jsdel + jspi
4. jsdel -> smp

*The `icra` command includes `svfinder get_sample_map` which generates smp files*

When using the `icra` command on large datasets with a large reference database, significant time and memory are required.
By seperating those steps above and utilizing parallel computing, the process can be expedited.
However, considering the CPU core count, disk storage space, and memory size of the computer or server, each step requires different setting.
The arguments `bamfol` and `pmpfol` allow files to be distributed separately.
Here is an example demonstrating the whole ICRA algorithm (*30G reference database*).

| No. | Files | Multi-threads | Memoery Usage |  Time Usage|File size (Example) |
|---|---|---|---|---|---|
|0|fastq|-|-|-|16G|
|1|fastq -> bam|Yes|-|-|8.3G|
|2|bam -> pmp|No|Low (<1G)|1.5h|11G|
|3|pmp -> jsdel + jspi|No|High (30G~100G)|4h~2days|913M (jsdel) + 11K (jspi)|
|4|jsdel -> smp|No|High (32.92G)|4min|303M|

## Step 1: BAM file

The first step involves using bowtie2 and samtools to generate bam files from the fastq files.
This can be done using bowtie2 and samtools directly or using the `generate_bam` command in `icra` with multiple threads.
With the argument `--generate_bam`, the `icra` command will only generate bam files.
The default command of the `generate_bam` is this:

```bash
bowtie2 -x ${db}/${db_name} --very-sensitive -k 20 -U ${dir}/${filename}.fastq --quiet -p ${threads} \
    | samtools view -bS -@ ${threads} - > ${bam_dir}/${filename}.bam
```

## Step 2: PMP file

The second step uses the function `sam2pmp` in the script `SGVFinder2\helpers\sam2pmp.py`.
This step cannot use multithreading but consumes minimal memory (less than 1 GB) and takes a long time (1~6 hours for 16G fastq).
Therefore, it is feasible and recommended to analyze multiple samples concurrently.
This step can be executed separately using the `bam_to_pmp` argument in the `icra` command.
With the argument `--bam_to_pmp`, the `icra` command will only convert bam files to pmp files.

## Step 3~4: ICRA

With the bam and pmp files available and their paths properly specified, the `icra` command will bypass step 1-2 and continue to generate jsdel files.
Step 3 consumes a large amount of memory (30GB-100G for 16Gfastq) and considerable time (3min-20ming per iteration for 16Gfastq, 4h-2days total) with the default argument (epsilon=1e-6,
max_iterations=100).
Step 4 is the same as `svfinder get_sample_map`, which consumes large memory but not much time. So it is not necessary to seperate step 4 from step 3.

## Conclusion

Step 1 can utilize multiple threads and/or analyze multiple samples concurrently. Step 2 can analyze multiple samples concurrently.
This makes full use of the computer or server resources and shortens the analysis time.
Step 3 consumes an uncertain amount of memory and is the bottleneck of the entire ICRA analysis process, requiring careful handling.
