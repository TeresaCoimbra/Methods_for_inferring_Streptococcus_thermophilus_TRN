# *Streptococcus thermophilus* - read alignment to reference genome

Pulling Bowtie2:
```
$podman pull biocontainers/bowtie2: 2.4.3--py38h72fc82f_0
```

Building of a bowtie samtools image:
```
podman build -f Dockerfile -t bowtie_samtools 
```

Example of command used to run the bowtie_samtools image:
```
$ podman run -it -v ./PRJNA395470/:/home/PRJNA395470:z  -v ./genome_index/:/home/genome_index:z bowtie_samtools bash
```

Script used to convert align and convert directly from sam to bam files (mapping.sh):
```
#!/usr/bin/env bash

for file in *.gz
do
        bamfile=`echo $file | sed 's/.gz/.bam/g'`
        bowtie2 -x ./../genome_index/sth_index -U $file | samtools view -Sb - > $bamfile
done
```

Command used to make the script executable:
```
$ chmod +x mapping.sh
```

In some files the name after alignment was ".fastq.bam", which was renamed into ".bam":
```
for file in *.bam; do mv "$file" "${file/.fastq.bam/.bam}"; done
```

Generating alignment statistics usign flagstat:
```
samtools flagstat <filename.bam> > <filename.txt>
```

Example of commands used to generate the statistics for the dataset PRJNA412475:
```
samtools flagstat SRR6112518_1.fastq.bam > SRR6112518_1.txt
samtools flagstat SRR6112518_2.fastq.bam > SRR6112518_2.txt
samtools flagstat SRR6112538_2.fastq.bam > SRR6112538_2.txt
samtools flagstat SRR6112538_1.fastq.bam > SRR6112538_1.txt
samtools flagstat SRR6112863_1.fastq.bam > SRR6112863_1.txt
samtools flagstat SRR6112863_2.fastq.bam > SRR6112863_2.txt
samtools flagstat SRR6112864_1.fastq.bam > SRR6112864_1.txt
samtools flagstat SRR6112864_2.fastq.bam > SRR6112864_2.txt
samtools flagstat SRR6112865_2.fastq.bam > SRR6112865_2.txt
samtools flagstat SRR6112865_1.fastq.bam > SRR6112865_1.txt
samtools flagstat SRR6112866_1.fastq.bam > SRR6112866_1.txt
samtools flagstat SRR6112866_2.fastq.bam > SRR6112866_2.txt
```


## Sorting the reads per genomic coordinates:
Another script was generated to sort the reads for genomic coordinates for all datasets ("sort.sh"):
```
#!/usr/bin/env bash
for file in *.bam
do
  sorted_file=`echo $file | sed 's/.bam/_sorted.bam/g'`
	stats_file=`echo $file | sed 's/.bam/.txt/g'`
	samtools sort $file -o $sorted_file
  samtools flagstat $sorted_file > $stats_file
done
```
The statistics files were compared and these files were removed since they didn't add any visible improvements. 

