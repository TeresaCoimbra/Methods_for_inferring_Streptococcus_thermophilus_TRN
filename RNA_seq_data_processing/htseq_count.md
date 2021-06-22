
# Using HTSeq count - count the overlapp of reads with genes

Pulling the container:
```
podman pull quay.io/biocontainers/htseq:0.13.5--py39h70b41aa_1
```

HTseq script used for counting (counting.sh):
```
#!/usr/bin/env bash

for file in ./*.bam
do
        resfile=`echo $file | sed 's/.bam/_count.txt/g'`
        htseq-count -f bam -s no -t gene -i locus_tag $file ./genome_cds.gff >> $resfile
done
```
The feature type selected was 'gene' and the GFF attribute used as feature ID was the locus tag. 

Making the file executable:
```
$chmod +x counting.sh
```

Running the HTSeq container and the script for all datasets: 
```
$podman run -v <dir>/:/data:z -it quay.io/biocontainers/htseq:0.13.5--py39h70b41aa_1 bash
```

Example for PRJNA395470:
```
$podman run -v ./PRJNA395470/aligned_files/:/data:z -it quay.io/biocontainers/htseq:0.13.5--py39h70b41aa_1 bash
#./counting.sh
```

