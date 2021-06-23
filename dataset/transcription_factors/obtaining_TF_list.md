# *Streptococcus thermophilus* - List of Transcription Factors

Transcription factors names were obtained by searching the gff file that was used for the alignment step:
```
$grep -wi transcription genome_cds.gff | cat > tfs.txt
```
The resulting file was processed to retrieve only the relevant information. 

Other transcription factors were found by searching on UniProt and NCBI protein databases. 
To obtain the correct locus tag notation (new locus tag) information was retrieved from: 
https://regprecise.lbl.gov/genome.jsp?genome_id=36

Is is important to mention that all transcription factors are anotated for the *Streptococcus thermophilus* strain LMD-9, which is the one used in the reference genome, as well as the one used on the largest dataset PRJNA395470.

In this folder it is possible to find the resulting file: "transcription_factors.csv".
