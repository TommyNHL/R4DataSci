fastqc *.gz -t 8

fastp -i *.gz -t 8

kallisto index -i Mus_musculus.GRCm39.cdna.all.index Mus_musculus.GRCm39.cdna.all.fa

kallisto quant -i Mus_musculus.GRCm39.cdna.all.index -o M0000 -t 8 M0000-R1.fastq.gz M0000-R2.fastq.gz &> M0000.log

multiqc -d .

echo "Finished"