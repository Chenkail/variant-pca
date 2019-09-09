#!/bin/bash
# Generate vcf file for each barcode
BAM=$1
samtools view -H $BAM > header.txt

barcode_file=$2

barcodes=$(<$barcode_file)
minqual=50

if [ ! -d "./OneCell" ]; then
    mkdir OneCell
fi

for cell in $barcodes
do
    echo $cell
    
    folder=./OneCell/$cell
    prefix=$folder/$cell
    
    if [ ! -d "$folder" ]; then
        mkdir $folder
    fi
    
    infodir="../Info"
    if [ ! -d "$infodir" ]; then
        mkdir $infodir
    fi
    
    samfile=$prefix.sam
    vcffile=$prefix.vcf
    vcfbest=$prefix.best.vcf
    outfile=$prefix.info.tsv

    # Make files
    cp header.txt $samfile
    samtools view $BAM | egrep $cell >> $samfile
    echo "samtools complete"
    freebayes -f ../../genome.fa $samfile > $vcffile
    echo "freebayes complete"
    bcftools filter -i "%QUAL>$minqual" $vcffile > $vcfbest
    echo "bcftools complete"
    
    # Find where table starts
    top=$(awk '/#CHROM/{ print NR; exit }' $vcfbest)
    # Remove header, reduce to important columns, and output to file
    cut -f 1,2,4,5 $vcfbest | tail -n +$top > $outfile
    cp $outfile ../Info
    echo "Info file generated"
    
    echo
    
done
