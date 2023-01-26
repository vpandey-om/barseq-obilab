# barseq-obilab

# how to run
Install: Need some Library such as screed regex.

Step1 :
python barseq.py -i data/sequence -b barcode_gene_file.csv -r testRes

arg1: barseq.py is script
arg2 with -i: path for all fastq files containing folder
arg3 with -b : barcode to gene mapping (for example see: barcode_gene_file.csv)
arg4 with -r : user based result folder (testRes); you have to remove folder if script fails.

Step 2: remove all zero counts

python remove_zero_barseq_count.py testRes/barcode_counts_table.csv Result_after_removing_zero.csv


arg1: This file will be genrated after step 1.
arg2: Desired file.
### md5sum	to test all fastq files are correctly downloaded.
for i in `ls *md5`;do md5sum -c $i; done
