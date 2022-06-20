# barseq-obilab

# how to run 
Step1 :
python barseq.py -i ./OBI_fertility_screen_2/allFastaFiles -b barcode_gene_file.csv -r result_170620_pool2

arg1: barseq.py is script 
arg2 with -i: path for all fastq files folder for diffrent smaples 
arg3 with -b : barcode to gene mapping (for eaxmple see: barcode_gene_file.csv)
arg4 with -r : result folder 

Step 2: remove all zero counts 





### md5sum	to test all fastq files are correctly downloaded. 
for i in `ls *md5`;do md5sum -c $i; done





