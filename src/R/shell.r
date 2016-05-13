system(awk '{ print $1 "\t" $4 "\t" $5 "\t" $10 }' UTR.gtf > genes.bed)		#UTR.gtf has already been pre-made by the user beforehand (e.g., all UTR regions for human).  Here we preprocess this general transfer format (gtf) file. 
sed("'-i 's/"//g' genes.bed"')		#preprocessing step
sed("-i 's/;//g' genes.bed")		#preprocessing step
system(bedtools getfasta -fi /hihg/ref/genomes/bwa_hg19/human_g1k_v37.fasta -bed genes.bed -split -name -fo UTRome.fa	)	#uses publicly available human_g1k_v37.fasta. For more info on bedtools command usage, please read: http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html

       download.file(wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz	)
        unzip(mature.fa.gz)
        system("sed -i 's/U/T/g' mature.fa")	#substitute uracils for thymines
        system( wget "/share/apps/hihg/src/fastx_toolkit/fastx_trimmer -Q33 -f 2 -l 9 -i mature.fa -o mir_seeds.fa")	

