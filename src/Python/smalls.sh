# Copyright (C) 2015 Bohdan Khomtchouk and Derek Van Booven

# This shell script file is part of Smalls.

# Smalls is a Python and shell scripting software program that efficiently aligns miRNA seed sequences to three-prime and five-prime untranslated regions in DNA (3'-UTR and 5'-UTR gene regions).
# Smalls' efficiency and flexibility stem from its utilization of advanced string algorithms and user-specified options for approximate matching (e.g., 2 mismatches).
# miRNA binding at UTR regions has been shown to influence biological processes as diverse as embryonic development or disease progression.  It is, therefore, imperative
# to be able to computationally evaluate a candidate list of miRNAs in relation to their potential UTR binding targets.  By evaluating which specific miRNAs align most often
# in the genome (as quantified by raw count of successful alignment at a given approximate mismatch rate), we can determine which specific miRNAs are most promising to 
# experimentally validate in the lab.  More importantly, we can determine the precise genomic position(s) of these miRNA-UTR alignments and, therefore, the identity of the respective gene(s)
# involved in the biological process at hand.

# Current work is underway to expand Smalls to cover other types of small RNAs in addition to miRNAs (e.g., snoRNAs, snRNAs, siRNAs, piRNAs, lncRNAs, vlincRNAs).           

# BK wishes to acknowledge financial support by the United States Department of Defense (DoD) through the National Defense Science and Engineering 
# Graduate Fellowship (NDSEG) Program: this research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense 
# Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

# This research was conducted in collaboration with and using the resources of the University of Miami Center for Computational Science.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ------------------------------------------------------------------------------------

# N.B.: this shell script serves merely to demonstrate the procedure of extraction of all UTRs from a given genome (e.g., human), as well as the procedure of data mining for mature miRNA seed sequences (from all known organisms).
# As such, this general shell script is easily customizable and expandable.  Its output (UTRome.fa and mir_seeds.fa) are to be preprocessed and used as input into the main Python scripts of Smalls. 

awk '{ print $1 "\t" $4 "\t" $5 "\t" $10 }' UTR.gtf > genes.bed		#UTR.gtf has already been pre-made by the user beforehand (e.g., all UTR regions for human).  Here we preprocess this general transfer format (gtf) file. 
sed -i 's/"//g' genes.bed		#preprocessing step
sed -i 's/;//g' genes.bed		#preprocessing step
bedtools getfasta -fi /hihg/ref/genomes/bwa_hg19/human_g1k_v37.fasta -bed genes.bed -split -name -fo UTRome.fa		#uses publicly available human_g1k_v37.fasta. For more info on bedtools command usage, please read: http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html
wget ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz		#download all mature miRNA sequences from miRBase database
gunzip mature.fa.gz		#unzip this file
sed -i 's/U/T/g' mature.fa	#substitute uracils for thymines
/share/apps/hihg/src/fastx_toolkit/fastx_trimmer -Q33 -f 2 -l 9 -i mature.fa -o mir_seeds.fa		#generate fasta file of 8-mer miRNA seed sequences using FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/index.html)