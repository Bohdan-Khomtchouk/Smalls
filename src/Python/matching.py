#!/usr/bin/env python

# This file is currently under development for project "Smalls" by Bohdan Khomtchouk and Derek Van Booven.

# Smalls is a Python and shell scripting software program that efficiently aligns miRNA seed sequences to three-prime and five-prime untranslated regions in DNA (3'-UTR and 5'-UTR gene regions).
# Smalls' efficiency and flexibility stem from its utilization of advanced string algorithms and user-specified options for approximate matching (e.g., 2 mismatches).
# miRNA binding at UTR regions has been shown to influence biological processes as diverse as embryonic development or disease progression.  It is, therefore, imperative
# to be able to computationally evaluate a candidate list of miRNAs in relation to their potential UTR binding targets.  By evaluating which specific miRNAs align most often
# in the genome (as quantified by raw count of successful alignment at a given approximate mismatch rate), we can determine which specific miRNAs are most promising to 
# experimentally validate in the lab.  More importantly, we can determine the precise genomic position(s) of these miRNA-UTR alignments and, therefore, the identity of the respective gene(s)
# involved in the biological process at hand.  Current work is underway to expand Smalls to cover other types of small RNAs in addition to miRNAs (e.g., snoRNAs, snRNAs, siRNAs, piRNAs, lncRNAs, vlincRNAs). 

# ------------------------------------------------------------------------------------	

"""Algorithms for DNA Sequencing - Johns Hopkins University"""

		import preliminaries #imports preliminaries.py file

		def readGenome():
        	genome = ''
        	with open('filename', 'r') as f:
            	for line in f:
                	# ignore header line with genome information
                	if not line[0] == '>':
                    	genome += line.rstrip()
        	return genome
		
		def boyer_moore(p, p_bm, t):
        	""" Do Boyer-Moore matching """
        	i = 0
        	occurrences = []
        	while i < len(t) - len(p) + 1:
            	shift = 1
            	mismatched = False
            	for j in range(len(p)-1, -1, -1):
                	if p[j] != t[i+j]:
                    	skip_bc = p_bm.bad_character_rule(j, t[i+j])
                    	skip_gs = p_bm.good_suffix_rule(j)
                    	shift = max(shift, skip_bc, skip_gs)
                    	mismatched = True
                    	break
            	if not mismatched:
                	occurrences.append(i)
                	skip_gs = p_bm.match_skip()
                	shift = max(shift, skip_gs)
            	i += shift
        	return occurrences
    

		def approximate_match(p, t, n):
        	segment_length = int(round(len(p) / (n+1)))
        	all_matches = set()
        	for i in range(n+1):
            	start = i*segment_length
            	end = min((i+1)*segment_length, len(p))
            	p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
            	matches = boyer_moore(p[start:end], p_bm, t)
    
            	# Extend matching segments to see if whole p matches
            	for m in matches:
                	if m < start or m-start+len(p) > len(t):
                    	continue
    
                	mismatches = 0
    
                	for j in range(0, start):
                    	if not p[j] == t[m-start+j]:
                        	mismatches += 1
                        	if mismatches > n:
                            	break
                	for j in range(end, len(p)):
                    	if not p[j] == t[m-start+j]:
                        	mismatches += 1
                        	if mismatches > n:
                            	break
    
                	if mismatches <= n:
                    	all_matches.add(m - start)
        	return list(all_matches)


# ------------------------------------------------------------------------------------