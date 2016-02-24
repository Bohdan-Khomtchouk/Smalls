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
	
	import string
    
    def z_array(s):
        """ Use Z algorithm (Theorem 1.4.1 from Gusfield's book 'Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology') to preprocess s """
        assert len(s) > 1
        z = [len(s)] + [0] * (len(s)-1)
    
        # Initial comparison of s[1:] with prefix
        for i in range(1, len(s)):
            if s[i] == s[i-1]:
                z[1] += 1
            else:
                break
    
        r, l = 0, 0
        if z[1] > 0:
            r, l = z[1], 1
    
        for k in range(2, len(s)):
            assert z[k] == 0
            if k > r:
                # Case 1
                for i in range(k, len(s)):
                    if s[i] == s[i-k]:
                        z[k] += 1
                    else:
                        break
                r, l = k + z[k] - 1, k
            else:
                # Case 2
                # Calculate length of beta
                nbeta = r - k + 1
                zkp = z[k - l]
                if nbeta > zkp:
                    # Case 2a: Zkp wins
                    z[k] = zkp
                else:
                    # Case 2b: Compare characters just past r
                    nmatch = 0
                    for i in range(r+1, len(s)):
                        if s[i] == s[i - k]:
                            nmatch += 1
                        else:
                            break
                    l, r = k, r + nmatch
                    z[k] = r - k + 1
        return z
    
    
    def n_array(s):
        """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
        return z_array(s[::-1])[::-1]
    
    
    def big_l_prime_array(p, n):
        """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
            L'[i] = largest index j less than n such that N[j] = |P[i:]| """
        lp = [0] * len(p)
        for j in range(len(p)-1):
            i = len(p) - n[j]
            if i < len(p):
                lp[i] = j + 1
    
        return lp
    
    
    def big_l_array(p, lp):
        """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
            L[i] = largest index j less than n such that N[j] >= |P[i:]| """
        l = [0] * len(p)
        l[1] = lp[1]
        for i in range(2, len(p)):
            l[i] = max(l[i-1], lp[i])
        return l
    
    
    def small_l_prime_array(n):
        """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
        small_lp = [0] * len(n)
        for i in range(len(n)):
            if n[i] == i+1:  # prefix matching a suffix
                small_lp[len(n)-i-1] = i+1
        for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
            if small_lp[i] == 0:
                small_lp[i] = small_lp[i+1]
        return small_lp
    
    
    def good_suffix_table(p):
        """ Return tables needed to apply good suffix rule. """
        n = n_array(p)
        lp = big_l_prime_array(p, n)
        return lp, big_l_array(p, lp), small_l_prime_array(n)
    
    
    def good_suffix_mismatch(i, big_l_prime, small_l_prime):
        """ Given a mismatch at offset i, and given L/L' and l' arrays,
            return amount to shift as determined by good suffix rule. """
        length = len(big_l_prime)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if big_l_prime[i] > 0:
            return length - big_l_prime[i]
        return length - small_l_prime[i]
    
    
    def good_suffix_match(small_l_prime):
        """ Given a full match of P to T, return amount to shift as
            determined by good suffix rule. """
        return len(small_l_prime) - small_l_prime[1]
    
    
    def dense_bad_char_tab(p, amap):
        """ Given pattern string and list with ordered alphabet characters, create
            and return a dense bad character table.  Table is indexed by offset
            then by character. """
        tab = []
        nxt = [0] * len(amap)
        for i in range(0, len(p)):
            c = p[i]
            assert c in amap
            tab.append(nxt[:])
            nxt[amap[c]] = i+1
        return tab
    
    
    class BoyerMoore(object):
        """ Encapsulates pattern and associated Boyer-Moore preprocessing. """
    
        def __init__(self, p, alphabet='ACGT'):
            self.p = p
            self.alphabet = alphabet
    
            # Create map from alphabet characters to integers
            self.amap = {}
            for i in range(len(self.alphabet)):
                self.amap[self.alphabet[i]] = i
    
            # Make bad character rule table
            self.bad_char = dense_bad_char_tab(p, self.amap)
    
            # Create good suffix rule table
            _, self.big_l, self.small_l_prime = good_suffix_table(p)
    
        def bad_character_rule(self, i, c):
            """ Return # skips given by bad character rule at offset i """
            assert c in self.amap
            ci = self.amap[c]
            assert i > (self.bad_char[i][ci]-1)
            return i - (self.bad_char[i][ci]-1)
    
        def good_suffix_rule(self, i):
            """ Given a mismatch at offset i, return amount to shift
                as determined by (weak) good suffix rule. """
            length = len(self.big_l)
            assert i < length
            if i == length - 1:
                return 0
            i += 1  # i points to leftmost matching position of P
            if self.big_l[i] > 0:
                return length - self.big_l[i]
            return length - self.small_l_prime[i]
    
        def match_skip(self):
            """ Return amount to shift in case where P matches T """
            return len(self.small_l_prime) - self.small_l_prime[1]