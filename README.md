# Smalls

## About

Smalls is destined to become a noncoding RNA R package (i.e., I'm working on it!).  Currently, Smalls is a Python/R/shell scripting software program that efficiently aligns miRNA seed sequences to three-prime and five-prime untranslated regions in DNA (3'-UTR and 5'-UTR gene regions).
Smalls' efficiency and flexibility stem from its utilization of advanced string algorithms and user-specified options for approximate matching (e.g., 2 mismatches).

miRNA binding at UTR regions has been shown to influence biological processes as diverse as embryonic development or disease progression.  It is, therefore, imperative
to be able to computationally evaluate a candidate list of miRNAs in relation to their potential UTR binding targets.  By evaluating which specific miRNAs align most often
in the genome (as quantified by raw count of successful alignment at a given approximate mismatch rate), we can determine which specific miRNAs are most promising to 
experimentally validate in the lab.  More importantly, we can determine the precise genomic position(s) of these miRNA-UTR alignments and, therefore, the identity of the respective gene(s)
involved in the biological process at hand.

Current work is underway to expand Smalls to cover other types of small RNAs in addition to miRNAs (e.g., snoRNAs, snRNAs, siRNAs, piRNAs, lncRNAs, vlincRNAs).           

BK wishes to acknowledge financial support by the United States Department of Defense (DoD) through the National Defense Science and Engineering 
Graduate Fellowship (NDSEG) Program: this research was conducted with Government support under and awarded by DoD, Army Research Office (ARO), National Defense 
Science and Engineering Graduate (NDSEG) Fellowship, 32 CFR 168a.

This research was conducted in collaboration with and using the resources of the University of Miami Center for Computational Science.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Example usage:

### R:

<img width="1067" alt="screen shot 2016-05-13 at 12 13 02 pm" src="https://cloud.githubusercontent.com/assets/9893806/15254326/40fee134-1904-11e6-9908-a5bc97be1c47.png">

### Python:

If filename is 'chr1.GRCh38.excerpt.fasta', then:

* `p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'` 47 base pairs

* `t = readGenome()` MANY base pairs

* `print(approximate_match(p, t, 2))`
* `[160162, 147558, 364263, 717706, 465647, 429299, 657496, 160729, 56922, 191452, 724927]`

Prints all starting positions where p has matched t (with 0, 1, or 2 mismatches allowed).  For sake of example, here are the first 24 base pairs of the respective matches (using IPython QtConsole 3.2.0): 

* `In [1]: list = [t[160162:160162+24], t[147558:147558+24], t[364263:364263+24], t[717706:717706+24], t[465647:465647+24], t[429299:429299+24], t[657496:657496+24], t[160729:160729+24], t[56922:56922+24], t[191452:191452+24], t[724927:724927+24]]`

* `In [2]: list`
* `Out[2]:` 
`['GGCACGGTGGCTCACGCATGTAAT',`
 `'GGCGCGGTGGCTCATGCCTGTAAT',`
 `'GGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCGCAGTGGCTCACGCCTGTAAT',`
 `'AGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCGCGGTGGCTCACACCTGTAAT',`
 `'GGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCGCGGTGGTTCACGCCTGTAAT',`
 `'GGCACGGTGGCTCACGCCTGTAAT']`

* `In [3]: len(list)`
* `Out[3]: 11`

* `In [4]: set(list) #unique 24-mers`
* `Out[4]:` 
`{'AGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCACGGTGGCTCACGCATGTAAT',`
 `'GGCACGGTGGCTCACGCCTGTAAT',`
 `'GGCGCAGTGGCTCACGCCTGTAAT',`
 `'GGCGCGGTGGCTCACACCTGTAAT',`
 `'GGCGCGGTGGCTCACGCCTGTAAT',`
 `'GGCGCGGTGGCTCATGCCTGTAAT',`
 `'GGCGCGGTGGTTCACGCCTGTAAT'}`
 
* `In [5]: len(set(list))`
* `Out[5]: 8`
