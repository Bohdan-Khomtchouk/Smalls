
library(stringi)
library(stringr)
 z_array<-function(s) {
#""" Use Z algorithm (Theorem 1.4.1 from Gusfield's book 'Algorithms on Strings, Trees and Sequences: Computer Science and Computational Biology') to preprocess s """
                     stopifnot(length(s) > 1)
                     z = [length(s)] + [0] * (length(s)-1)
                     
                     # Initial comparison of s[1:] with prefix
                     for(i in i:length(s)){
                      if (s[i] == s[i-1]){
                       #increment by one
                        z[1] =z[1]+1
                       }
                     }
                     else{
                        break
                     }
                     r=0
                     l=0
                     if(z[1] > 0){
                         l=1
                        r=z[1]
                     }
                     for (k in range(2, length(s))){
                        stopifnot( z[k] == 0)
                     }
                     if (k > r){
                     # Case 1
                     for(i in range(k, length(s)){
                        if (s[i] == s[i-k]){
                          z[k] += 1
                     }
                     else{
                        break
                     }
                       r= k+z[k]-1
                       l=k
                     else{
                     # Case 2
                     # Calculate length of beta
                     nbeta = r - k + 1
                     zkp = z[k - l]
                     if (nbeta > zkp){
                     # Case 2a: Zkp wins
                     z[k] = zkp
                     }
                     }
                     else{
                     # Case 2b: Compare characters just past r
                     nmatch = 0
                     for(i in r+1:length(s)){
                     if (s[i] == s[i - k]){
                     nmatch += 1
                     }
                     else{
                     break
                     }
                     }
                     l = k
                     r=r+nmatch
                     z[k] = r - k + 1
                     return(z)
                     }
                     
                      n_array<-function(s){
                     #""" Compile the N array (Gusfield theorem 2.2.2) from the Z array """
                     return (z_arrays[::-1])[::-1])
                      }
                     
                      big_l_prime_array<-function(p, n){
                     #""" Compile L' array (Gusfield theorem 2.2.2) using p and N array.
                    # L'[i] = largest index j less than n such that N[j] = |P[i:]| """
                     lp = [0] * length(p)
                     for(j in range(length(p)-1)){
                     i = length(p) - n[j]
                     if (i < length(p)){
                     lp[i] = j + 1
                     }
                     }
                     return (lp)
                      }
                     
                      big_l_array<-function(p, lp){
                     #""" Compile L array (Gusfield theorem 2.2.2) using p and L' array.
                     #L[i] = largest index j less than n such that N[j] >= |P[i:]| """
                     l = l[1 ] * length(p)
                     l[1] = lp[1]
                     for (i in 2:length(p)){
                     l[i] = max(l[i-1], lp[i])
                     return(l)
                     }
                     }
                     
                      small_l_prime_array<-function(n){
                     #""" Compile lp' array (Gusfield theorem 2.2.4) using N array. """
                     small_lp = [0] * length(n)
                     for (i in range(length(n))){
                     if(n[i] == i+1){  # prefix matching a suffix
                     small_lp[length(n)-i-1] = i+1
                     }
                     }
                     for (i in length(n)-2, -1){  # "smear" them out to the left
                     if (small_lp[i] == 0){
                     small_lp[i] = small_lp[i+1]
                     }
                     return(small_lp)
                     }
                     }
                      good_suffix_table<-function(p){
                     #""" Return tables needed to apply good suffix rule. """
                     n = n_array(p)
                     lp = big_l_prime_array(p, n)
                     return (lp, big_l_array(p, lp), small_l_prime_array(n))
                      }
                     
                      good_suffix_mismatch<-function(i, big_l_prime, small_l_prime){
                    # """ Given a mismatch at offset i, and given L/L' and l' arrays,
                    # return amount to shift as determined by good suffix rule. """
                     lengthObj = length(big_l_prime)
                     stopifnot(i < length)
                     if(i == lengthObj - 1){
                        return (0)
                     }
                     i =i+ 1  # i points to leftmost matching position of P
                     if (big_l_prime[i] > 0){
                     return (lengthObj - big_l_prime[i])
                     return (lengthObj - small_l_prime[i])
                     }
                      }
                     
                      good_suffix_match<-function(small_l_prime){
                     #""" Given a full match of P to T, return amount to shift as
                     #determined by good suffix rule. """
                     return(len(small_l_prime) - small_l_prime[1])
                      }
                     
                      dense_bad_char_tab<-function(p, amap){
                     #""" Given pattern string and list with ordered alphabet characters, create
                     #and return a dense bad character table.  Table is indexed by offset
                     #then by character. """
                     tab = []
                     nxt = [0] * length(amap)
                     for (i in 0:length(p)){
                     c = p[i]
                     stopifnot( c in amap)
                     tab.append(nxt[:])
                     nxt[amap[c]] = i+1
                     return (tab)
                     }
                      }
}

                     class BoyerMoore(object):
                     """ Encapsulates pattern and associated Boyer-Moore preprocessing. """
                     
                      __init__<-function( p, alphabet='ACGT'){
                     p = p
                     alphabet = alphabet
                      
                     # Create map from alphabet characters to integers
                    amap = {}
                     for(i in range(length(alphabet))){
                       amap[alphabet[i]] = i
                     }
                     # Make bad character rule table
                     bad_char = dense_bad_char_tab(p, amap)
                     
                     # Create good suffix rule table
                     _, big_l, small_l_prime = good_suffix_table(p)
                      }
                      bad_character_rule<-function( i, c){
                     #""" Return # skips given by bad character rule at offset i """
                    stopifnot( c in amap)
                     ci = amap[c]
                     stopifnot( i > (bad_char[i][ci]-1))
                     return (i - (bad_char[i][ci]-1))
                      }
                      good_suffix_rule<-function(i){
                   #  """ Given a mismatch at offset i, return amount to shift
                    # as determined by (weak) good suffix rule. """
                     lengthObj = length(big_l)
                     stopifnot( i < lengthObj)
                     if(i == lengthObj - 1){
                     return(0)
                     }
                     i = i+1  # i points to leftmost matching position of P
                     if (big_l[i] > 0){
                     return(lengthObj - big_l[i])
                     return(lengthObj - small_l_prime[i])
                     }
                      }
                      match_skip<-function(){
                     #""" Return amount to shift in case where P matches T """
                     return(length(small_l_prime) - small_l_prime[1])
                     }
                     }
 
                    