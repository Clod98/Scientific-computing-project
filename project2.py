"""
Claudiu Acsinte 233207
"""
#import matplotlib.pyplot as plt
from time import time

def kmer_hist(S,k):
    """
    Args:
        S (string): the input string
        k (int): the length of the k_mers

    Returns:
        h (list): histogram list, in h[l] we have how many unique 
                k_mers have frequency l
        mfkmers (list): all the kmers having maximum frequency
    """
    #dictionary containing frequency of each kmer
    freq = {}
    
    #maximum frequency found so far
    max = 0
    
    #scanning linealry all the kmers in string S
    for i in range(len(S)-k+1):
        kmer = S[i:i+k]
        #updates frequency
        freq[kmer]=freq.get(kmer,0)+1
        #updates maximum frequency
        if freq[kmer] > max:
            max = freq[kmer]
            
    #creating the histogram list
    h=[0]*(max+1)
    for k in freq:
        h[freq[k]] = h[freq[k]]+1 
    #creating the list of most frequent kmers
    mfkmer = [kmer for kmer in freq if freq[kmer]==max]

    return h,mfkmer



def rol_hash(old, new, h , prec):
    """computes the new hash for the string S[i+1:i+k+1] 
        updating the hash value of S[i:i+k].
        The numerical basis chosen is 4 since we have 4 nucleotids 

    Args:
        old (int): the charachter to be eliminated from 
                    the old kmer, S[i] (already turned into a number)
        new (int): the charachter to be added from 
                    the new kmer, S[i+k+1] (already turned into a number)
        h (int): the old hash value
        prec (int): the value (basis)**(len(kmer)-1)

    Returns:
        int: the new hash value, we choose a big prime number for the modulus,
            a Marsenne prime number
    """
    return ((h-prec*old)*4+new)%(2**13-1)


def init_hash(s):
    """computes hash value for strings at the beginning of the procedure

    Args:
        s (string): the string to be hashed

    Returns:
        int: the hash value
    """
    res = 0

    #We assigned values 3,2,1,0 to A,C,G,T respectively
    for i in range(len(s)):
        if s[i]=="A":
            res+=3*(4**(len(s)-i-1))
        if s[i]=="C":
            res+=2*(4**(len(s)-i-1))
        if s[i]=="G":
            res+=1*(4**(len(s)-i-1))
        if s[i]=="T":
            res+=0*(4**(len(s)-i-1))
            
    return res%(2**13-1)
  
def kmer_search(S,L):
    """Uses Rabin-Karp in the multi search version to find the
        frequencies of all kmers and their first appearence

    Args:
        S (string): the DNA sequence to be studied
        L (list): all the kmers to be looked up

    Returns:
        pos (int): the left-most position found for the most 
                    frequent(s) kmers, None if nothing was found
        max (int): the maximum frequence
     """
    
    #the lenght of the kmers in list L
    k = len(L[0])
    
    #An hashtable containing in each bucket a list of pairs.
    #Each pair is composed by a kmer and its index in list L  
    hash_table = [[] for _ in range(2**13-1)]
    
    #The first position of each kmer in list L
    first_pos = [-1]*(len(L))
    #The frequency of each kmer in list L
    freq = [0]*(len(L))
    
    #Updating the hashtable
    for i in range(len(L)):
        hash_table[init_hash(L[i])].append((L[i],i))
        
    #The hash of the beginning of the dna sequence, it will be updated in next steps
    h1 = init_hash(S[0:k])
    
    #Rabin-Karp in multi search version
    for i in range(len(S)-k+1):
        for string,index in hash_table[h1]:
            if S[i:i+k]==string:
                freq[index]+=1
                if first_pos[index]==-1: 
                    first_pos[index]=i
            
            #Updating the hash function
        if i+k < len(S):
            #the char to be eliminated from the old kmer
            old = 0
            if S[i]=="A":
                old = 3
            if S[i]=="C":
                old = 2
            if S[i]=="G":
                old = 1
            if S[i]=="T":
                old = 0
            #the char to be added to obtein the next kmer
            new = 0
            if S[i+k]=="A":
                new = 3
            if S[i+k]=="C":
                new = 2
            if S[i+k]=="G":
                new = 1
            if S[i+k]=="T":
                new = 0
                
            h1 = rol_hash(old,new,h1,4**(k-1))
        
    max = 0
    pos = -1
    for i in range(len(freq)):
        if freq[i] > max:
            max = freq[i]
            pos = first_pos[i]
    if pos == -1:
        pos = None
    return pos,max
        
    


#h,km = kmer_hist(S,8)

#print(f"Biggest frequency: {len(h)-1}\nMost frequent kmers: {km}")

with open(r"C:\Users\39345\Desktop\scientificcomputing\CleanedHuman.txt") as f:
    S = f.read()
start = time()
h,k = kmer_hist(S,8)
length = time()-start
print(f"It took {length} seconds")
"""
k = 3, 8, 15, 25, 35
"""

"""
fig,ax = plt.subplots()
ax.bar(range(len(h)),h)
ax.set_xlabel("Frequencies")
ax.set_ylabel("Number of 8-mers")
ax.set_title("8-spectrum of human chromosome 21")
ax.set_xlim(-100,10000)
ax.set_ylim(-10,250)
plt.savefig("",dpi=350)
plt.show()
"""


