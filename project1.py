"""
Claudiu Acsinte 233207
"""


"""
from time import time
import sys
"""

def m_shuffling(m,k,l,first,freq):
   
    """
    This function performs the m-shuffling
    
    Parameters:
    m (int): the number of the suffling
    k (int): the dimension of the input list "l" (actually half of it)
    l (list): the list to be shuffled
    first (list): vector such that first[h] is the first time h+1 has been on top
    freq (list): vector such that freq[h] is the number of times h+1 has been on top
    
    Returns:
    l (list): the list after the m-shuffling
    first (list): vector first after the m-shuffling
    freq (list): vector freq after the m-shuffling
    
    """
    
    """ here we store the first 2m elements of list "l" """
    temp = [0]*(2*m) 
    for i in range(2*m):
        temp[i] = l[i]
        
    """
    We noticed that element i of list "l" will be send
    in position 2*i+1 if i<m or in position (i-m)*2 if m<=i<2*m
    """    
    for i in range(m):
        l[2*i+1] = temp[i]
    for i in range(m,2*m):
        l[(i-m)*2] = temp[i]
        
    """
    If the element now on top has never been on top before we 
    update the list "first" using the current shuffling.
    Notice how we need to subtract 1 from the elements in list "l" in 
    order to access list "first" since the elements in "l" are the numbers from 
    1 to 2*k
    """
    if first[l[0]-1] == 0:
        first[l[0]-1] = m
    """
    We also update the list of frequencies
    """
    freq[l[0]-1] = freq[l[0]-1]+1
    
    return 
    
def complete_shuffle(k):
    
    """The actual function of our project

    Parameters:
    k (int): half of the dimension of our list
    
    Returns:
    L (list): the three first elements of our complete list after the complete shuffling
    s (int): the first time the top element has been on top
    f (int): number of times the top element has been on top
    
    """
    l = list(range(1,2*k+1))
    first = [0]*(2*k)
    freq = [0]*(2*k)
    for m in range(1,k+1):
        m_shuffling(m,k,l,first,freq)
    L = [l[0],l[1],l[2]]
    s = first[l[0]-1]
    f = freq[l[0]-1]

    return L,s,f


"""
Test


def main():
    k = 50000
    start = time()
    L,s,f = complete_shuffle(k)
    speed = time()-start
    print(f" For k = {k} we obteined L,s,f = {L} {s} {f}\n It took {speed} seconds")


if __name__ == "__main__":
    main()
"""

