# !/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import re
import scipy
from itertools import combinations

# http://rosalind.info/problems/dna/
def count_base(dna):
    dna = dna.upper()
    return [dna.count(b) for b in ['A', 'C', 'G', 'T']]
count_base('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC')

# http://rosalind.info/problems/rna/
def transcribe(dna):
    dna = dna.upper()
    return dna.replace('T', 'U')
transcribe('GATGGAACTTGACTACGTAAATT')

# http://rosalind.info/problems/revc/
def reverse_complement(dna):
    base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}  
    dna = dna.upper()
    return ''.join(map(lambda x: base_pair[x], dna))[::-1]
reverse_complement('AAAACCCGGT')

# http://rosalind.info/problems/fib/
def fibonacci(n, k):
    '''
    Given: Positive integers n≤40 and k≤5
    Return: The total number of rabbit pairs that will be present after n months, 
            if we begin with 1 pair and in each generation, 
            every pair of reproduction-age rabbits (>= 2 months old) produces 
            a litter of k rabbit pairs (instead of only 1 pair).
    '''
    total_number = [1]*n
    for i in range(2, n):
        total_number[i] = total_number[i-1] + k * total_number[i-2]
    return total_number[-1]
        

fibonacci(5,3)

# http://rosalind.info/problems/fibd/

def mortal_fibonacci(n, m): # some problems (integer overflow)
    '''
    Given: Positive integers n≤100 and m≤20
    Return: The total number of pairs of rabbits that will remain after the nth month 
            if all rabbits live for m months. (each rabbit reproduces only when >= 2 month old)
    '''
    total_number = np.zeros((m+1, n))
    # the entry at jth row and ith column represents the number of pair of rabbit 
    # that is j month old at the i+1 th month
    total_number[0, 0] = 1
    total_number[1, 1] = 1
    for i in range(2, n):
        for j in range(1, min(i, m)+1):
            #print(i, j)
            #print(total_number[j-1,i-1])
            total_number[j, i] = total_number[j-1, i-1]
            # update the age of rabbit in each month
        total_number[0, i] = sum(total_number[2:(m+1), i])
        # only the rabbits >= 2 months old can reproduce
    #print(total_number)
    return int(sum(total_number[:m, n-1]))


def mortal_fibonacci(n, m): # use python list, no integer overflow problem
    '''
    Given: Positive integers n≤100 and m≤20
    Return: The total number of pairs of rabbits that will remain after the nth month 
            if all rabbits live for m months. (each rabbit reproduces only when >= 2 month old)
    '''
    total_number = []
    for i in range(n):
        lst = []
        for j in range(m+1):
            lst.append(0)
        total_number.append(lst)
    # the entry at jth column and ith row represents the number of pair of rabbit 
    # that is j month old at the i+1 th month
    total_number[0][0] = 1
    total_number[1][1] = 1
    for i in range(2, n):
        for j in range(1, min(i, m)+1):
            total_number[i][j] = total_number[i-1][j-1]
            # update the age of rabbit in each month
        total_number[i][0] = sum(total_number[i][2:(m+1)])
        # print(total_number)
        # only the rabbits >= 2 months old can reproduce
    #print(total_number)
    return int(sum(total_number[n-1][:m]))

# http://rosalind.info/problems/gc/
def GC_content(dna):
    return np.mean([s == 'G' or s == 'C' for s in dna if s in ['A', 'T', 'G', 'C']])*100
GC_content('CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT')

def read_fasta(file):
    with open(file, 'r') as f:
        raw = f.read().split('>')
        raw.pop(0)
    dnas = {}
    for entry in raw:
        entry = entry.replace('\n', '')
        match = re.search(r'Rosalind_\d+', entry)
        dnas[match.group(0)] = entry[(match.span()[1]):]
    return dnas

def get_highest_GC(file):
    dnas = read_fasta(file)
    GCs = {x: GC_content(dnas[x]) for x in dnas}
    winner = max(GCs)
    print(winner+'\n'+str(round(GCs[winner], 7)))
get_highest_GC('rosalind_gc.txt')

# http://rosalind.info/problems/hamm/
def hamming_dist(dna1, dna2):
    assert len(dna1)==len(dna2)
    return sum([dna1[i] != dna2[i] for i in range(len(dna1))])

# http://rosalind.info/problems/perm/
def rearrange(n):
    '''
    > rearrange(3)
    [[3, 1, 2], [3, 2, 1], [1, 3, 2], [1, 2, 3], [2, 1, 3], [2, 3, 1]]
    '''
    permutations = {1: [[1]], 2: [[1,2], [2,1]]}
    
    def switch(item, i, j):
        '''
        switch i to j and j to i and keep all other entries as the same
        '''
        if item == i:
            return j
        elif item == j:
            return i
        else:
            return item
        
    for i in range(3, n+1):
        total_arrangement = []
        temp_permut_minus_1 = [[i] + x for x in permutations[i-1]] 
        total_arrangement.extend(temp_permut_minus_1.copy())
        # e.g. for i = 3, this is [3,1,2] and [3,2,1]
        for j in range(1, i):
            switched_permut = [list(map(lambda item: switch(item, i, j), lst)) for lst in temp_permut_minus_1]
            total_arrangement.extend(switched_permut)
            # switch between i and j e.g. for j = 1, this is [1,3,2] and [1,2,3]
        permutations[i] = total_arrangement
    return permutations[n]

def map_lsts_to_strings(lst):
    '''
    > map_lsts_to_strings([[3, 1, 2], [3, 2, 1], [1, 3, 2], [1, 2, 3], [2, 1, 3], [2, 3, 1]])
    6
    3 1 2
    3 2 1
    1 3 2
    1 2 3
    2 1 3
    2 3 1
    '''
    print(len(lst))
    for item in lst:
        print(' '.join([str(x) for x in item]))
map_lsts_to_strings(rearrange(3))

# http://rosalind.info/problems/sign/
def signed_permutation(n):
    '''
    >2
    -1 -2
    -1 2
    1 -2
    1 2
    -2 -1
    -2 1
    2 -1
    2 1
    '''
    def map_to_sign(item, i):
        '''
        convert i to -i in a lst
        '''
        if item == i:
            return -i
        else:
            return item
        
    unsigned_permutation = rearrange(n)
    for i in range(1, n+1):
        # double the current list by adding the possibility of converting each i to -i for all i = 1, 2, 3, ..., n
        new_permutation = [list(map(lambda item: map_to_sign(item, i), lst)) for lst in unsigned_permutation]
        unsigned_permutation.extend(new_permutation)
    return unsigned_permutation

def map_lst_to_strings_in_file(lst, filename):
    with open(filename, 'w+') as f:
        f.writelines([str(len(lst))+'\n'])
        f.writelines([' '.join([str(x) for x in item])+'\n' for item in lst])
map_lst_to_strings_in_file(signed_permutation(3), 'test.txt')

# http://rosalind.info/problems/lexf/
def lexicograph(lst, n):
    '''
    Given: A collection of at most 10 symbols defining an ordered alphabet, and a positive integer n(n≤10).
    Return: All strings of length n that can be formed from the alphabet, 
            ordered lexicographically (use the standard order of symbols in the English alphabet).
    '''
    
    lst = [s.upper() for s in lst]
    lst.sort()
    possible_strings = {1:lst.copy()}
    for i in range(2, n+1):
        temp_strings = []
        for s in lst:
            temp_strings.extend([s + string for string in possible_strings[i-1]])
        possible_strings[i] = temp_strings
    for string in possible_strings[n]:
        print(string)
    return possible_strings[n]
lexicograph(['A','G', 'T', 'C'], 2)

# http://rosalind.info/problems/kmer/
def kmer_composition(dna, k):
    kmer_lst = lexicograph(['A', 'T', 'G', 'C'], k)
    freq = [0]*4**k
    for i in range(4**k):
        freq[i] = len(find_motif(dna, kmer_lst[i], 'list'))
    print(" ".join([str(x) for x in freq]))
    return freq

# http://rosalind.info/problems/iprb/
def mendel_prob(k, m, n):
    '''
    Given: Three positive integers k, m, and n , representing a population containing k+m+n organisms:
           k individuals are homozygous dominant for a factor, m are heterozygous, and n are homozygous recessive.
    Return: The probability that two randomly selected mating organisms will produce an individual possessing a dominant allele 
            (and thus displaying the dominant phenotype). Assume that any two organisms can mate.
    '''
    # two dominant: choose(k, 2) p = 1
    # two heterozygous: choose(m, 2) p = 3/4
    # two recessive: choose(n, 2) p = 0
    # dominant + heterozygous: k*m, p = 1
    # dominant + recessive: k*n, p = 1
    # heterozygous + recessive: m*n, p = 1/2
    # total choices: choose(k+m+n, 2)
    return (k*(k-1)/2 * 1 + m*(m-1)/2 * 3/4 + n*(n-1)/2 * 0 + k*m*1 + k*n*1 + m*n/2)/((k+m+n)*(k+m+n-1)/2)
mendel_prob(2,2,2)

# http://rosalind.info/problems/prot/
from Bio.Seq import Seq
def translate(mrna):
    return str(Seq(mrna).translate())[:-1] #the last one is a stop codon

# http://rosalind.info/problems/mrna/
import Bio.Data.CodonTable
codon_table = Bio.Data.CodonTable.standard_rna_table.forward_table
reverse_codon_table = {}
for codon in codon_table:
    aa = codon_table[codon]
    if aa not in reverse_codon_table:
        reverse_codon_table[aa] = [codon]
    else:
        reverse_codon_table[aa].append(codon)

def infer_mrna(protein):
    '''
    Given: A protein string of length at most 1000 aa.
    Return: The total number of different RNA strings from which the protein could have been translated, 
            modulo 1,000,000. (Don't neglect the importance of the stop codon in protein translation.)
    '''
    total_number = 3 # 3 options for stop codon
    for aa in protein:
        total_number = total_number * len(reverse_codon_table[aa]) % 1000000
    return total_number

infer_mrna('MA')

# http://rosalind.info/problems/subs/
def find_motif(s, t, output = 'string'):
    '''
    Given: Two DNA strings s and t (each of length at most 1 kbp).
    Return: All locations of t as a substring of s.
    '''
    possible_indices = [str(i+1) for i in range(len(s)-len(t)+1) if s[i:i+len(t)] == t]
    if output == 'list':
        return possible_indices
    return ' '.join(possible_indices)
find_motif('GATATATGCATATACTT',
'ATAT')

# http://rosalind.info/problems/cons/
def consensus_string(file):
    dnas = read_fasta(file)
    consensus_dict = {'A':[], 'C':[], 'G':[], 'T':[]}
    for i in range(len(dnas[list(dnas.keys())[0]])):
        for base in consensus_dict:
            consensus_dict[base].append(sum([dna[i] == base for dna in dnas.values()]))
    consensus_matrix = pd.DataFrame(consensus_dict).transpose()
    print(''.join(list(consensus_matrix.idxmax())))
    for base in consensus_dict:
        print(base + ': ' + ' '.join([str(count) for count in consensus_dict[base]]))
consensus_string('rosalind_cons.txt')


# http://rosalind.info/problems/kmp/
# kmf algorithm reference: https://www.inf.hs-flensburg.de/lang/algorithmen/pattern/kmpen.htm
def failure_array(dna, filename = False):
    max_lengths = [0]
    for k in range(1, len(dna)):
        max_lengths.append(0)
        # the sup for prefix = suffix length at position k is reached if we can extend the max length at k-1 by 1
        # i.e. if dna[max_lengths[k-1]] = dna[k]
        if dna[max_lengths[k-1]] == dna[k]:
            max_lengths[k] = max_lengths[k-1] + 1
        else:
            for j in range(max_lengths[k-1], 0, -1):
                if dna[:j] == dna[(k-j+1):(k+1)]:
                    max_lengths[k] = j
                    break
    if filename:
        with open(filename, 'w+') as f:
            f.write(' '.join([str(l) for l in max_lengths]))
    return max_lengths
failure_array('CAGCATGGTATCACAGCAGAG')  # failure array = lps

#reference: https://towardsdatascience.com/pattern-search-with-the-knuth-morris-pratt-kmp-algorithm-8562407dba5b


def find_lps(pattern):

    # name lps indicates longest proper prefix which is also suffix.
    prefix = 0
    lps_array = [0]*len(pattern)

    for i in range(1, len(pattern)):

        while pattern[prefix] != pattern[i] and prefix > 0:
            prefix = lps_array[prefix-1]
            # use proof by contradiction to show that there is no possible match between
            # lps_array[prefix-1] and prefix

        if pattern[prefix] == pattern[i]:
            prefix += 1
            lps_array[i] = prefix

    return lps_array

def kmp(pattern, text, unique):

    # Knuth-Morris-Pratt algorithm

    prefix = 0 # use prefix to keep track of matches
    lps_array = find_lps(pattern)
    matched_indices = []

    for i, ch in enumerate(text):

        while pattern[prefix] != ch and prefix > 0:
            prefix = lps_array[prefix-1]
        
        if pattern[prefix] == ch:

            if prefix == len(pattern) - 1: # full match
                matched_indices.append(i - prefix)

                if unique: # if we know the match is uqniue, we can return the match as soon as we find one
                    return matched_indices

                prefix = lps_array[prefix] # we turn the suffix into prefix and find the next match

            else:
                prefix += 1

    return matched_indices

def index_transform(matched_indices):
    return " ".join([str(i+1) for i in matched_indices]) # python index starts from 0

# http://rosalind.info/problems/splc/
def exon_protein(file):
    '''
    first dna in file is the entire sequence
    the rest dnas are introns
    
    Return: a protein made by concatanating the codons
    '''
    dnas = read_fasta(file, 'list')
    text = dnas[0]
    patterns = dnas[1:]
    for pattern in patterns:
        start = kmp(pattern, text, True)[0]
        end = start + len(pattern) - 1
        text = text[0:start] + text[end+1:]
    raw_protein = translate(text)
    return raw_protein.replace('*', '') # remove stop codons

exon_protein('rosalind_splc.txt')

# http://rosalind.info/problems/sseq/
def find_spliced_motif(s, t):
    '''
    Given: Two DNA strings s and t (each of length at most 1 kbp) in FASTA format.
    Return: One collection of indices of s in which the symbols of t appear as a subsequence of s. 
            If multiple solutions exist, you may return any one.
    '''
    s_index = 0
    t_index = 0
    subsequence = []
    while t_index < len(t):
        if s[s_index] == t[t_index]:
            subsequence.append(s_index)
            s_index += 1
            t_index += 1
        else:
            s_index += 1
    return index_transform(subsequence)

find_spliced_motif('ACGTACGTGACG', 'GTA')

# http://rosalind.info/problems/revp/
def check_palindrome(dna):
    if len(dna) % 2 != 0:
        return False
    forward_pointer = 0
    backward_pointer = len(dna)-1
    base_pair = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    while forward_pointer < backward_pointer:
        if dna[forward_pointer] == base_pair[dna[backward_pointer]]:
            forward_pointer += 1
            backward_pointer -= 1
        else:
            return False
    return True

def check_palindrome_brute_force(dna):
    return reverse_complement(dna) == dna
            
    
def find_palindrome(dna):

    '''
    Return a list of tuples:
    each tuple: (i, l)
    where i is the position of start base in the dna (index starts from 1)
    and l is the length of reverse palindromes
    if there are mutlple l for oen particular i, return the largest l
    '''
    
    palindrome_list = []
    
    for i in range(len(dna)):
        
        match_indices = []
        
        candidate_palindromes = [] # candidate palindromes of all lengths
        
        j = i + 2
        
        while j < len(dna):
            
            candidate_j = 0
            
            pattern = reverse_complement(dna[i:j]) 
            '''
            for example, if the ith index is 'A', and pattern = 'AG',
            the reverse complement of 'AG' is 'CT'
            we want to see if 'CT' is in dna,
            if not, there is no palindrome that starts from 'A'
            if yes, we find the largest length of palindrome based on all the 'CT' matches
            and increase j until there is no palindrome that starts from 'A'
             '''
            #print(dna[i:j], pattern)
            matched_indices = [x for x in kmp(pattern, dna) if x > i] 
            #print(matched_indices)
            
            if matched_indices:
                k = len(matched_indices) - 1
                while k >= 0:
                    candidate = dna[i:(matched_indices[k]+(j-i))]
                    if check_palindrome(candidate):
                        
                        candidate_j = len(candidate)
                        #print('yes', len(candidate))
                        break
                    else:
                        k -= 1
            
            else:
                break
                
            if candidate_j > 0:
                candidate_palindromes.append(candidate_j)
            
            #print(candidate_j)
            
            j += 1
        #print(candidate_palindromes)
            
        if candidate_palindromes:
            palindrome_list.append((i+1, max(candidate_palindromes)))
            # i+1 to account for python index start from 0
        
    return palindrome_list

def find_palindrome_brute_force(dna):
    '''
    Return a list of tuples:
    each tuple: (i, l)
    where i is the position of start base in the dna (index starts from 1)
    and l is the length of reverse palindromes
    there could be multiple l for each i
    '''
    palindrome_candidates = []
    for i in range(len(dna)-3):
        for j in range(i+4, i+13):
            if j > len(dna):
                break
            candidate = dna[i:j]
            if check_palindrome_brute_force(candidate):
                palindrome_candidates.append((i+1, len(candidate)))
    return palindrome_candidates

            
find_palindrome('TCAATGCATGCGGGTCTATATGCAT')  

def convert_tuples_to_strings(lst_of_tuples):
    for tu in lst_of_tuples:
        print(str(tu[0])+' '+str(tu[1]))
        
convert_tuples_to_strings(find_palindrome('TCAATGCATGCGGGTCTATATGCAT'))

class SuffixTree: # still some problem

    def __init__(self, lst_of_strings):
        self.strings = lst_of_strings
        self.terminators = {lst_of_strings[i]:str(i+1) for i in range(len(lst_of_strings))}
        self.root = Node()
        for string in self.strings:
            self.createTree(string)
    
    def createTree(self, string): # constructor for one string, need to modify to add later strings
        ter = self.terminators[string]
        #self.root.edge = {string+ter: Node()}
        #string = string[1:]
        while len(string) > 0:
            lst = self.getLastNode(string, self.root)
            if lst[2] == 'end at node':
                return
            elif lst[2] == 'end within edge':
                node = lst[0]
                edge = node.getEdge()
                label = lst[1]
                new_label_1 = lst[3][0]
                new_label_2 = lst[3][1]
                edge[new_label_1+ter] = Node(edge = {new_label_2:edge[label]})
                edge.pop(label, None)
            else:
                node = lst[0]
                edge = node.getEdge()
                new_label = lst[3]
                edge[new_label+ter] = Node()
            string = string[1:]


    
    def getLastNode(self, string, node):
        edge = node.getEdge()
        for label in node.getEdge():
            l = min([len(label), len(string)])
            i = 0
            while i < l:
                if label[i] != string[i]:
                    break
                else:
                    i+=1
            '''
            some problems with this implementation:
            consider a suffix tree with xab and we want to add xac
            we have reached xa, which has length smaller than 3
            we want to go to the next node of xa 
            if xa is a leaf node, we add c to the leaf node to make it a non-leaf node
            if xa is not a leaf node, we continue to the next node,
            if xa does not have any edges containing c, we add c to the node xa
            if xa has an edge c, we do nothing
            if xa has an edge that starts with c but includes more, we split this edgeinto two edges
            connected by a node


            '''
            if i == 1:
                continue
            elif i < len(string):
                return (self.getLastNode(string[len(label):], edge[label]))
            elif i == len(string) and i == len(label):
                return [edge[label], None, 'end at node', '']
            else: #(i == len(string) < len(label))
                return [node, label, 'end within edge', [string, label[i:]]]
        return [node, None, 'no matches', string]
    
class Node:

    def __init__(self, edge = {}):
        # edge is a dictionary with key-value pair:
        # key: string, value: next node
        self.edge = edge


    def getEdge(self):
        return self.edge

# http://rosalind.info/problems/lcsm/

from suffix_trees import STree

def dna_lcs(file):
    # longest common string
    dnas = read_fasta(file, 'list')
    st = STree.STree(dnas)
    return st.lcs()

# http://rosalind.info/problems/mprt/

import requests
from bs4 import BeautifulSoup

def glycolysation_motif(file):
    with open(file, 'r+') as f:
        protein_ids = f.read().split('\n')
    returned = {}
    for protein in protein_ids:
        url = "https://www.uniprot.org/uniprot/" + protein + ".fasta"
        page = requests.get(url)
        soup = BeautifulSoup(page.content, 'html.parser')
        seq = ''.join(str(soup).split('\n')[1:-1])
        matches = re.finditer(r'(?=(N[ARNDCQEGHILKMFSTWYV][ST]))', seq) # r'(?=(pattern))' --> find overlapping matches
        # N-glycosylation motif: Asn-X-Serine/Threonine where X is any amino acid other than Proline
        returned[protein] = []
        for m in matches:        
            returned[protein].append(m.start())
    return returned

def convert_dict_to_strings(dictionary):
    for key in dictionary:
        matched_indices = dictionary[key]
        if len(matched_indices) > 0:
            print(key)
            print(index_transform(matched_indices))
            
convert_dict_to_strings(glycolysation_motif('rosalind_mprt.txt'))   

# http://rosalind.info/problems/orf/
def find_ORFs(dna):
    
    mrna = dna.replace('T', 'U')
    complement_mrna = reverse_complement(dna).replace('T', 'U')
    start_strand1 = list(re.finditer(r'(AUG)', mrna))
    start_strand2 = list(re.finditer(r'(AUG)', complement_mrna))
    # both strands (strand2 is the reverse complement of strand1) can be coding strand --> 6 ORs
    
    proteins = set()
    
    for starts, seq in [(start_strand1, mrna), (start_strand2, complement_mrna)]:
        
        for m in starts:
            
            i = m.start()
            while seq[i:i+3] not in ['UAG', 'UAA', 'UGA'] and i < len(seq):
                i += 3
            if i < len(seq):
                p = translate(seq[m.start():i+3])
                proteins.add(p)
    
    for p in proteins:
        
        print(p)

find_ORFs('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
        
# http://rosalind.info/problems/grph/
def overlap_graph(file):
    '''
    For a collection of strings and a positive integer k, the overlap graph for the strings is a directed graph Ok
    in which each string is represented by a node, and string s is connected to string t with a directed edge 
    when there is a length k suffix of s that matches a length k prefix of t, as long as s≠t; we demand s≠t
    to prevent directed loops in the overlap graph (although directed cycles may be present).
    
    Given: A collection of DNA strings in FASTA format having total length at most 10 kbp.
    Return: The adjacency list corresponding to O3. You may return edges in any order.
    '''
    dnas = read_fasta(file)
    overlap_dnas = set()
    for dna_id1 in dnas:
        for dna_id2 in dnas:
            if dna_id1 != dna_id2:
                if dnas[dna_id1][-3:] == dnas[dna_id2][:3]:
                    overlap_dnas.add(dna_id1 + ' ' + dna_id2)
                elif dnas[dna_id2][-3:] == dnas[dna_id1][:3]:
                    overlap_dnas.add(dna_id2 + ' ' + dna_id1)
    for item in overlap_dnas:
        print(item)

# http://rosalind.info/problems/pmch/
from scipy.special import factorial

def perfect_matching(rna):
    '''
    return total number of perfect matching basepair edges
    '''
    assert rna.count('A') == rna.count('U') and rna.count('G') == rna.count('C')
    return factorial(rna.count('A'), exact = True) * factorial(rna.count('G'), exact = True)
    # for i = 1, 2, ... rna.count('A'):
    # for each Ui, we have n-i+1 ways to choose Aj to pair with Ui
    

perfect_matching('AGCUAGUCAU')

# http://rosalind.info/problems/cat/
total_numbers = {'':1,'A':0, 'G':0, 'U':0, 'C':0, 'AA':0, 'GG':0, 'UU':0, 'CC':0,
                    'AG':0, 'AC':0, 'CA':0, 'CU':0, 'GA':0, 'GU':0, 'UC':0, 'UG':0,
                    'AU':1, 'UA':1, 'CG':1, 'GC':1} 
# we need to add '' to account for the edge that connects the first base to the last base
def noncrossing(rna, memo):
    '''
    return total number of noncrossing perfect matching basepair edges
    '''
    if rna not in memo:
        if rna.count('A') != rna.count('U') or rna.count('G') != rna.count('C'):
            memo[rna] = 0
        else:
            memo[rna] = sum([noncrossing(rna[0]+rna[i], memo)*(noncrossing(rna[1:i], memo))*(noncrossing(rna[i+1:], memo)) for i in range(1, len(rna), 2)])
            # we fix the first base, and try to connect it with every other base to create a separating line --> rna[0]+rna[i]
    return (memo[rna]) % 1000000

noncrossing('AUAU', total_numbers)

# http://rosalind.info/problems/motz/
def noncrossing_all(rna, memo):
    '''
    return total number of noncrossing matching basepair (could be non-perfect) edges
    '''
    if rna not in memo:
        if len(rna) <= 1:
            memo[rna] = 1
        else:
            memo[rna] = noncrossing_all(rna[1:], memo) # base 1 is not involved in any base pair
            base_pair = {'A':'U', 'U':'A', 'C':'G', 'G':'C'}  
            for i in range(1, len(rna)):
                if base_pair[rna[0]] == rna[i]:
                    memo[rna] += noncrossing_all(rna[1:i], memo) * noncrossing_all(rna[i+1:], memo)
            # we fix the first base, and try to connect it with every other base to create a separating line --> rna[0]+rna[i]
    return (memo[rna]) % 1000000
noncrossing_all('AUAU', {})

# http://rosalind.info/problems/rnas/
def noncrossing_with_wobble(rna, memo):
    '''
    return total number of noncrossing matching basepair (could be non-perfect) edges (including U-G basepair)
    two bases can pair only if they are at least 4 bases apart
    '''
    if rna not in memo:
        if len(rna) <= 3:
            memo[rna] = 1
        else:
            memo[rna] = noncrossing_with_wobble(rna[1:], memo) # base 1 is not involved in any base pair
            base_pair = {'A':['U'], 'U':['A', 'G'], 'C':['G'], 'G':['U','C']}
            for i in range(4, len(rna)):
                if rna[i] in base_pair[rna[0]]:
                    memo[rna] += noncrossing_with_wobble(rna[1:i], memo) * noncrossing_with_wobble(rna[i+1:], memo)
            # we fix the first base, and try to connect it with every other base to create a separating line --> rna[0]+rna[i]
    return (memo[rna])
noncrossing_with_wobble('GCGGUUCGAGUGCCGUCAUCCGAGGAAUCGUAUGCAUAGUAGUAGGAUUGGACCUGCAAGUAGCAGCGUUUCAACCGUUGUUACCGGAGUAUACUGGAGGGGGUCGCCACUCGGACGUGCAUCCAGUAUUUCAGUUGUUAGCAGUGAUGAUCGUUUGGAUUUAUGGUCUUUUUGUU', {})


# https://towardsdatascience.com/solving-tsp-using-dynamic-programming-2c77da86610d

def overlap_dist(dna1, dna2):
    candidate_start_index = [i for i,val in enumerate(dna1) if val==dna2[0]]
    for i in candidate_start_index:
        j1 = i+1
        j2 = 1
        while j1 < len(dna1):
            if dna1[j1] == dna2[j2]:
                j1+=1
                j2+=1
            else:
                break
        if j1 == len(dna1):
            return j2
    return 1
def distance_array(file, output = 'list'):
    dnas = read_fasta(file, 'list')
    if output == 'numpy':
        dist = np.zeros((len(dnas), len(dnas)))
        for i in range(len(dnas)):
            for j in range(len(dnas)):
                dist[i,j] = overlap_dist(dnas[i], dnas[j])
    else:
        dist = [[0]*len(dnas) for i in range(len(dnas))]
        for i in range(len(dnas)):
            for j in range(len(dnas)):
                dist[i][j] = overlap_dist(dnas[i], dnas[j])
    return (dist, dnas)

def shortest_common_superstring(file): # needs an approximate solution, complexity: https://math.stackexchange.com/questions/190313/time-complexity-of-the-travelling-salesman-problem
    dist, dnas = distance_array(file)
    memo = {(tuple([i]), i): tuple([None, 0]) for i in range(len(dnas))} 
    # key: prev_visited_indices, prev_last_index; value: next_to_last_index, total dist thus far
    queue = [(tuple([i]), i) for i in range(len(dnas))]
    all_string_indices = set(range(len(dnas)))

    while queue:
        prev_visited_indices, prev_last_index = queue.pop(0)
        _, prev_dist = memo[(prev_visited_indices, prev_last_index)]
        to_visit = all_string_indices.difference(set(prev_visited_indices))
        
        for next_last_index in to_visit:
            new_visited_indices = tuple(sorted(list(prev_visited_indices)+[next_last_index]))
            new_dist = prev_dist + dist[prev_last_index][next_last_index]
            
            if (new_visited_indices, next_last_index) not in memo:
                memo[(new_visited_indices, next_last_index)] = (prev_last_index, new_dist)        
                queue.append((new_visited_indices, next_last_index))
            elif new_dist > memo[(new_visited_indices, next_last_index)][1]: 
                # we want maximum overlap and thus greatest distance between each string
                memo[(new_visited_indices, next_last_index)] = (prev_last_index, new_dist)
                # we do not need to append this one to queue, because as in BFS, 
                # we always trace all paths to a smaller subset before heading toward a larger subset
                # so (new_visited_indices, next_last_index) should still be in the queue
    
    shortest_superstring = retrace_optimal_path(memo, dnas)
    return shortest_superstring

def retrace_optimal_path(memo, dnas):
    superstring = ""
    points_to_retrace = tuple(range(len(dnas)))
    full_path_memo = dict((k,v) for k,v in memo.items() if k[0] == points_to_retrace)
    # only includes keys with full paths
    optimal_choice = max(full_path_memo.keys(), key = lambda x: full_path_memo[x][1])
    last_index = optimal_choice[1]
    superstring += dnas[last_index]
    next_to_last_index, curr_dist = memo[optimal_choice]
    points_to_retrace = tuple(sorted(set(points_to_retrace).difference({last_index})))

    while next_to_last_index is not None:
        last_index = next_to_last_index
        next_to_last_index, prev_dist = memo[(points_to_retrace, last_index)]
        overlap = curr_dist - prev_dist
        superstring = dnas[last_index][:-overlap] + superstring
        curr_dist = prev_dist
        points_to_retrace = tuple(sorted(set(points_to_retrace).difference({last_index})))
    
    return superstring
        

# http://rosalind.info/problems/long/ :The dataset is guaranteed to satisfy the following condition: there exists a unique way to reconstruct the entire chromosome from these reads by gluing together pairs of reads that overlap by more than half their length.

# create a graph 
# topologically sort the graph
# hamiltonian path

def create_graph(file):
    dist, dnas = distance_array(file)
    adjacency_list = {i: [] for i in range(len(dnas))} 
    # key: vertex; value: each vertex where there is a directed edge from key to that vertex
    reverse_adjacency_list = {i: [] for i in range(len(dnas))}
    # key: vertex; value: each vertex where there is a directed edge from that vertex to key

    for i in range(len(dnas)):
        for j in range(len(dnas)):
            if i != j and dist[i][j] >= max(len(dnas[i]), len(dnas[j])) / 2:
                adjacency_list[i].append(j)
                reverse_adjacency_list[j].append(i)
    return (adjacency_list, reverse_adjacency_list, dist, dnas)

def topological_sort(adjacency_list, reverse_adjacency_list):
    
    visited = []
    order = []
    indegree_0_vertices = [v for v in reverse_adjacency_list.keys() if not reverse_adjacency_list[v]]
    
    def dfs(v):
        visited.append(v)
        for v2 in adjacency_list[v]:
            if v2 not in visited:
                dfs(v2)
        order.append(v)
    
    for v in indegree_0_vertices:
        dfs(v)
    
    order = order[::-1]
    return order

def scs(file):
    
    adjacency_list, reverse_adjacency_list, dist, dnas = create_graph(file)
    order = topological_sort(adjacency_list, reverse_adjacency_list)
    prev_vertex = order[0]
    superstring = dnas[prev_vertex]
    
    i = 1
    while i < len(dnas):
        next_vertex = order[i]
        overlap = dist[prev_vertex][next_vertex]
        superstring += dnas[next_vertex][overlap:]
        prev_vertex = next_vertex
        i += 1
    
    return superstring
        
        
shortest_common_superstring('rosalind_long-2.txt') == scs('rosalind_long-2.txt')

# http://rosalind.info/problems/corr/
def correct_reads(file):
    dnas = read_fasta(file, 'list')
    dna_counts = {}
    reads = {'correct': [], 'incorrect': []}
    changes = []
    for dna in set(dnas):
        rc_dna = reverse_complement(dna)
        if dna not in dna_counts and rc_dna not in dna_counts:
            dna_counts[dna] = dnas.count(dna) + dnas.count(rc_dna)
            if dna_counts[dna] > 1:
                reads['correct'].append(dna)
            else:
                reads['incorrect'].append(dna)
                
    for error_dna in reads['incorrect']:
        for correct_dna in reads['correct']:
            if hamming_dist(error_dna, correct_dna) == 1:
                changes.append(error_dna + '->' + correct_dna)
    for change in changes:
        print(change)

correct_reads('rosalind_corr.txt')


# http://rosalind.info/problems/tree/
def find_set_index(lst_of_sets, elem):
    '''
    Given a list of sets, return the index of the set in which the element is in 
    '''
    for i in range(len(lst_of_sets)):
        if elem in lst_of_sets[i]:
            return i
    print(lst_of_sets, elem)

def compress(lst_of_sets, index1, index2):
    '''
    Given a list of sets, combine the set at index1 and the set at index2 into one set
    '''
    lst_of_sets[index1] = lst_of_sets[index1].union(lst_of_sets[index2])
    lst_of_sets.pop(index2)
    return lst_of_sets

def connected_components(file):
    '''
    return the number of connected components in a graph
    '''
    with open(file, 'r+') as f:
        texts = f.readlines()
    texts = [re.sub('\n', '', text) for text in texts]
    #print(texts)
    n = int(texts[0])
    connected = [{i} for i in range(1, n+1)]
    for text in texts[1:]:
        t1, t2 = re.findall(r'(\d+)\s(\d+)', text)[0]
        index1 = find_set_index(connected, int(t1))
        index2 = find_set_index(connected, int(t2))
        connected = compress(connected, index1, index2)
    return len(connected)

def edges_to_add(file):
    '''
    The minimum number of edges to be added so that a graph becomes a tree (a connected graph) 
    is the number of connected components - 1
    '''
    return connected_components(file) - 1

# an even simpler way is to realize that a mimimum of n-1 edges can turn n vertices into a connected graph
# and once the graph is connected, adding additional edges will create cycles

# http://rosalind.info/problems/inod/
def internal_nodes(n):
    '''
    Return: the number of internal nodes (k) in an unrooted binary tree with n leaves
    Each leaf has degree 1
    All internal nodes in an unrooted binary tree have degree 3
    2*edge number = 3k + n
    k+n total vertices --> k+n-1 total edges
    2*(k+n-1) = 3k+n --> k = n - 2
    '''
    return n-2

# http://rosalind.info/problems/lcsq/
def longest_common_subsequence(s, t):
    longest_string = [['']*len(t) for i in range(len(s))]
    longest_string_length = [[0]*len(t) for i in range(len(s))]

    for i in range(len(s)):
        if i == 0:
            for j in range(len(t)):
                if s[i] == t[j]:
                    longest_string_length[i][j] = 1
                    longest_string[i][j] = s[0]
                else:
                    if j > 0:
                        longest_string_length[i][j] = longest_string_length[i][j-1]
                        longest_string[i][j] = longest_string[i][j-1]
        else:
            for j in range(len(t)):
                if j == 0:
                    if s[i] == t[j]:
                        longest_string_length[i][j] = 1
                        longest_string[i][j] = t[0]
                else:
                    candidate_lengths = {(i-1, j-1): 0, (i-1, j): 0, (i, j-1): 0}
                    candidate_strings = {(i-1, j-1): '', (i-1, j): '', (i, j-1): ''}
                    if s[i] == t[j]:
                        candidate_lengths[(i-1, j-1)] = 1
                        candidate_strings[(i-1, j-1)] = s[i]
                    for index_s, index_t in candidate_lengths.keys():
                        candidate_lengths[(index_s, index_t)] += longest_string_length[index_s][index_t]
                        candidate_strings[(index_s, index_t)] = longest_string[index_s][index_t]+candidate_strings[(index_s, index_t)]
                    max_index = max(candidate_lengths.keys(), key = lambda x: candidate_lengths[x])
                    longest_string_length[i][j] = candidate_lengths[max_index]
                    longest_string[i][j] = candidate_strings[max_index]
    return longest_string[-1][-1]

longest_common_subsequence('AACCTTGG', 'ACACTGTGA')

# http://rosalind.info/problems/edit/
def edit_distance(s, t):
    min_edit = [[0]*len(t) for i in range(len(s))]
    for i in range(len(s)):
        if i == 0:
            for j in range(len(t)):
                min_edit[i][j] = 1-int(s[i] in t[:j+1]) + j # if s[i] in t[:j+1], return 0+j, else, return 1+j
                
        else:
            for j in range(len(t)):
                if j == 0:
                    min_edit[i][j] = 1-int(t[j] in s[:i+1]) + i
                else:
                    if s[i] == t[j]:
                        min_edit[i][j] = min_edit[i-1][j-1]
                    else:
                        min_edit[i][j] = min([min_edit[i-1][j-1]+1,
                                              min_edit[i-1][j]+1,
                                              min_edit[i][j-1]+1
                                             ])
    print(min_edit[-1][-1])
    return min_edit

edit_distance('PLEASANTLY', 'MEANLY')

# http://rosalind.info/problems/edta/
def edit_reconstruct(s, t):
    '''
    it is easier to reconstruct according to the min_edit table than 
    reconstructing during creating the min_edit table
    '''
    min_edit = edit_distance(s, t)
    s_reconstruct = ''
    t_reconstruct= ''
    i = len(s) - 1
    j = len(t) - 1
    while i > 0 and j > 0:
        if min_edit[i][j] == min_edit[i-1][j-1] + 1-int(s[i] == t[j]):
            s_reconstruct = s[i] + s_reconstruct
            t_reconstruct = t[j] + t_reconstruct
            i -= 1
            j -= 1
        elif min_edit[i][j] == min_edit[i][j-1] + 1:
            s_reconstruct = '-' + s_reconstruct
            t_reconstruct = t[j] + t_reconstruct
            j -= 1
        else: # min_edit[i][j] == min_edit[i-1][j] + 1
            s_reconstruct = s[i] + s_reconstruct
            t_reconstruct = '-' + t_reconstruct
            i -= 1
    if i == 0 and j == 0:
        s_reconstruct = s[0] + s_reconstruct
        t_reconstruct = t[0] + t_reconstruct
    elif i == 0:
        if s[0] in t[:j+1]:
            match_index = t[:j+1].index(s[0])
            s_reconstruct = '-'*match_index + s[0] + '-'*(j-match_index) + s_reconstruct
            # insert s[0] at exactly the place where it matches t
        else:
            s_reconstruct = '-'*j + s[0] + s_reconstruct
        t_reconstruct = t[:j+1] + t_reconstruct
    else:
        if t[0] in s[:i+1]:
            match_index = s[:i+1].index(t[0])
            t_reconstruct = '-'*match_index + t[0] + '-'*(i-match_index) + t_reconstruct
        else:
            t_reconstruct = '-'*i + t[0] + t_reconstruct
        s_reconstruct = s[:i+1] + s_reconstruct
    
    print(s_reconstruct)
    print(t_reconstruct)
    with open('result.txt', 'w+') as f:
        f.writelines([str(min_edit[-1][-1])+'\n', s_reconstruct+'\n', t_reconstruct+'\n'])
    assert (not any([s_reconstruct[i] == '-' and t_reconstruct == '-' for i in range(len(s_reconstruct))]))
    assert (hamming_dist(s_reconstruct, t_reconstruct) == min_edit[-1][-1])
        
edit_reconstruct('PRETTY', 'PRTTEIN')

# http://rosalind.info/problems/glob/
with open('BLOSUM62.txt', 'r') as f:
        raw_BLOSUM62 = f.read().split()
raw_BLOSUM62 = raw_BLOSUM62[20:]
BLOSUM62 = {}

for i in range(20):
    BLOSUM62[raw_BLOSUM62[21*i]] = [int(x) for x in raw_BLOSUM62[21*i+1:21*i+21]]
BLOSUM62 = pd.DataFrame(BLOSUM62, index = BLOSUM62.keys())

def global_alignment_score_linear(s, t, scoring_matrix = BLOSUM62, gap_penalty = 5):
    '''
    linear gap penalty of -5 per gap
    '''
    max_score = [[0]*(len(t)+1) for i in range((len(s)+1))]
    # max_score[i][j] stores the max score when aligning s[:i] against s[:j]
    # hence max_score[0][0] is always 0
    
    for i in range(len(s)+1):
        max_score[i][0] = - i * gap_penalty
        
    for j in range(len(t)+1):
        max_score[0][j] = - j * gap_penalty
                
    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            max_score[i][j] = max([max_score[i-1][j-1] + scoring_matrix.loc[s[i-1], t[j-1]],
                                              max_score[i-1][j] - gap_penalty,
                                              max_score[i][j-1] - gap_penalty
                                             ])
    print(max_score[-1][-1])
    return max_score

global_alignment_score_linear('PLEASANTLY', 'MEANLY')

proteins = read_fasta('rosalind_glob.txt', 'list')
x = global_alignment_score_linear(proteins[0], proteins[1])
        
# http://rosalind.info/problems/gaff/
def global_alignment_affine(s, t, scoring_matrix = BLOSUM62, start_penalty = 11, extend_penalty = 1):
    '''
    reference: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
    returns not only score but also alignment results
    '''
    Sn = [[0]*(len(t)+1) for i in range(len(s)+1)] # no gap at both ends
    Ss = [[0]*(len(t)+1) for i in range(len(s)+1)] # gap at end of s
    St = [[0]*(len(t)+1) for i in range(len(s)+1)] # gap at end of t
    # a gap cannot be at both the end of s and the end of t
    
    Pn = [[0]*(len(t)+1) for i in range(len(s)+1)]
    Ps = [[0]*(len(t)+1) for i in range(len(s)+1)]
    Pt = [[0]*(len(t)+1) for i in range(len(s)+1)]
    # each P matrix traces the matrix of the previous indices that leads to the max score at the current index
    # 0 represents the matrix Sn
    # 1 represents the matrix Ss
    # 2 represents the matrix St
    
    lst_match_matrices = [Sn, Ss, St]
    lst_prev_matrices = [Pn, Ps, Pt]
    
    for i in range(1, len(s)+1):
        Sn[i][0] = -float('inf') # we try to match s[:i] with "", can only match xxxxxxx with -------
        Ss[i][0] = -float('inf')
        St[i][0] = -start_penalty -(i-1)*extend_penalty
        
    for j in range(1, len(t)+1):
        Sn[0][j] = -float('inf')
        Ss[0][j] = -start_penalty -(i-1)*extend_penalty
        St[0][j] = -float('inf')
    
    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            curr_score = scoring_matrix.loc[s[i-1], t[j-1]]
            alignment_scores_no_gap = {0: Sn[i-1][j-1] + curr_score,
                                       1: Ss[i-1][j-1] + curr_score,
                                       2: St[i-1][j-1] + curr_score}
            alignment_scores_gap_s = {0: Sn[i][j-1] - start_penalty,
                                      2: St[i][j-1] - start_penalty,
                                      1: Ss[i][j-1] - extend_penalty}
            alignment_scores_gap_t = {0: Sn[i-1][j] - start_penalty,
                                      1: Ss[i-1][j] - start_penalty,
                                      2: St[i-1][j] - extend_penalty}
                    
            Pn[i][j] = max([0,1,2], key = lambda x: alignment_scores_no_gap[x])
            Ps[i][j] = max([0,1,2], key = lambda x: alignment_scores_gap_s[x])
            Pt[i][j] = max([0,1,2], key = lambda x: alignment_scores_gap_t[x])
                    
            Sn[i][j] = alignment_scores_no_gap[Pn[i][j]]
            Ss[i][j] = alignment_scores_gap_s[Ps[i][j]]
            St[i][j] = alignment_scores_gap_t[Pt[i][j]]
                    
    max_score_s_t = max([Sn[-1][-1], Ss[-1][-1], St[-1][-1]])
    print(max_score_s_t)
    
    # add trace step
    
    s_reconstruct = ''
    t_reconstruct= ''
    i = len(s) 
    j = len(t) 
    
    curr_matrix_index = max([0,1,2], key = lambda x: lst_match_matrices[x][i][j])
    prev_matrix_index = lst_prev_matrices[curr_matrix_index][i][j]
                  
    while i > 0 and j > 0:
        if curr_matrix_index == 0:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1
        elif curr_matrix_index == 1:
            s_reconstruct = '-' + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1
        else:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = '-' + t_reconstruct
            i-=1
        curr_matrix_index = prev_matrix_index
        prev_matrix_index = lst_prev_matrices[curr_matrix_index][i][j]


    if i == 0 and j == 0:
        x = "nothing" # do nothing
    elif i == 0:
        s_reconstruct = '-'*j + s_reconstruct
        t_reconstruct = t[:j] + t_reconstruct
    else: # j == 0
        s_reconstruct = s[:i] + s_reconstruct
        t_reconstruct = '-'*i + t_reconstruct
            
    print(s_reconstruct)
    print(t_reconstruct)

# http://rosalind.info/problems/gcon/
def global_alignment_score_constant(s, t, scoring_matrix = BLOSUM62, gap_penalty = 5):
    '''
    constant gap penalty of -5 per contiguous stretch of gaps
    constant gap penalty = affine gap penalty with start penalty = constant and extend penalty = 0
    reference: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf
    '''
    global_alignment_affine(s, t, scoring_matrix, start_penalty=gap_penalty, extend_penalty=0)
    

global_alignment_score_constant('PLEASANTLY', 'MEANLY')

proteins = read_fasta('rosalind_gcon.txt', 'list')
global_alignment_score_constant(proteins[0], proteins[1])

# http://rosalind.info/problems/loca/
def local_alignment_linear(s, t, scoring_matrix = PAM250, gap_penalty = 5):
    '''
    linear gap penalty
    Smith–Waterman algorithm
    reference: https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/local.pdf
    '''        
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    
    
    # we want to track the best score and the location where that best score is taken 
    best = 0
    opt_indices = (0, 0)
    
    # the first row and first column have entry 0
    
    for i in range(1, len(s)+1):
        for j in range(1, len(t)+1):
            max_score[i][j] = max([max_score[i-1][j-1] + scoring_matrix.loc[s[i-1], t[j-1]],
                                              max_score[i-1][j] - gap_penalty,
                                              max_score[i][j-1] - gap_penalty,
                                             0])
            if max_score[i][j] > best:
                best = max_score[i][j]
                opt_indices = (i, j)
                
    print(best)
    # trace back:
    i, j = opt_indices
    s_reconstruct = ''
    t_reconstruct = ''
    
    while max_score[i][j] > 0:
        
        if max_score[i][j] == max_score[i-1][j-1] + scoring_matrix.loc[s[i-1], t[j-1]]:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i -= 1
            j -= 1
        
        elif max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = '-' + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j -= 1
        
        else:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = '-' + t_reconstruct
            i -= 1
            
    print(s_reconstruct.replace('-', ''))
    print(t_reconstruct.replace('-', ''))
    print(s_reconstruct)
    print(t_reconstruct)


local_alignment_linear('MEANLYPRTEINSTRING', 'PLEASANTLYEINSTEIN')

# http://rosalind.info/problems/laff/

def local_alignment_affine(s, t, scoring_matrix = BLOSUM62, start_penalty = 11, extend_penalty = 1):
    '''
    Affine gap penalty
    Smith–Waterman algorithm
    reference: http://cseweb.ucsd.edu/classes/fa09/cse182/slides/L3.pptx.pdf
    This function is slow for proteins above 1000 aa lengths
    Converted this function to laff.java
    '''    
    def local_alignment_affine(s, t, scoring_matrix = BLOSUM62, start_penalty = 11, extend_penalty = 1):
    '''
    linear gap penalty
    Smith–Waterman algorithm
    reference: http://cseweb.ucsd.edu/classes/fa09/cse182/slides/L3.pptx.pdf
    '''    
    Sn = np.zeros((len(t)+1)) # no gap
    Ss = np.zeros((len(t)+1)) # gap at end of s
    St = np.zeros((len(t)+1)) # gap at end of t

    backtrack = 3 * np.ones((len(s)+1, len(t)+1))

    # we no longer need to keep prev_matrix for all 3 matrices
    # since gap_at_end_of_s and gap_at_end_of_t are guaranteed no beter than no_gap_at_end for local alignment
    
    best, opt_indices = 0, (0,0) 
        
                
    for i in range(1, len(s)+1):
        new_Sn = np.zeros((len(t)+1))
        new_Ss = np.zeros((len(t)+1))
        new_St = np.zeros((len(t)+1))
        
        # first entries are all 0
                
        for j in range(1, len(t)):

            curr_score = scoring_matrix.loc[s[i-1], t[j-1]]
            
            # max_score[0][i][j-1] is guaranteed to be better than max_score[2][i][j-1]
            # max_score[0][i-1][j] is guaranteed to be better than max_score[1][i-1][j]
            
            new_Ss[j] = max([new_Sn[j-1] - start_penalty,
                             new_Ss[j-1] - extend_penalty])
            new_St[j] = max([Sn[j] - start_penalty,
                             St[j] - extend_penalty])
            candidates = [Sn[j-1] + curr_score,
                          new_Ss[j],
                          new_St[j],
                          0]
            
            new_Sn[j] = max(candidates)
            backtrack[i,j] = candidates.index(new_Sn[j])
            
            if new_Sn[j] > best:
                best, opt_indices = new_Sn[j], (i, j)
                
        Sn = new_Sn
        Ss = new_Ss
        St = new_St
                
    print(int(best))
    
    # trace back:
    
    max_i, max_j = opt_indices
    i = max_i
    j = max_j
    
    while backtrack[i,j] != 3:
        
        if backtrack[i,j] == 0:
            i -= 1
            j -= 1
        
        elif backtrack[i,j] == 1:
            j -= 1
        
        else:
            i -= 1
    
    print(s[i:max_i])
    print(t[j:max_j])
    
local_alignment_affine('PLEASANTLY', 'MEANLY')

# http://rosalind.info/problems/smgb/
def semiglobal_alignment_linear_dna(s, t, match_score = 1, substitution_penalty = 1, gap_penalty = 1):
    '''
    1. same initialization as in local alignment to account for starting gap
    2. same iteration as in global alignment
    3. best = max(last row and last column of score matrix) to account for ending gap
    reference: http://www.cs.cmu.edu/~durand/03-711/2015/Lectures/PW_sequence_alignment_2015.pdf
    4. write a smgb.java equivalent to deal with large (1000-10000 bp) dnas
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    
    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty
            max_score[i][j] = max([max_score[i-1][j-1] + one_match,
                                   max_score[i-1][j] - gap_penalty,
                                   max_score[i][j-1] - gap_penalty])
    
    best = max([max([max_score[i][-1] for i in range(len(s)+1)]),
                max(max_score[-1])])
    print(best)
    
    # traceback:
    
    s_reconstruct = ""
    t_reconstruct = ""
    
    if best == max([max_score[i][-1] for i in range(len(s)+1)]):
        # gap at end of t
        max_i = [max_score[i][-1] for i in range(len(s)+1)].index(best)
        s_reconstruct = s[max_i:]
        t_reconstruct = '-'*(len(s)-max_i)
            
        i = max_i
        j = len(t)
        
    else: # best == max(max_score[-1])
        max_j = max_score[-1].index(best)
        t_reconstruct = t[max_j:]
        s_reconstruct = '-'*(len(t)-max_j)
            
        j = max_j
        i = len(s)
            
    while i*j > 0:
        one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty

        if max_score[i][j] == max_score[i-1][j-1] + one_match:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1

        elif max_score[i][j] == max_score[i-1][j] - gap_penalty:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = "-" + t_reconstruct
            i-=1

        else: # max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = "-" + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1
                
    if i == 0 and j == 0: # no gap at start
        x = "nothing" # do nothing
        
    elif i == 0: # gap at start of s
        s_reconstruct = '-'*j + s_reconstruct
        t_reconstruct = t[:j] + t_reconstruct
        
    else: # j == 0, gap at start of t
        t_reconstruct = '-'*i + t_reconstruct
        s_reconstruct = s[:i] + s_reconstruct   
    
    print(s_reconstruct)
    print(t_reconstruct)
    
semiglobal_alignment_linear_dna('CAGCACTTGGATTCTCGG', 'CAGCGTGG')

# http://rosalind.info/problems/oap/
def overlap_alignment_linear_dna(s, t, match_score = 1, substitution_penalty = 2, gap_penalty = 2):
    '''
    1. same initialization as in local alignment to account for trivial suffix of s ('') match with ''
    and '' match with trivial prefix of t ('')
    2. same iteration as in global alignment, so that makes sure always prefix of t match with s
    3. best = max(last row of score matrix) to make sure end at i == -1, so that always suffix of s match with t
    4. traceback until j = 0
    5. write a oap.java equivalent to deal with large (1000-10000 bp) dnas
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    
    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty
            max_score[i][j] = max([max_score[i-1][j-1] + one_match,
                                   max_score[i-1][j] - gap_penalty,
                                   max_score[i][j-1] - gap_penalty])
    
    best = max(max_score[-1])
    
    print(best)
    
    # traceback:
    
    s_reconstruct = ""
    t_reconstruct = ""
    

    max_j = max_score[-1].index(best)    
    j = max_j
    i = len(s)
            
    while j > 0:
        one_match = match_score if s[i-1] == t[j-1] else - substitution_penalty

        if max_score[i][j] == max_score[i-1][j-1] + one_match:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            i-=1
            j-=1

        elif max_score[i][j] == max_score[i-1][j] - gap_penalty:
            s_reconstruct = s[i-1] + s_reconstruct
            t_reconstruct = "-" + t_reconstruct
            i-=1

        else: # max_score[i][j] == max_score[i][j-1] - gap_penalty:
            s_reconstruct = "-" + s_reconstruct
            t_reconstruct = t[j-1] + t_reconstruct
            j-=1
                
    print(s_reconstruct)
    print(t_reconstruct)

overlap_alignment_linear_dna('CTAAGGGATTCCGGTAATTAGACAG', 'ATAGACCATATGTCAGTGACTGTGTAA')

# two slightly different approaches for: http://rosalind.info/problems/ctea/

def count_optimal_alignments(s, t):
    '''
    recursion + memoization approach
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    for i in range(0, len(s)+1):
        max_score[i][0] = -i
    for j in range(0, len(t)+1):
        max_score[0][j] = -j
    index_with_tied_matches = {}
        
    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = 0 if s[i-1] == t[j-1] else -1
            '''
            remember edit distance scoring scheme is that matching symbol does not contribute to score,
            but substitutions and gap contributes -1
            '''
            candidates = [max_score[i-1][j-1] + one_match,
                          max_score[i-1][j] - 1,
                          max_score[i][j-1] - 1]
            max_score[i][j] = max(candidates)
            
            is_max = [cand == max_score[i][j] for cand in candidates]
            candidate_indices = [(i-1, j-1), (i-1, j), (i, j-1)]
            index_with_tied_matches[(i,j)] = [candidate_indices[k] for k in range(3) if is_max[k]]
            
    count_matrix = [[0]*(len(t)+1) for i in range(len(s)+1)]

    modulus = 2**27 - 1
        
    def count(indices):
        i, j = indices
        if i == 0 or j == 0:
            count_matrix[i][j] = 1
            return 1
        elif count_matrix[i][j] != 0:
            return count_matrix[i][j]
        else:
            prev_indices = index_with_tied_matches[(i,j)]
            count_matrix[i][j] = sum([count(prev) for prev in prev_indices]) % modulus
            return count_matrix[i][j]
        
    return count((len(s),len(t))) 
    
count_optimal_alignments('PLEASANTLY', 'MEANLY')

def count_optimal_alignments(s, t):
    '''
    dynamic programming approach
    '''
    max_score = [[0]*(len(t)+1) for i in range(len(s)+1)]
    count_matrix = [[0]*(len(t)+1) for i in range(len(s)+1)]
    
    modulus = 2**27 - 1
    
    for i in range(0, len(s)+1):
        max_score[i][0] = -i
        count_matrix[i][0] = 1
    for j in range(0, len(t)+1):
        max_score[0][j] = -j
        count_matrix[0][j] = 1
        
    for i in range(1,len(s)+1):
        for j in range(1, len(t)+1):
            one_match = 0 if s[i-1] == t[j-1] else -1
            '''
            remember edit distance scoring scheme is that matching symbol does not contribute to score,
            but substitutions and gap contributes -1
            '''
            candidates = [max_score[i-1][j-1] + one_match,
                          max_score[i-1][j] - 1,
                          max_score[i][j-1] - 1]
            max_score[i][j] = max(candidates)
            
            is_max = [cand == max_score[i][j] for cand in candidates]
            candidate_indices = [(i-1, j-1), (i-1, j), (i, j-1)]
            
            for cand_indices in [candidate_indices[k] for k in range(3) if is_max[k]]:
                i0, j0 = cand_indices
                count_matrix[i][j] += count_matrix[i0][j0]
            count_matrix[i][j] = count_matrix[i][j] % modulus

    return count_matrix[len(s)][len(t)]

dnas = read_fasta('rosalind_ctea.txt', 'list')
count_optimal_alignments(dnas[0], dnas[1])

# http://rosalind.info/problems/mult/
def generalized_needleman_wunsch(dnas):
    '''
    Generalized Needleman-Wunsch algorithm for n dimension alignment
    Iteration rule: we consider all possible prev_index (see http://readiab.org/book/0.1.2/2/3)
    score_to_add for each prev_index = sum of all possible pairs of -int(s[i-1] != t[j-1]) - number of gaps * non_gaps
    '''
    dims = [len(dna)+1 for dna in dnas]
    n = len(dnas)
    indices = list(range(n))
    S = np.zeros(tuple(dims))
    num_entries = 1
    for dim in dims:
        num_entries *= dim
    traceback = np.array([" "*n for i in range(num_entries)]).reshape(tuple(dims)) 
    # when initlaizing strings, make the length of string equal to the length of string to be put during iteration
    
    # initialization, we actually only need to initialize S[0,0,...,0] = 0
    '''
    for k in range(n):
        loc = [0] * n
        for i in range(dims[k]):
            loc[k] = i
            S[tuple(loc)] = -(n-1)*i
            
    '''

    # iteration:
    
    # generate all combinations of indices
    all_possible_indices = np.meshgrid(*[list(range(0, dims[k])) for k in range(n)])
    for k in range(n):
        while len(all_possible_indices[k].shape) > 1:
            all_possible_indices[k] = np.hstack(all_possible_indices[k]) # flatten the indices
    all_possible_indices = [[all_possible_indices[k][i] for k in range(n)] for i in range(len(all_possible_indices[0]))]
    
    '''
    Note that in the 2-dimensional case, we only care about (i,0) and (0,j) during initialization
    But in higher dimensions, we have multiple borderline hyperplanes to consider:
    i.e. for 4-dimensional case: we have to consider hyperplanes like (i,j,0,0)
    '''

    for loc in all_possible_indices[1:]: # skip (0,0,0,0)
        candidates = {}
        for k in range(n, 0, -1):
            # we get nCk combinations of choosing index to minus 1
            possible_combinations = combinations(indices, k)
            
            for comb in possible_combinations:
                
                new_comb = [i for i in range(n) if (loc[i] > 0 and i in comb)]
                # i.e. if we have comb = (0,1,2) and loc = [0,0,1,1], we can only have prev index as [0,0,0,1]
                if len(new_comb) == 0:
                    continue # i.e. if we have comb = (0,1,2) and loc = [0,0,0,1], we can't minus the first 3 entries
                num = len(new_comb)
                #print(loc, comb, new_comb)
                cand = [loc[i]-1*int(i in new_comb) for i in range(n)]
                score_to_add = 0
                
                for i in new_comb:
                    for j in new_comb:
                        if i < j:
                            score_to_add -= int(dnas[i][cand[i]] != dnas[j][cand[j]])
                            # add all the scores between different pairs
                score_to_add -= (n-num)*num # gap penalty: (n-num) gaps match with num non-gaps
                
                #print(loc, cand, score_to_add)
                #print(S[tuple(cand)])
                candidates[tuple(cand)] = S[tuple(cand)] + score_to_add
                
        #print(loc, candidates)
        traceback[tuple(loc)] = "".join([str(s) for s in list(max(candidates.keys(), key = lambda x: candidates[x]))])
        # store the prev_index as string
        S[tuple(loc)] = candidates[max(candidates.keys(), key = lambda x: candidates[x])]
        
    #print(S)
    print(S[tuple(all_possible_indices[-1])])
    
    loc = all_possible_indices[-1]
    reconstructs = ["" for k in range(n)]
    
    while not any([index == 0 for index in loc]):
        prev_loc = [int(x) for x in traceback[tuple(loc)]] # unpack indices stored as string
        for k in range(n):
            if prev_loc[k] == loc[k] - 1:
                reconstructs[k] = dnas[k][loc[k] - 1] + reconstructs[k]
            else:
                reconstructs[k] = '-' + reconstructs[k]
        loc = prev_loc
    
    gap_length = 0
    for k in range(n):
        if loc[k] != 0:
            gap_length = max(gap_length, loc[k]) # get the max gap length to be put at the start
            
    for k in range(n):
        reconstructs[k] = '-' * (gap_length - loc[k]) + dnas[k][:loc[k]] + reconstructs[k]
    
    max_score = 0
    for i in range(n):
        for j in range(i+1, n):
            max_score -= sum([reconstructs[i][k] != reconstructs[j][k] for k in range(len(reconstructs[0]))])
    print(max_score)
    
    for reconstruct in reconstructs:
        print(reconstruct)

dnas = read_fasta('rosalind_mult.txt', 'list')
generalized_needleman_wunsch(dnas)