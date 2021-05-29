import re
import pandas as pd

# for reading files
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

# for formatting outputs:
def index_transform(matched_indices):
    return " ".join([str(i+1) for i in matched_indices]) # python index starts from 0

def convert_tuples_to_strings(lst_of_tuples):
    for tu in lst_of_tuples:
        print(str(tu[0])+' '+str(tu[1]))

def convert_dict_to_strings(dictionary):
    for key in dictionary:
        matched_indices = dictionary[key]
        if len(matched_indices) > 0:
            print(key)
            print(index_transform(matched_indices))

# scoring matrices:

with open('BLOSUM62.txt', 'r') as f:
        raw_BLOSUM62 = f.read().split()
raw_BLOSUM62 = raw_BLOSUM62[20:]
BLOSUM62 = {}

for i in range(20):
    BLOSUM62[raw_BLOSUM62[21*i]] = [int(x) for x in raw_BLOSUM62[21*i+1:21*i+21]]
BLOSUM62 = pd.DataFrame(BLOSUM62, index = BLOSUM62.keys())

with open('PAM250.txt', 'r') as f:
        raw_PAM250 = f.read().split()
raw_PAM250 = raw_PAM250[20:]
PAM250 = {}

for i in range(20):
    PAM250[raw_PAM250[21*i]] = [int(x) for x in raw_PAM250[21*i+1:21*i+21]]
PAM250 = pd.DataFrame(PAM250, index = PAM250.keys())

# for alignment sanity check
def check_score(s, t, f1, f2, max_score, multiple):
    '''
    s and t are reconstructed strings for global alignment
    f1 is a match score function
    f2 is a gap penalty function
    max_score is the maximum score for global alignment
    if MSA: multiple = true
    '''
    assert len(s) == len(t), "s and t has different lengths"
    i = 0
    gap_s = 0
    gap_t = 0
    score = 0
    while i < len(s):
        if s[i] == '-':
            assert s[i] != t[i], "two gaps cannot be at the same position"
            gap_s += 1
            score += f2(gap_t)
            gap_t = 0
        elif t[i] == '-':
            assert s[i] != t[i], "two gaps cannot be at the same position"
            gap_t += 1
            score += f2(gap_s)
            gap_s = 0
        else:
            score += f1(s[i], t[i])
            score += f2(gap_s)
            gap_s = 0
            score += f2(gap_t)
            gap_t = 0
        print(score)
        i+=1
        
    score += f2(gap_s)
    score += f2(gap_t)
    
    if not multiple:
        assert score == max_score, "s and t does not give rise to maximum alignment score"
    else:
        return score