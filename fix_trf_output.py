import sys

# Indices for various fields used in pattern correction
pattern_index   = 13
period_index    = 2
nrepeat_index   = 3
cons_size_index = 4

# Indices for fields used in overlap resolution 
start_index = 0
stop_index  = 1
score_index = 7

# Number of tokens in a correct input line
num_tokens = 17

# Run the Z-algorithm on the provided sequence. Returns an array of length len(NMER), where arr[i] corresponds
# to the length of the sequence starting at position i that is a prefix of NMER
def zalgorithm(nmer):
    res = len(nmer)*[0]
    l   = -1
    r   = -1

    for i in range(1, len(nmer), 1):
        if i > r:
            j = i
            while j < len(nmer) and nmer[j] == nmer[j-i]:
                j = j + 1

            res[i] = j-i
            l      = i
            r      = j-1
        else:
            i_prime = i-l
            beta    = r-i+1

            if res[i_prime] < beta:
                res[i] = res[i_prime]
            elif res[i_prime] == beta:
                j = r + 1
                while j < len(nmer) and nmer[j] == nmer[j-i]:
                    j = j + 1

                res[i] = j-i
                l      = i
                r      = j-1                
            else:
                res[i] = beta
    res[0] = len(nmer)
    return res

# Returns the repeating subunit of the provided sequence if it consists of exactly 2 or more copies. Otherwise, returns false
def has_subseq(nmer):
    prefix_lengths = zalgorithm(nmer)

    for k in range(1, len(nmer), 1):
        if len(nmer) % k == 0:
            match = True

            for segment in range(len(nmer)/k):
                coord = segment*k
               
                if prefix_lengths[coord] < k:
                    match = False
                    break

            if match:
                return nmer[0:k]
    return False

# Create a new file containing the same information as contained in INPUT_FILE,  except that the period, pattern and number of repeats for entries 
# with incorrect patterns are corrected. Also filters out entries with a corrected period greater than MAX_PERIOD
def correct_pattern_errors(input_file, output_file, max_period=6):
    data    = open(input_file, "r")
    output  = open(output_file, "w")

    # Read the chromosome heading and write it out to the filtered file
    output.write(data.readline())
 
    fail_count    = 0
    skip_count    = 0
    correct_count = 0

    for line in data:
        tokens = line.strip().split()
        
        if len(tokens) != num_tokens:
            exit("ERROR: Malformed input file")
            skip_count = skip_count + 1
            continue

        pattern = tokens[pattern_index]
        subseq  = has_subseq(pattern)
        if subseq:
            new_patt = subseq
            tokens[nrepeat_index]   = str(float(tokens[nrepeat_index])*len(pattern)/len(new_patt))
            tokens[pattern_index]   = new_patt
            tokens[period_index]    = str(len(new_patt))
            tokens[cons_size_index] = str(len(new_patt))
            new_line = ' '.join(tokens) + '\n'
            
            if len(new_patt) <= max_period:
                correct_count += 1
                output.write(new_line)
            else:
                fail_count += 1
        else:
            if int(tokens[cons_size_index]) <= max_period:
                output.write(line)
            else:
                fail_count += 1

    data.close()
    output.close()
    print("Modified %d records whose patterns were incorrect"%(correct_count))
    print("Removed  %d records whose patterns had a PERIOD > %d"%(fail_count, max_period))

def main():
    if len(sys.argv) != 3:
        exit("ERROR: This program accepts exactly two arguments: the input and output file names. Exiting...")
    
    print("FIXING FILE " + sys.argv[1])
    correct_pattern_errors(sys.argv[1], sys.argv[2], max_period=6)


if __name__ == "__main__":
    main()
