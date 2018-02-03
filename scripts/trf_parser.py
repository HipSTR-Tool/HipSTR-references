import sys

num_tokens     = 17
score_index    = 7
period_index   = 4  # Actually the size of the consensus sequence, but this is more accurate
motif_index    = 13
sequence_index = 14
start_index    = 0
stop_index     = 1
nrepeat_index  = 3

def min_perm(seq):
    min_perm = seq
    for i in xrange(len(seq)):
        other = seq[i:]+seq[0:i]
        if other < min_perm:
            min_perm = other
    return min_perm

def rev_complement(seq):
    rev     = seq[::-1]
    mapping = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    res     = ""
    for i in xrange(len(rev)):
        res = res+mapping[rev[i]]
    return res

def canonical_motif(motif):
    return min(min_perm(motif), min_perm(rev_complement(motif)))

def create_filtered_trf_bed_file(input_files):
    max_period        = 6
    period_vals       = [1,  2,  3,  4,  5,  6]
    period_thresholds = [20, 22, 28, 28, 32, 34]
    filt_count        = 0
    keep_count        = 0

    for i in range(len(input_files)):
        data  = open(input_files[i], "r")
        chrom = data.readline().strip()[1:]
        for line in data:
            tokens = line.strip().split()
            period = int(tokens[period_index])
            if period > max_period:
                filt_count += 1
                continue

            score = int(tokens[score_index])
            if score <= period_thresholds[period-1]:
                filt_count += 1
                continue

            new_tokens    = [chrom] + list(map(lambda x: tokens[x], [start_index, stop_index, period_index, motif_index, nrepeat_index, score_index, sequence_index]))
            new_tokens[4] = min_perm(new_tokens[4])
            print("\t".join(new_tokens))
            keep_count += 1

        data.close()
    sys.stderr.write("Kept %d out of %d records (%.2f%%) with sufficiently high scores\n"%(keep_count, keep_count+filt_count, 100.0*keep_count/(keep_count+filt_count)))
    return

def main():
    if len(sys.argv) != 2:
        exit("ERROR: Program requires exactly one arguments: a comma-separated list of file paths")
    files  = sys.argv[1].strip().split(",")
    create_filtered_trf_bed_file(files)

if __name__ == "__main__":
    main()
