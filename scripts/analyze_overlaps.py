import sys

def main():
    data       = open(sys.argv[1], "r")
    merge_pass = open(sys.argv[2], "w")
    merge_fail = open(sys.argv[3], "w")
    pass_count = 0
    fail_count = 0
    regions    = []
    cur_chrom  = ""
    min_start  = 0
    max_stop   = -1        

    for line in data:
        tokens              = line.strip().split("\t")
        chrom               = tokens[0]
        start, stop, period = map(int, tokens[1:4])
        motif               = tokens[4]
        num_repeats         = float(tokens[5])
        score               = float(tokens[6])
        sequence            = tokens[7]

        if chrom != cur_chrom or start > max_stop:
            # Analyze the previous set of overlapping regions
            region_size = max_stop - min_start + 1
            max_score   = 0
            max_index   = -1
            cov_frac    = 0.0
            for index,region in enumerate(regions):
                if region[5] >= max_score:
                    max_index = index
                    max_score = region[5]
                    cov_frac  = 1.0*(region[2]-region[1]+1)/region_size

            if max_index != -1:
                region = regions[max_index]
                if cov_frac > 0.85:
                    merge_pass.write("%s\t%s\t%d\t%d\t%.1f\t%s\n"%(region[0], region[1], region[2], region[3], region[4], region[6]))
                    pass_count += (len(regions) > 1)
                elif len(set(list(map(lambda x: x[3], regions)))) == 1:
                    all_motifs   = [reg[6] for reg in regions]
                    merge_pass.write("%s\t%s\t%d\t%d\t%.1f\t%s\n"%(region[0], min_start, max_stop, region[3], 1.0*(max_stop-min_start+1)/int(region[3]), "/".join(all_motifs)))
                else:
                    for region in regions:
                        merge_fail.write("%s\t%s\t%d\t%d\t%.1f\t%d\t%s\t%s\n"%(region[0], region[1], region[2], region[3], region[4], region[5], region[6], region[7]))
                    fail_count += (len(regions) > 1)

            # Reset the region info
            regions   = []
            min_start = start
            max_stop  = stop
            cur_chrom = chrom
        else:
            max_stop = max(max_stop, stop)

        regions.append((chrom, start, stop, period, num_repeats, score, motif, sequence))

    data.close()
    merge_pass.close()
    merge_fail.close()

if __name__ == "__main__":
    main()
