# Building a mouse HipSTR reference
Here, we provide a step-by-step outline of how we built the **mm10** HipSTR reference.
For more detail about each step, please reference to the corresponding file in the
**human** subdirectory of this repository.

First, let's download the mm10 files
```
cd mm10
mkdir raw_fastq trf_results fixed_trf_results
cd raw_fastq
for chrom in $(seq 1 19) X Y
do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr$chrom.fa.gz
    gunzip chr$chrom.fa.gz &
done
cd ../../
```

Next, we use Tandem Repeats Finder to identify repeats
```
for chrom in $(seq 1 19) X Y
do
    echo mm10/raw_fastq/chr$chrom.fa mm10/trf_results 5
done | xargs -L 1 -P 30 ./run_TRF.sh
```

We then filter out repeats with a period longer than 6 and fix a few issues with incorrect TRF entries:
```                                                                                                                                                                                         
for chrom in $(seq 1 19) X Y
do
    echo fix_trf_output.py mm10/trf_results/chr$chrom.fa mm10/fixed_trf_results/chr$chrom.fa
done | xargs -L 1 -P 40 python
```

We reformat the TRF entries and filter to only include repeats with a sufficiently high score:
```
files=""
for chrom in $(seq 1 19) X Y
do
    files="$files,mm10/fixed_trf_results/chr$chrom.fa"
done
files=`echo $files | sed "s/,//"`
python trf_parser.py $files > filtered_repeats.mm10.bed
bedtools sort -i filtered_repeats.mm10.bed > filtered_repeats.mm10.sorted.bed
```

Next, we merge overlapping STRs into single entries and filter repeats that fail merging:
```
python analyze_overlaps.py filtered_repeats.mm10.sorted.bed pass.mm10 fail.mm10
```

We then remove any entries within 10bp of a failed merge region:
```
bedtools window -w 10 -a pass.mm10 -b fail.mm10 -v > pass.mm10.r2
```

To minimize the effects of nearby STRs on genotyping errors, we extract entries that aren't within 10bp
of another entry or are within 10bp of one or more entries that all share the same period
```
bedtools merge -i pass.mm10.r2 -c 4 -o collapse -d 10 | grep -v "," > pass.mm10.r3
bedtools merge -i pass.mm10.r2 -c 4,4,4 -o collapse,count_distinct,distinct -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6}' >> pass.mm10.r3
```

Lastly, we construct the final reference for mm10
```
cat pass.mm10.r3 | bedtools sort | awk -v OFS="\t" '{print $0, ($3-$2+1)/$4, "MOUSE_STR_"NR}' > mm10.hipstr_reference.bed
```
