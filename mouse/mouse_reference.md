# Building a mouse HipSTR reference
Here, we provide a step-by-step outline of how we built the **mm10** HipSTR reference.
For more detail about each step, please reference to the corresponding file in the
**human** subdirectory of this repository.

This script begins by using [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) to identify repeats on each chromosome and then
extensively filters and processes the resulting repeats using standard unix utilities such as [grep](http://linuxcommand.org/man_pages/grep1.html) and [awk](http://linuxcommand.org/man_pages/awk1.html).
It requires an executable version of *Tandem Repeats Finder* called **trf409.legacylinux64**. Please download it from the TRF link above and ensure the other utilities are installed before
proceeding to the steps below

First, let's download this repository as it contains the scripts we need:
```
git clone https://github.com/HipSTR-Tool/HipSTR-references.git
cd HipSTR-references
```

Next, let's download the FASTA files for the mm10 reference genome:
```
mkdir mm10
cd mm10
mkdir raw_fasta trf_results fixed_trf_results
cd raw_fasta
for chrom in $(seq 1 19) X Y
do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr$chrom.fa.gz
    gunzip chr$chrom.fa.gz &
done
cd ../../
```

We now use Tandem Repeats Finder to identify repeats on each chromosome:
```
chmod 755 trf409.legacylinux64
chmod 755 scripts/run_TRF.sh
for chrom in $(seq 1 19) X Y
do
    echo mm10/raw_fasta/chr$chrom.fa mm10/trf_results 5
done | xargs -L 1 -P 30 ./scripts/run_TRF.sh
```

We then filter out repeats with a period longer than 6 and fix a few issues with incorrect TRF entries:
```                                                                                                                                                                                         
for chrom in $(seq 1 19) X Y
do
    echo scripts/fix_trf_output.py mm10/trf_results/chr$chrom.fa mm10/fixed_trf_results/chr$chrom.fa
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
python scripts/trf_parser.py $files > filtered_repeats.mm10.bed
bedtools sort -i filtered_repeats.mm10.bed > filtered_repeats.mm10.sorted.bed
```

Next, we merge overlapping STRs into single entries and filter repeats that fail merging:
```
python scripts/analyze_overlaps.py filtered_repeats.mm10.sorted.bed pass.mm10 fail.mm10
```

We then remove any entries within 10bp of a failed merge region:
```
bedtools window -w 10 -a pass.mm10 -b fail.mm10 -v > pass.mm10.r2
```

To minimize the effects of nearby STRs on genotyping errors, we extract entries that aren't within 10bp
of another entry or are within 10bp of one or more entries that all share the same period
```
bedtools merge -i pass.mm10.r2 -c 4,6 -o collapse -d 10 | grep -v "," > pass.mm10.r3
bedtools merge -i pass.mm10.r2 -c 4,4,4,6 -o collapse,count_distinct,distinct,collapse -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6, $7}' | sed "s/,/\//g" >> pass.mm10.r3
```

Lastly, we construct the final reference for mm10 and delete the temporary files we created along the way
```
cat pass.mm10.r3 | bedtools sort | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "MOUSE_STR_"NR, $5}' > mm10.hipstr_reference.bed
rm -r mm10
rm fail.mm10 filtered_repeats.mm10.bed filtered_repeats.mm10.sorted.bed
rm pass.mm10 pass.mm10.r2 pass.mm10.r3
```
