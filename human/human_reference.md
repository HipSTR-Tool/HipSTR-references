# Building an STR reference file for HipSTR
Here we describe how we built the HipSTR reference for the hg19 and hg38 versions of the human reference genome. We begin by scanning each
version of the reference genome for repeats using [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html). We then
extensively filter and process the resulting repeats using standard unix utilities such as [grep](http://linuxcommand.org/man_pages/grep1.html) and [awk](http://linuxcommand.org/man_pages/awk1.html), 
as well as three command line utilities: [datamash](https://www.gnu.org/software/datamash/), [bedtools](http://bedtools.readthedocs.io/en/latest/) and [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).

## Scanning the genome for STRs
First, let's download the relevant reference files for these two version of the human reference genome:
```
mkdir hg19 hg38                                                                                                                                     
# Download the hg19 files
cd hg19
mkdir raw_fastq trf_results fixed_trf_results
cd raw_fastq
for chrom in $(seq 1 22) X Y
do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr$chrom.fa.gz
    gunzip chr$chrom.fa.gz &
done
cd ../../

# Download the hg38 files                                                                                                                                                                                          
cd hg38
mkdir raw_fastq trf_results fixed_trf_results
cd raw_fastq
for chrom in $(seq 1 22) X Y
do
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr$chrom.fa.gz
    gunzip chr$chrom.fa.gz &
done
cd ../../
```

Next, we used [Tandem Repeats Finder (TRF)](https://tandem.bu.edu/trf/trf.html) to identify repeats in each of these chromosomes.
We first downloaded an executable version of this program called **trf409.legacylinux64**. We then scanned each sequence using the **run_TRF.sh** script included in this repository as follows:
```
for chrom in $(seq 1 22) X Y
do
    echo hg19/raw_fastq/chr$chrom.fa hg19/trf_results 5
    echo hg38/raw_fastq/chr$chrom.fa hg38/trf_results 5
done | xargs -L 1 -P 30 ./run_TRF.sh
```

As STRs are classically defined as having 1-6bp repeat motifs, we filtered the TRF output to only include these repeats. The script **fix_trf_output.py** performs this filtering and also corrects TRF entries in which the repeat size is incorrectly reported:
```
for chrom in $(seq 1 22) X Y
do
    echo fix_trf_output.py hg19/trf_results/chr$chrom.fa hg19/fixed_trf_results/chr$chrom.fa
    echo fix_trf_output.py hg38/trf_results/chr$chrom.fa hg38/fixed_trf_results/chr$chrom.fa
done | xargs -L 1 -P 40 python
```


## Filtering the STRs

```
files=""
for chrom in $(seq 1 22) X Y
do
    files="$files,hg19/fixed_trf_results/chr$chrom.fa"
done
files=`echo $files | sed "s/,//"`
python trf_parser.py $files > filtered_repeats.hg19.bed
bedtools sort -i filtered_repeats.hg19.bed > filtered_repeats.hg19.sorted.bed
```

```
files=""
for chrom in $(seq 1 22) X Y
do
    files="$files,hg38/fixed_trf_results/chr$chrom.fa"
done
files=`echo $files | sed "s/,//"`
python trf_parser.py $files > filtered_repeats.hg38.bed
bedtools sort -i filtered_repeats.hg38.bed > filtered_repeats.hg38.sorted.bed
```


## Merging STRs within and across assemblies


Merge overlapping STRs into single entries and filter those repeats that fail merging
```
python analyze_overlaps.py filtered_repeats.hg19.sorted.bed pass.hg19 fail.hg19
```

Remove any entries within 10bp of a failed merge region
```
bedtools window -w 10 -a pass.hg19 -b fail.hg19 -v > pass.hg19.r2
```

Extract entries that aren't within 10bp of another entry or are within 10bp of one or more entries that all share the same period
```
bedtools merge -i pass.hg19.r2 -c 4 -o collapse -d 10 | grep -v "," > pass.hg19.r3
bedtools merge -i pass.hg19.r2 -c 4,4,4 -o collapse,count_distinct,distinct -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6}' >> pass.hg19.r3
```

liftOver coordinates to hg38 and remove entries that failed to lift
```
./liftOver -minMatch=1.0 pass.hg19.r3  hg19ToHg38.over.chain.gz pass.hg19.r3.lifted_to_hg38 unlifted
grep -v ^# unlifted | cat - pass.hg19.r3 | sort -k 1,1V -k 2,2n | uniq -u > pass.hg19.r4
```

Remove if multiple entries lifted to the same hg38 region
```
./liftOver -minMatch=1.0 pass.hg19.r4 hg19ToHg38.over.chain.gz pass.hg19.r4.lifted_to_hg38 unlifted
paste pass.hg19.r4.lifted_to_hg38 pass.hg19.r4 | datamash -g 1,2,3,4 -s collapse 5,6,7,8 | grep -v "," | cut -f 5-8 > pass.hg19.r5
```

Remove entries that lift over to hg38 regions within 10bp of another region
```
./liftOver -minMatch=1.0 pass.hg19.r5 hg19ToHg38.over.chain.gz pass.hg19.r5.lifted_to_hg38 unlifted
paste pass.hg19.r5.lifted_to_hg38 pass.hg19.r5 | bedtools sort | bedtools merge -i - -c 4,5,6,7,8 -o collapse -d 10 | grep "," | cut -f 5- | awk -v OFS="\t" '{split($1,a,","); split($2,b,","); split($3, c, ",");
cat pass.hg19.r5 to_remove.bed | sort -k 1,1 -k 2,2n | uniq -u > pass.hg19.r6
./liftOver -minMatch=1.0 pass.hg19.r6 hg19ToHg38.over.chain.gz pass.hg19.r6.lifted_to_hg38 unlifted
```

Perform a similar set of analyses for GRCh38, but without the complications of liftOver filtering
```
python analyze_overlaps.py filtered_repeats.hg38.sorted.bed pass.hg38 fail.hg38
bedtools window -w 10 -a pass.hg38 -b fail.hg38 -v > pass.hg38.r2
bedtools merge -i pass.hg38.r2 -c 4 -o collapse -d 10 | grep -v "," > pass.hg38.r3
bedtools merge -i pass.hg38.r2 -c 4,4,4 -o collapse,count_distinct,distinct -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6}' >> pass.hg38.r3
```

Remove STRs that don't exactly lift to the detected hg38 coordinates
```
bedtools intersect -a pass.hg19.r6.lifted_to_hg38 -b pass.hg38.r3 -v | cut -f 1-4 | bedtools sort | uniq > bad_markers.hg19.bed
bedtools intersect -a pass.hg19.r6.lifted_to_hg38 -b pass.hg38.r3 -wa -wb | awk '$2 != $6 || $3 != $7 || $4 != $8' | cut -f 1-4 | bedtools sort | uniq >> bad_markers.hg19.bed
bedtools intersect -a pass.hg19.r6.lifted_to_hg38 -b pass.hg38.r3 -wa -wb | awk '$2 != $6 || $3 != $7 || $4 != $8' | cut -f 5-8 | bedtools sort | uniq > bad_markers.hg38.bed
paste pass.hg19.r6.lifted_to_hg38 pass.hg19.r6 | bedtools intersect -a - -b bad_markers.hg19.bed -v | cut -f 5-8 > pass.hg19.r7
bedtools intersect -a pass.hg38.r3 -b bad_markers.hg38.bed -v > pass.hg38.r4
```

Remove entries for STRs that overlap annotated loci
```
bedtools intersect -a pass.hg19.r7 -b all_annot_markers.hg19.fixed.bed -v > pass.hg19.r8
bedtools intersect -a pass.hg38.r4 -b all_annot_markers.hg38.fixed.bed -v > pass.hg38.r5
./liftOver -minMatch=1.0 pass.hg19.r8 hg19ToHg38.over.chain.gz pass.hg19.r8.lifted_to_hg38 unlifted
```

Perform various sanity checks on the results
i)  Verify that all of hg19's entries are shared with GRCh38
ii) Verify that the number of entries unique to GRCh38 matches the number of non-intersecting entries
```
nlines_a=`cat pass.hg19.r8.lifted_to_hg38 | wc -l`
nlines_b=`cat pass.hg38.r5 | wc -l`
nshared=`cat pass.hg19.r8.lifted_to_hg38 pass.hg38.r5 | sort | uniq -d | wc -l`
nunique=`bedtools intersect -a pass.hg38.r5 -b pass.hg19.r8.lifted_to_hg38 -v | wc -l`
if [ $nlines_a -ne $nshared ]
then
    echo "Logical error"
    exit 1
fi
nextra=`expr $nlines_b - $nlines_a`
if [ $nextra -ne $nunique ]
then
    echo "Logical error"
    exit 1
fi
```

Construct the final two references while ensuring that the numbering of STRs shared across references is consistent

hg19
```
cat pass.hg19.r8 | awk -v OFS="\t" '{print $0, ($3-$2+1)/$4, "STR_"NR}' > hg19.hipstr_reference.bed
cat all_annot_markers.hg19.fixed.bed | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, $6}' >> hg19.hipstr_reference.bed
bedtools sort -i hg19.hipstr_reference.bed > hg19.hipstr_reference.sorted.bed 
```

hg38
```
cat pass.hg19.r8.lifted_to_hg38 | awk -v OFS="\t" '{print $0, ($3-$2+1)/$4, "STR_"NR}' > hg38.hipstr_reference.bed
bedtools intersect -a pass.hg38.r5 -b pass.hg19.r8.lifted_to_hg38 -v | awk -v OFFSET=$nshared -v OFS="\t" '{print $0, ($3-$2+1)/$4, "STR_"(NR+OFFSET)}' >> hg38.hipstr_reference.bed
cat all_annot_markers.hg38.fixed.bed >> hg38.hipstr_reference.bed
bedtools sort -i hg38.hipstr_reference.bed > hg38.hipstr_reference.sorted.bed 
```

