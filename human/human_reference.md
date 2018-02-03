# Building an STR reference file for HipSTR
This script describes how we built the HipSTR reference for the hg19 and hg38 versions of the human reference genome. 

It begins by using [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) to identify repeats on each chromosome and then
extensively filters and processes the resulting repeats using standard unix utilities such as [grep](http://linuxcommand.org/man_pages/grep1.html) and [awk](http://linuxcommand.org/man_pages/awk1.html),
as well as three command line utilities: [datamash](https://www.gnu.org/software/datamash/), [bedtools](http://bedtools.readthedocs.io/en/latest/) and [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver). This script also requires an executable version of *Tandem Repeats Finder* called **trf409.legacylinux64**. Please download it from the TRF link above and ensure the other utilities are installed before proceeding to follow the steps below

## Scanning the genome for STRs
First, lest's download this github repository to obtain the required scripts:
```
git clone https://github.com/HipSTR-Tool/HipSTR-references.git
cd HipSTR-references
```

Next, let's download the FASTA files for these two versions of the human reference genome:
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

# Download the liftOver chain file
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

Next, we use [Tandem Repeats Finder (TRF)](https://tandem.bu.edu/trf/trf.html) to identify repeats in each of these chromosomes.
We scan each sequence using the **run_TRF.sh** script included in this repository:
```
chmod 755 trf409.legacylinux64
chmod 755 scripts/run_TRF.sh
for chrom in $(seq 1 22) X Y
do
    echo hg19/raw_fastq/chr$chrom.fa hg19/trf_results 5
    echo hg38/raw_fastq/chr$chrom.fa hg38/trf_results 5
done | xargs -L 1 -P 30 ./scripts/run_TRF.sh
```

As STRs are classically defined as having 1-6bp repeat motifs, we filtered the TRF output to only include these repeats. The script **fix_trf_output.py** performs this filtering and also corrects TRF entries in which the repeat size is incorrectly reported:
```
for chrom in $(seq 1 22) X Y
do
    echo scripts/fix_trf_output.py hg19/trf_results/chr$chrom.fa hg19/fixed_trf_results/chr$chrom.fa
    echo scripts/fix_trf_output.py hg38/trf_results/chr$chrom.fa hg38/fixed_trf_results/chr$chrom.fa
done | xargs -L 1 -P 40 python
```

## Filtering the STRs
We now filter out STRs with very low TRF scores, effectively keeping those repeats that are much more repetitive than one would expect by random chance. 
```
files=""
for chrom in $(seq 1 22) X Y
do
    files="$files,hg19/fixed_trf_results/chr$chrom.fa"
done
files=`echo $files | sed "s/,//"`
python scripts/trf_parser.py $files > filtered_repeats.hg19.bed
bedtools sort -i filtered_repeats.hg19.bed > filtered_repeats.hg19.sorted.bed
```

```
files=""
for chrom in $(seq 1 22) X Y
do
    files="$files,hg38/fixed_trf_results/chr$chrom.fa"
done
files=`echo $files | sed "s/,//"`
python scripts/trf_parser.py $files > filtered_repeats.hg38.bed
bedtools sort -i filtered_repeats.hg38.bed > filtered_repeats.hg38.sorted.bed
```


## Merging STRs within and across assemblies


Merge overlapping STRs into single entries and filter those repeats that fail merging
```
python scripts/analyze_overlaps.py filtered_repeats.hg19.sorted.bed pass.hg19 fail.hg19
```

Remove any entries within 10bp of a failed merge region
```
bedtools window -w 10 -a pass.hg19 -b fail.hg19 -v > pass.hg19.r2
```

Extract entries that aren't within 10bp of another entry or are within 10bp of one or more entries that all share the same period
```
bedtools merge -i pass.hg19.r2 -c 4,6 -o collapse -d 10 | grep -v "," > pass.hg19.r3
bedtools merge -i pass.hg19.r2 -c 4,4,4,6 -o collapse,count_distinct,distinct,collapse -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6, $7}' | sed "s/,/\//g"  >> pass.hg19.r3
```

liftOver coordinates to hg38 and remove entries that failed to lift
```
./liftOver -minMatch=1.0 pass.hg19.r3  hg19ToHg38.over.chain.gz pass.hg19.r3.lifted_to_hg38 unlifted
grep -v ^# unlifted | cat - pass.hg19.r3 | sort -k 1,1V -k 2,2n | uniq -u > pass.hg19.r4
```

Remove an hg19 region if multiple entries lifted to the same hg38 region
```
./liftOver -minMatch=1.0 pass.hg19.r4 hg19ToHg38.over.chain.gz pass.hg19.r4.lifted_to_hg38 unlifted
paste pass.hg19.r4.lifted_to_hg38 pass.hg19.r4 | datamash -g 1,2,3,4,5 -s collapse 6,7,8,9,10 | grep -v "," | cut -f 6-10 > pass.hg19.r5
```

Remove hg19 entries that lift over to hg38 regions within 10bp of another region
```
./liftOver -minMatch=1.0 pass.hg19.r5 hg19ToHg38.over.chain.gz pass.hg19.r5.lifted_to_hg38 unlifted
paste pass.hg19.r5.lifted_to_hg38 pass.hg19.r5 | bedtools sort | bedtools merge -i - -c 6,7,8,9,10 -o collapse -d 10 | grep "," | cut -f 4- | awk -v OFS="\t" '{split($1,a,","); split($2,b,","); split($3, c, ","); split($4,d,","); split($5,e,","); for (i in a) print a[i], b[i], c[i], d[i], e[i]}' > to_remove.bed
cat pass.hg19.r5 to_remove.bed | sort -k 1,1 -k 2,2n | uniq -u > pass.hg19.r6
./liftOver -minMatch=1.0 pass.hg19.r6 hg19ToHg38.over.chain.gz pass.hg19.r6.lifted_to_hg38 unlifted
```

Perform a similar set of analyses for GRCh38, but without the complications of liftOver filtering
```
python scripts/analyze_overlaps.py filtered_repeats.hg38.sorted.bed pass.hg38 fail.hg38
bedtools window -w 10 -a pass.hg38 -b fail.hg38 -v > pass.hg38.r2
bedtools merge -i pass.hg38.r2 -c 4,6 -o collapse -d 10 | grep -v "," > pass.hg38.r3
bedtools merge -i pass.hg38.r2 -c 4,4,4,6 -o collapse,count_distinct,distinct,collapse -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6, $7}' | sed "s/,/\//g"  >> pass.hg38.r3
```

Remove STRs that don't exactly lift to the detected hg38 coordinates
```
bedtools intersect -a pass.hg19.r6.lifted_to_hg38 -b pass.hg38.r3 -v | cut -f 1-5 | bedtools sort | uniq > bad_markers.hg19.bed
bedtools intersect -a pass.hg19.r6.lifted_to_hg38 -b pass.hg38.r3 -wa -wb | awk '$2 != $7 || $3 != $8 || $4 != $9 || $5 != $10' | cut -f 1-5  | bedtools sort | uniq >> bad_markers.hg19.bed
bedtools intersect -a pass.hg19.r6.lifted_to_hg38 -b pass.hg38.r3 -wa -wb | awk '$2 != $7 || $3 != $8 || $4 != $9 || $5 != $10' | cut -f 6-10 | bedtools sort | uniq > bad_markers.hg38.bed
paste pass.hg19.r6.lifted_to_hg38 pass.hg19.r6 | bedtools intersect -a - -b bad_markers.hg19.bed -v | cut -f 6-10 > pass.hg19.r7
bedtools intersect -a pass.hg38.r3 -b bad_markers.hg38.bed -v > pass.hg38.r4
```

Remove entries for STRs that overlap annotated loci
```
bedtools intersect -a pass.hg19.r7 -b human/all_annot_markers.hg19.fixed.bed -v > pass.hg19.r8
bedtools intersect -a pass.hg38.r4 -b human/all_annot_markers.hg38.fixed.bed -v > pass.hg38.r5
./liftOver -minMatch=1.0 pass.hg19.r8 hg19ToHg38.over.chain.gz pass.hg19.r8.lifted_to_hg38 unlifted
```

Perform various sanity checks on the results
 1.  Verify that all of hg19's entries are shared with GRCh38
 2. Verify that the number of entries unique to GRCh38 matches the number of non-intersecting entries
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

Construct the final two references while ensuring that the names of STRs shared across references is consistent

hg19
```
cat pass.hg19.r8 | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "Human_STR_"NR, $5}' > hg19.hipstr_reference.bed
cat human/all_annot_markers.hg19.fixed.bed | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, $6, "N/A"}' >> hg19.hipstr_reference.bed
bedtools sort -i hg19.hipstr_reference.bed > hg19.hipstr_reference.sorted.bed 
```

hg38
```
cat pass.hg19.r8.lifted_to_hg38 | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "Human_STR_"NR, $5}' > hg38.hipstr_reference.bed
bedtools intersect -a pass.hg38.r5 -b pass.hg19.r8.lifted_to_hg38 -v | awk -v OFFSET=$nshared -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "Human_STR_"(NR+OFFSET), $5}' >> hg38.hipstr_reference.bed
cat human/all_annot_markers.hg38.fixed.bed | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, $6, "N/A"}' >> hg38.hipstr_reference.bed
bedtools sort -i hg38.hipstr_reference.bed > hg38.hipstr_reference.sorted.bed 
```

Generate the final compressed human references for HipSTR
```
cat hg19.hipstr_reference.sorted.bed                  | gzip -c >   hg19.hipstr_reference.bed.gz
cat hg19.hipstr_reference.sorted.bed | sed "s/^chr//" | gzip -c > GRCh37.hipstr_reference.bed.gz
cat hg38.hipstr_reference.sorted.bed                  | gzip -c >   hg38.hipstr_reference.bed.gz
cat hg38.hipstr_reference.sorted.bed | sed "s/^chr//" | gzip -c > GRCh38.hipstr_reference.bed.gz
```

Delete the temporary files created along the way
```
rm -r hg19 hg38
rm pass.hg19* pass.hg38* fail.hg19 fail.hg38 bad_markers.hg19.bed bad_markers.hg38.bed to_remove.bed
rm filtered_repeats.hg*
rm hg19.hipstr_reference.bed hg38.hipstr_reference.bed unlifted
rm hg19.hipstr_reference.sorted.bed hg38.hipstr_reference.sorted.bed
```
