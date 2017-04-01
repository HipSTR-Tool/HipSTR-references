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

Next, we used [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) to identify repeats in each of these chromosomes.
We first downloaded an executable version of this program called **trf409.legacylinux64**. We then scanned each sequence using the **run_TRF.sh** script included in this repository as follows:
```
for chrom in $(seq 1 22) X Y
do
    echo hg19/raw_fastq/chr$chrom.fa hg19/trf_results 5
    echo hg38/raw_fastq/chr$chrom.fa hg38/trf_results 5
done | xargs -L 1 -P 30 ./run_TRF.sh
```



## Filtering and merging the STRs
