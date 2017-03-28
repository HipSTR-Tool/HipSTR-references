# Building an STR reference file for HipSTR
Here we describe how we built the HipSTR reference for the GRCh37 and GRCh38 versions of the human reference genome. We begin by scanning each
version of the reference genome for repeats using [Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html). We then
extensively filter and process the resulting repeats using standard unix utilities such as [grep](http://linuxcommand.org/man_pages/grep1.html) and [awk](http://linuxcommand.org/man_pages/awk1.html), 
as well as three command line utilities: [datamash](https://www.gnu.org/software/datamash/), [bedtools](http://bedtools.readthedocs.io/en/latest/) and [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver).




## Scanning the genome for STRs
First, let's download the relevant reference files for these two version of the human reference genome:
```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

## Filtering and merging the STRs
[Tandem Repeats Finder](https://tandem.bu.edu/trf/trf.html) frequently reports multiple repeats that overlap one another.
