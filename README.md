# HipSTR-references
Reference files for running [HipSTR](https://hipstr-tool.github.io/HipSTR/) with various organisms

## Getting Started
HipSTR requires a BED file containing the STR regions to genotype. We've generated these files for various organisms and assemblies as listed in the table below:

| Organism | Assembly | URL  |
| :------  | :------  | :--- |
| Human    | hg19     |  https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg19.hipstr_reference.bed.gz |
| Human    | hg38     |  https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg38.hipstr_reference.bed.gz |
| Human    | GRCh37   |  https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/GRCh37.hipstr_reference.bed.gz |
| Human    | GRCh38   |  https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/GRCh38.hipstr_reference.bed.gz |
| Mouse    | mm10     |  https://github.com/HipSTR-Tool/HipSTR-references/raw/master/mouse/mm10.hipstr_reference.bed.gz |

To use a particular file with HipSTR, download it and decompress it using gzip. For instance, if you are interested in the BED file for **humans** for the **hg19** assembly:
    
    wget https://github.com/HipSTR-Tool/HipSTR-references/raw/master/human/hg19.hipstr_reference.bed.gz
    gunzip hg19.hipstr_reference.bed.gz

The file *hg19.hipstr_reference.bed* is then ready to use as input to HipSTR's **--regions** option.

## Building the references
For both mouse and human, we've provided detailed information about how we built each reference in this repository.
If you're interested in learning more, please checkout the **human** or **mouse** subdirectories.

## The human HipSTR reference

## Other organisms
If you're interested in a file for a different model organism, send us an email at hipstrtool@gmail.com and we'll be happy to help build one

