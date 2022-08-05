[![DOI](https://zenodo.org/badge/140651301.svg)](https://zenodo.org/badge/latestdoi/140651301)


# MolecularMethodsInGenomeResearch
This directory contains teaching materials for a practical course about lab methods in genome research. Besides a set of slides for presentation in class, there are scripts which are needed to perform basic bioinformatic analysis. The concept behind this course and additonal details about the content are included in a publication: https://doi.org/10.1515/jib-2019-0005


<img alt="Molecular Methods in Genome Research course content overview" src="https://github.com/bpucker/figures/blob/main/MolecularMethodsInGenomeReserachCourseOverview.png" width="50%" height="50%">


Feel free to use any of the provided materials in your own courses.

## Extraction of sequences for primer design
This script allows the extraction of a region of interest for primer design and other applications.

```
Usage
python3 seqex3.py --in <FILE> --out <FILE> --contig <STR> --start <INT> --end <INT>

Mandatory:
--in       STR   Input FASTA file.
--out      STR   Output FASTA file.
--contig   STR   Sequence ID
--start    INT   Start position
--end      INT   End position
```

`--in` specifies the input FASTA file for the sequence extraction. Sequence IDs will be split at the first space.

`--out` specifies the output FASTA file. Extracted sequence parts will be stored in this file.

`--contig` specifies the sequence ID of a target sequence. A part of this sequence will be extracted.

`--start` specifies the start position of the region of interest.

`--end` specifies the end position of the region of interest.


# References:

https://doi.org/10.1371/journal.pone.0164321

https://www.biorxiv.org/content/early/2018/09/06/407627

http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030196

