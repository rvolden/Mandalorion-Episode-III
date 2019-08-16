# Mandalorion-Episode-III #
*Revenge of the Isoforms*

Takes R2C2/C3POa data and defines high confidence isoforms.

## Dependencies ##

- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)
- [emtrey](https://github.com/rvolden/emtrey) ([go](https://golang.org/dl/))
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)

The paths to these will need to be put into your config file [like this.](example_config) If you have the program installed or in your path already, replace the path with the name of the program.

## Usage ##
```bash
python3 defineAndQuantifyWrapper.py [OPTIONS]
```

Required options:
```
-c  config file containing paths to required dependencies (above)
-p  output path
-g  annotation file (gtf)
-G  genome file (fasta)
-a  adapter file (fasta)
-f  R2C2 read file (fasta)
-b  R2C2 subread file (fastq)
```


Tweakable parameters:
```
-m  score matrix file (defaults to NUC.4.4.mat)
-u  upstream buffer, defines leniency window for TSS definition (default 10)
-d  downstream buffer, defines leniency window for polyA definition (default 50)
-s  subsample consensus, defines how many random subreads are used to make isoforms (default 500)
-r  minimum ratio, percentage of reads that align to a locus required for an isoform (default 0.05)
-i  minimum internal ratio (default 0.125)
-R  minimum number of reads for an isoform (default 5)
-O  overhangs, defines bounds for unaligned bases on ends, format: min5',max5',min3',max3' (default 0,40,0,40)
-t  number of threads to use for minimap2 (default 4)
```


