# Mandalorion-Episode-III #
*Revenge of the Isoforms*

Takes R2C2/C3POa data and defines high confidence isoforms.

## Dependencies ##

- [minimap2](https://github.com/lh3/minimap2)
- [racon](https://github.com/isovic/racon)
- [emtrey](https://github.com/rvolden/emtrey) ([go](https://golang.org/dl/))
- [blat source](https://users.soe.ucsc.edu/~kent/src/blatSrc35.zip) or [blat executable](http://hgdownload.soe.ucsc.edu/admin/exe/)

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
