# make_metaplot

## Purpose 

Generate metaplot for NGS data


## Getting Started

```
# 1. Generate template config file
$ python make_metaplot.py -c template.yaml -t

# 2. Modify the following data
# region_list
# bam_list
# bw_list
# bw_fwd_list, bw_rev_list 
# matrix

# 3. Run program
$ python make_metaplot.py -c template.yaml
```



## How-to

1. bamCoverage

2. computeMatrix

3. plotProfile


Strand-specific: 

sense: bw_fwd + bed_fwd; bw_rev + bed_rev 
anti:  bw_fwd + bed_rev; bw_rev + bed_fwd


## Dependencies

[deeptools](https://github.com/deeptools/deepTools) - using the submodules: `bamCoverage`, `computeMatrix` and `plotProfile` to create metaplot file

[pysam](https://github.com/pysam-developers/pysam) - reading and processing SAM files 

[xopen](https://github.com/pycompression/xopen) - read/write gzipped files

[pyBigwig](https://github.com/deeptools/pyBigWig) - quick access to bigWig and bigBed files

[json](https://docs.python.org/3/library/json.html) - JSON encoder and decoder

[yaml](https://pyyaml.org/wiki/PyYAMLDocumentation) - a YAML parser and emitter for Python

[toml](https://github.com/uiri/toml) - a TOML parser and emitter for Python (not required)


## Alternative tools 

[ngs.plot.r](https://github.com/shenlab-sinai/ngsplot) - a convinent tool for the same purpose. lack support for custome BED files, strand-specific libraries.