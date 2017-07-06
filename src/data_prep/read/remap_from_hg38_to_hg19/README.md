# remap_from_hg38_to_hg19

The read data is mapped on GRCh38. Since the population reference is constructed with hg19, we want to remap the reads to hg19.


# Step-wise procedure

## 1) split the bam file into small sam files
To maximize the parallel computation, we divide the input bam file into small sam files with the following scripts:
 
```
split-bam-exec.sh
split-bam.sh
```

## 2) map each file to hg19
Map the reads to hg19. To utilize the array-job function of cluster system, one need to prepare a tsv file.

```
remap_from_hg38_to_hg19.sh
splitSam.tsv
```

