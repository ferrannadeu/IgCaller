### Demo dataset to test IgCaller

Tumor and normal BAM files of the IG loci of a CLL patient might be found in the "BAM" folder. 
The expected output is placed in the "Outputs" folder. 

In order to generate this demo, IgCaller was run using the following command line:

```
python3 IgCaller_v1.py -I IgCaller_reference_files/ -N BAM/1344-01-01ND.bam -T BAM/1344-01-03TD.bam -V hg19 -C ensembl
```
