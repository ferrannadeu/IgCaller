# IgCaller

#### Reconstructing immunoglobulin gene rearrangements and oncogenic translocations from WGS

IgCaller is a python program designed to fully characterize the immunoglobulin gene rearrangements and oncogenic translocations from whole-genome sequencing in lymphoid neoplasms. For more information: https://rdcu.be/b5tsw.

![alt text](https://github.com/ferrannadeu/IgCaller/blob/master/IgCaller_workflow.jpg)

### Requirements

IgCaller is based on python3 and requires the following modules: statistics, regex (v2.5.29 and v2.5.30), argparse (v1.1), numpy (1.16.2 and v1.16.3), and scipy (v1.2.1 and v1.3.0). Although providing the versions of the modules tested, we are not aware about any specific version requirement for running IgCaller. Other modules used by IgCaller but already included in base python are: subprocess, sys, os, itertools, operator, collections, gzip, pickle.   

The only required non-python program is [samtools](http://www.htslib.org) (v1.6 and v1.9 have been tested).

### Download

Download and uncompress the ZIP file OR clone the IgCaller repository:

```
git clone https://github.com/ferrannadeu/IgCaller
```

### Running IgCaller

#### Mandatory arguments:
*	inputsFolder (-I): path to the folder containing the supplied IgCaller reference files.
*	genomeVersion (-V): version of the reference human genome used when aligning the WGS data (hg19 or hg38).
*	chromosomeAnnotation (-C): chromosome annotation [ensembl = without 'chr' (i.e. 1); ucsc = with 'chr' (i.e. chr1)].
*	bamT (-T): path to tumor BAM file aligned using [BWA-MEM](http://bio-bwa.sourceforge.net/) (v0.7.15 and v0.7.17 have been tested).
*	bamN (-N): path to normal BAM file aligned using BWA-MEM, if available.
*	refGenome (-R): path to reference genome FASTA file (not mandatory, but recommended, when specifying a normal BAM file. Mandatory when bamN is not specified).

#### Optional arguments:
*	pathToSamtools (-ptsam): path to the directory where samtools is installed. There is no need to specify it if samtools is found in “PATH” (default = ‘empty’, assuming it is in “PATH”).
*	outputPath (-o): path to the directory where the output should be stored. Inside the defined directory IgCaller will automatically create a folder named “tumorSample”_IgCaller where output files will be saved (default, current working directory).
*	mappingQuality (-mq): mapping quality cut off to filter out reads for IG V(D)J reconstruction (default = 0).
*	baseQuality (-bq): base quality cut off to consider a position in samtools mpileup when reconstructing both normal and tumor sequences (default = 13).
*	minDepth (-d): depth cut off to consider a position (default = 1).
*	minAltDepth (-ad): alternate depth cut off to consider a potential alternate nucleotide (default = 1)
*	vafCutoffNormal (-vafN): minimum variant allele frequency (VAF) to consider a nucleotide when reconstructing the germ line sequence using the supplied normal BAM file (if available) (default = 0.2).
*	vafCutoff (-vaf): minimum VAF to consider a nucleotide when reconstructing the tumor sequence (default = 0.1). Try to increase this value if only unproductive rearrangements are found due to stop codons. We have observed that relatively high coverage WGS (i.e. 100x) might carry many variants (likely sequencing artifacts) at VAFs around 10-20%.
*	tumorPurity (-p): purity of the tumor sample (i.e. tumor cell content) (default = 1). It is used to adjust the VAF of the mutations found in the tumor BAM file before filtering them using the vafCutoff, to adjust the score of each rearrangement, and to adjust the reduction of read depth in the CSR analysis.
*	minNumberReadsTumorOncoIg (-mntonco): minimum score supporting an IG rearrangement in order to be annotated (default = 4).
*	minNumberReadsTumorOncoIgPass (-mntoncoPass): minimum score supporting an IG rearrangement in the tumor sample in order to be considered as high confidence (default = 10).
*	maxNumberReadsNormalOncoIg (-mnnonco): maximum number of reads supporting an IG rearrangement in the normal sample in order to be considered as high confidence (default = 2).
*	mappingQualityOncoIg (-mqOnco): mapping quality cut off to filter out reads when analyzing oncogenic IG rearrangements (default = 15).
*	numThreads (-@): maximum number of threads to be used by samtools (default = 1).
* keepMiniIgBams (-kmb): should IgCaller keep (i.e. no remove) mini IG BAM files used in the analysis? (default = no).
* sequencing (-seq): sequencing technique (whole-genome sequencing (wgs) or whole-exome sequencing (wes)) (default = wgs).


#### Basic command line:
```
python3 path/to/IgCaller/IgCaller_v1.1.py -I path/to/IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T path/to/bams/tumor.bam -N path/to/bams/normal.bam -R path/to/reference/genome_hg19.fa -o path/to/IgCaller/outputs/
```

#### Tested on:
IgCaller was tested on a MacBook Pro (macOS Mojave), Ubuntu (16.04 and 18.04), and MareNostrum 4 (Barcelona Supercomputing Center, SUSE Linux Enterpirse Server 12 SP2 with python/3.6.1).

#### Running time:
IgCaller only requires 1 CPU, and it usually takes <2-5 minutes to characterize the complete immunoglobulin gene of one tumor sample. Execution time increases with sequencing coverage and number of potential rearrangements present in the sample.

### Demo dataset
A demo dataset to test IgCaller is provided under the "Demo" folder. If you downloaded IgCaller using the download ZIP file button, the following bam files should be manually downloaded and placed in their respective folders:
* Demo/BAM/[1344-01-01ND.bam](https://github.com/ferrannadeu/IgCaller/blob/master/Demo/BAM/1344-01-01ND.bam)
* Demo/BAM/[1344-01-03TD.bam](https://github.com/ferrannadeu/IgCaller/blob/master/Demo/BAM/1344-01-03TD.bam)

These files are stored on Git LFS and they are not included in the downloaded ZIP file (only a pointer is downloaded) but they are downloaded when "git clone"ing the repository.

### Outputs

IgCaller returns a set of tab-separated files:

*	tumor_sample_output_filtered.tsv: High confidence rearrangements passing the defined filters.
*	tumor_sample_output_IGH.tsv: File containing all IGH rearrangements.
*	tumor_sample_output_IGK.tsv: File containing all IGK rearrangements.
*	tumor_sample_output_IGL.tsv: File containing all IGL rearrangements.
*	tumor_sample_output_class_switch.tsv: File containing all CSR rearrangements.
*	tumor_sample_output_oncogenic_IG_rearrangements.tsv: File containing all oncogenic IG rearrangements (translocations, deletions, inversions, and gains) identified genome-wide.

### Other notes

An R script to help the study of mutational signatures in CLL is available under the "Mutational_signature_analysis_in_CLL" folder. This script aims to determine the presence/absence of non-canonical AID mutations (signature 9) in CLL patients using an already defined catalogue of single nucleotide variants.

### Citation

If you use IgCaller, please cite:

Nadeu, F., Mas-de-les-Valls, R., Navarro, A. et al. IgCaller for reconstructing immunoglobulin gene rearrangements and oncogenic translocations from whole-genome sequencing in lymphoid neoplasms. Nature Communications 11, 3390 (2020). https://doi.org/10.1038/s41467-020-17095-7.

### Contact

Bugs, comments and improvements can be send to *nadeu@clinic.cat*. They will be very much appreciated!

### Releases

* Current master branch (v1.1): 
  * Added the annotation whether the non-IG breakpoint of the IG rearrangements identified map within interspersed repeats and low complexity DNA sequences according to the [RepeatMasker UCSC track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk). repName, repClass and	repFamily are annotated. Note that most false positive translocations initially called by IgCaller map within simple repeats. This flag should allow to easily clean up the list of oncogenic IG translocations reported without spending too much time doing manual validation. 
  * Added the annotation of the gene closest to the non-IG breakpoint of the oncogenic IG rearrangements identified based on RefSeq annotations (maximum upstream and downstream distance considered: 250kb).
  * Minor bug that made IgCaller crash in 1 out of the >1,000 samples tested is now fixed.
  * Added compatibility for BAM files obtained from whole-exome sequencing (see optional argument --sequencing/-seq). The functionality of IgCaller is the same for WGS and WES samples with only a small difference in the pre-defined filtering step.
* v1.0:
  * First version of IgCaller as described in the original manuscript.
