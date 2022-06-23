# IgCaller

#### Reconstructing immunoglobulin gene rearrangements and oncogenic translocations from WGS, WES, and capture NGS data

IgCaller is a python program designed to fully characterize the immunoglobulin gene rearrangements and oncogenic translocations in lymphoid neoplasms. It was originally developed to work with WGS data but it has been extended to work with WES and high-coverage, capture-based NGS data. For more information, please read the original [publication](https://rdcu.be/b5tsw).

![alt text](https://github.com/ferrannadeu/IgCaller/blob/master/IgCaller_workflow.jpg)

### Requirements

IgCaller is based on python3 and requires the following modules: statistics, regex (v2.5.29 and v2.5.30), argparse (v1.1), numpy (1.16.2 and v1.16.3), and scipy (v1.2.1 and v1.3.0). Although providing the versions of the modules tested, we are not aware about any specific version requirement for running IgCaller. Other modules used by IgCaller but already included in base python are: subprocess, sys, os, itertools, operator, collections, gzip, pickle.   

The only required non-python program is [samtools](http://www.htslib.org) (v1.6 and v1.9 have been tested).

### Installation

Download and uncompress the ZIP file or clone the IgCaller repository:

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
* sequencing (-seq): sequencing technique (whole-genome sequencing (wgs), whole-exome sequencing (wes), or high coverage NGS (capture)) (default = wgs).
*	outputPath (-o): path to the directory where the output should be stored. Inside the defined directory IgCaller will automatically create a folder named “tumorSample”_IgCaller where output files will be saved (default, current working directory).
*	tumorPurity (-p): purity of the tumor sample (i.e. tumor cell content) (default = 1). It is used to adjust the VAF of the mutations found in the tumor BAM file before filtering them using the vafCutoff, to adjust the score of each rearrangement, and to adjust the reduction of read depth in the CSR analysis.
* keepMiniIgBams (-kmb): should IgCaller keep (i.e. no remove) mini IG BAM files used in the analysis? (default = no).
*	pathToSamtools (-ptsam): path to the directory where samtools is installed. There is no need to specify it if samtools is found in “PATH” (default = ‘empty’, assuming it is in “PATH”).
*	numThreads (-@): maximum number of threads to be used by samtools (default = 1).
*	mappingQuality (-mq): mapping quality cut off to filter out reads for IG V(D)J reconstruction (default = 0).
*	baseQuality (-bq): base quality cut off to consider a position in samtools mpileup when reconstructing both normal and tumor sequences (default = 13).
*	minDepth (-d): depth cut off to consider a position (default = 2).
*	minAltDepth (-ad): alternate depth cut off to consider a potential alternate nucleotide (default = 1)
*	vafCutoffNormal (-vafN): minimum variant allele frequency (VAF) to consider a nucleotide when reconstructing the germ line sequence using the supplied normal BAM file (if available) (default = 0.2).
*	vafCutoff (-vaf): minimum VAF to consider a nucleotide when reconstructing the tumor sequence when working with phased reads and without phased reads, respectively (two numbers between 0-1 separated by a comma, default=0.66,0.1). Try to increase these values if only unproductive rearrangements are found due to stop codons. We have observed that relatively high coverage WGS data (i.e. 100x) might carry many variants (likely sequencing artifacts) at VAFs around 10-20%.
* scoreCutoff (-s): minimum score supporting a gene rearrangement in order to be considered as high confidence (default will consider 3 for wgs/wes and 15 for capture. Usase: -s X where X is any number).
* shortReportedVseq (-shortV): make V sequence start at approx FR1 (yes/no, default=yes).
*	mappingQualityOncoIg (-mqOnco): mapping quality cut off to filter out reads when analyzing oncogenic IG rearrangements (default = 15).
*	minNumberReadsTumorOncoIg (-mntonco): minimum score supporting an IG rearrangement in order to be annotated (default = 5).
*	minNumberReadsTumorOncoIgPass (-mntoncoPass): minimum score supporting an IG rearrangement in the tumor sample in order to be considered as high confidence (default = 10).
*	maxNumberReadsNormalOncoIg (-mnnonco): maximum number of reads supporting an IG rearrangement in the normal sample in order to be considered as high confidence (default = 2).
* maxNumberCountInPoN (-mncPoN): maximum number of count in panel of normals (PoN) in order to be considered as high confidence (default=2).

#### Basic command line:
```
python3 path/to/IgCaller/IgCaller.py -I path/to/IgCaller/IgCaller_reference_files/ -V hg19 -C ensembl -T path/to/bams/tumor.bam -N path/to/bams/normal.bam -R path/to/reference/genome_hg19.fa -o path/to/IgCaller/outputs/
```

#### Tested on:
IgCaller was tested on a MacBook Pro (macOS Mojave and Big Sur), Ubuntu (16.04 and 18.04), and MareNostrum 4 (Barcelona Supercomputing Center, SUSE Linux Enterpirse Server 12 SP2 with python/3.6.1).

#### Running time:
IgCaller only requires 1 CPU, and it usually takes around 5 minutes to characterize the complete immunoglobulin gene of one tumor sample. Execution time increases with sequencing coverage and number of potential rearrangements present in the sample.

### Demo dataset
A demo dataset to test IgCaller is provided under the "Demo" folder.

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

Bugs, comments and improvements can be send to *nadeu@clinic.cat*.

### Releases
* Master branch (v1.2.1):
  * Improved phasing of mutations within the rearranged allele. When reads spanning the J-V rearrangement are not covering a fraction of the V gene, IgCaller now tries to phase reads on the fly based on the previous mutations identified in the V gene. This new functionality improves the capacity to phase mutations and SNPs while increases the sensitivity to detect mutations specially in heavily mutated V genes.
  * Corrected internal split-read annotation that made IgCaller crash in a subset of samples aligned to a reference genome that included "alt" contigs. See [issue #5](https://github.com/ferrannadeu/IgCaller/issues/5) for details. 
* v1.2:
  * General improvements (major):
    * Added compatibility for reference genomes containing lowercase nucleotides (lowercase nucleotides are automatically converted to uppercase in the fly).
    * Better tolerance of ambiguous nucleotides in the sequences.
    * Small imperfections when counting the nucleotides found in each genomic position during the reconstruction of the consensus V(D)J sequence have been solved/improved.
    * Increased sensitivity to detect V(D)J rearrangements by expanding the J gene region used to consider ‘insert-size’ reads spanning a given J-V rearrangement. Due to the small size of the J gene and insert size of the NGS libraries, some reads spanning a J-V rearrangement might align 5’ of the J gene. These reads are now recovered. Besides, shorter split reads are kept and tried to ‘re-align’ than in the previous version.
    * Improved (and simplified) pre-defined filtering of high-quality gene rearrangements. High-quality rearrangements found in the “output_filtered.tsv” file are now considered based on the score of the rearrangement, its mapping quality (new), mechanisms [deletion, inversion] (new), CDR3 sequence (new), and V(D)J genes annotated. As in the original version, it assumes a clonal/oligoclonal situation where the same CDR3 sequence and/or combination of V(D)J genes is not possible/real but caused by an inaccurate mapping of the reads leading to multiple solutions.
  * Phasing of mutations/SNPs within the V(D)J rearranged sequence:
    * IgCaller now considers only those reads spanning the J-V rearrangement of interest to better define mutations and SNPs present in the rearranged sequence/allele. If less than ‘minDepth’ reads spanning the J-V rearrangement are found in a given position, this position is analyzed considering all reads covering it, as performed in the previous versions. Note that a column with the details of the “phasing” step has been added (column labeled “V-PhaseInfo”). This column looks like “a/b - c”, where a,b,c are numbers: a = number of mutations phased, b = total number of mutations found in V gene, c = number of SNPs not phased.
    * 'vafCutoff / -vaf' default set to 0.66,0.1 for phased and unphased mutations/SNPs, respectively. This means that a minimum VAF of 0.66 is needed to consider a mutation when a position is analyzed using only “phased reads” (i.e. those spanning the J-V rearrangement), while a VAF of 0.1 is needed if a position is analyzed using all reads covering the position (i.e. “unphased reads”). 
    * Added ‘--shortReportedVseq / -shortV’ to make the V sequence start approx. at FR1 (active by default). Useful when interested only on the FR1-FR3 sequence so the reported number of phased mutations and SNPs is related to this shorter V gene sequence of interest.
    * The coordiantes of some V genes have been manually curated for a proper functionality of the “-shortV” method.
  * Genome-wide oncogenic IG rearrangements:
    * Incoherent filtering of split-read support when analyzing genome-wide IG translocations has been corrected (see issue #4: https://github.com/ferrannadeu/IgCaller/issues/4). It increases the sensitivity to detect IG oncogenic events (translocations, deletions, etc.).
    * The IGL locus considered when analyzing IG translocations has been slightly expanded.
  * General improvements (minor):
    * Added ‘-seq capture’ mode: added to improve compatibility of IgCaller with high coverage, capture-based NGS data (depth of >250/500x). By default, a minimum score of 15 is needed to keep a gene rearrangement as “high quality” in the “output_filtered.tsv”. Besides, to speed up the execution, under this ‘-seq capture’ mode, IgCaller removes in the first steps of the analysis those J-V rearrangements supported by less than two split reads (note that these rearrangements will not be included neither in the locus-specific output file). If the latter functionality is not desired, IgCaller could be run as “-seq wes” to analyze high coverage data, although it might remarkably increase the execution time and number of low-evidence (likely artefactual) rearrangements passing the default parameters.
    * Added '--scoreCutoff / -s': minimum score supporting a gene rearrangement in order to be considered as high confidence (by default will consider 3 for wgs/wes and 15 for capture). It might be useful when the default parameters are to permissive/restrictive. It should be adjusted based on the quality of the data and sequencing coverage.
    * ‘--minDepth / -d’ default set to 2 (it was set to 1 in the previous versions). A minDepth of 2 is necessary for a proper analysis considering “phased reads”.
    * Removed version number from the name of the scripts. It has been added within the scripts and accessible using ‘python3 IgCaller.py -v’.
    * Minor bug that made IgCaller crash in a specific scenario when annotating the gene near the breakpoints of the translocations is now fixed.
    * Other minor improvements related to strand orientation, filter of artefactual/really short V gene rearrangements, etc. have been implemented.
* v1.1: 
  * Added a panel of normals (PoN) with the oncogenic IG rearrangements found by IgCaller in normal WGS data, which can be considered as sequencing artifacts and filtered out from the tumoral samples. When filtering using these PoN, the breakpoint in the IG locus is not considered based on its position but considering the IG locus as a whole (i.e. if the break occurs in the IGH, IGK or IGL locus irrespectively of the exact position), while the exact position of the non-IG breakpoint is considered (with a +/- 1,000 bp window). The normal samples analyzed in the IgCaller manuscript were used to build a PoN for hg19 (n=243 samples) and hg38 (n=161). See optional argument --maxNumberCountInPoN (-mncPoN) to adjust this filter.
  * Added the annotation whether the non-IG breakpoint of the IG rearrangements identified map within interspersed repeats and low complexity sequences according to the [RepeatMasker UCSC track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=rmsk). repName, repClass and	repFamily are annotated.
  * Added the annotation of the gene closest to the non-IG breakpoint of the oncogenic IG rearrangements identified based on RefSeq annotations (maximum upstream and downstream distance considered: 250 kb).
  * Minor bug that made IgCaller crash in 1 out of the >1,000 samples tested is now fixed.
  * Added compatibility for BAM files obtained from whole-exome sequencing (see optional argument --sequencing/-seq). The functionality of IgCaller is the same for WGS and WES samples with only a small difference in the pre-defined filtering step.  
* v1.0:
  * First version of IgCaller as described in the original manuscript.
