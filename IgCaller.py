# Modules
import argparse
import subprocess
import sys

# Manage inputs
parser  = argparse.ArgumentParser(prog='IgCaller', description='''IgCaller v1.2.1 (https://github.com/ferrannadeu/IgCaller)''')

parser.add_argument('-v', '--version', action='version', version="%(prog)s v1.2.1")

parser.add_argument('-I', '--inputsFolder',
					dest = "inputsFolder",
					action = "store",
					required=True,
					help = "Path to folder containing IgCaller reference files") 

parser.add_argument('-V', '--genomeVersion',
					dest = "genomeVersion",
					action = "store",
					choices=['hg19', 'hg38'],
					default = None,
					required=True,
					help = "Reference genome version [hg19, hg38]") 

parser.add_argument('-C', '--chromosomeAnnotation',
					dest = "chrAnnot",
					action = "store",
					choices=['ensembl', 'ucsc'],
					default = None,
					required=True,
					help = "Chromosome annotation [ensembl = without 'chr' (i.e. 1); ucsc = with 'chr' (i.e. chr1)]") 

parser.add_argument('-T', '--bamT',
					dest = "bamT",
					action = "store",
					required=True,
					help = "Tumor bam file") 

parser.add_argument('-N', '--bamN',
					dest = "bamN",
					action = "store",
					default = None,
					help = "Normal bam file, if avilable") 

parser.add_argument('-R', '--refGenome',
					dest = "refGenome",
					action = "store",
					default=None,
					help = "Reference genome fasta") 

parser.add_argument('-seq', '--sequencing', 
					dest = "seq",
					action = "store",
					default = "wgs",
					choices=['wgs', 'wes', 'capture'],
					help = "Sequencing technique [wgs/wes/capture, default=wgs]")

parser.add_argument('-o', '--outputPath',
					dest = "outputPath",
					action = "store",
					default = None,
					help = "A folder inside this directory will be created with the output [default = current directory]") 

parser.add_argument('-p', '--tumorPurity', 
					dest = "tumorPurity",
					action = "store",
					default = "1",
					help = "Purity of the tumor sample (tumor cell contect) [0-1, default=1]") 

parser.add_argument('-kmb', '--keepMiniIgBams', 
					dest = "keepMiniIgBams",
					action = "store",
					default = "no",
					help = "Should IgCaller keep (ie no remove) mini IG BAM files used in the analysis? [yes/no, default=no]")

parser.add_argument('-ptsam', '--pathToSamtools',
					dest = "pathToSamtools",
					action = "store",
					default = "",
					help = "Path to the directory where samtools is installed. No need to specify it if samtools is found in 'PATH'") 

parser.add_argument('-@', '--threadsForSamtools', 
					dest = "threadsForSamtools",
					action = "store",
					default = "1",
					help = "Maximum number of threads used for samtools [default=1]")

parser.add_argument('-mq', '--mappingQuality',
					dest = "mapq",
					action = "store",
					default = "0",
					help = "Mapping quality cut off to filter out reads for IG reconstruction [default=0]") 

parser.add_argument('-bq', '--baseQuality',
					dest = "baseq",
					action = "store",
					default = "13",
					help = "Base quality cut off for samtools mpileup for mutation analysis [default=13]") 

parser.add_argument('-d', '--minDepth',
					dest = "depth",
					action = "store",
					default = "2",
					help = "Depth cut off to consider a position [0-inf, default=2]")
					
parser.add_argument('-ad', '--minAltDepth',
					dest = "altDepth",
					action = "store",
					default = "1",
					help = "Alt depth cut off to consider a nucleotide [0-inf, default=1]")
					
parser.add_argument('-vafN', '--vafCutoffNormal',
					dest = "vafCutoffNormal",
					action = "store",
					default = "0.20",
					help = "VAF cut off to consider a mutation [0-1, default=0.20]")

parser.add_argument('-vaf', '--vafCutoff',
					dest = "vafCutoff",
					action = "store",
					default = "0.66,0.1",
					help = "VAF cut off to consider a mutation when working with phased reads and without phased reads, respectively [two numbers between 0-1 separated by a comma, default=0.66,0.1]")	

parser.add_argument('-s', '--scoreCutoff',
					dest = "scoreCutoff",
					action = "store",
					default = "",
					help = "Minimum score supporting a gene rearrangement in order to be considered as high confidence [default='empty'; will consider 3 for wgs/wes and 15 for capture]")

parser.add_argument('-shortV', '--shortReportedVseq', 
					dest = "shortV",
					action = "store",
					default = "yes",
					choices=['yes', 'no'],
					help = "Make V sequence start at (approx) FR1 [yes/no, default=yes]")

parser.add_argument('-mqOnco', '--mappingQualityOncoIg',
					dest = "mapqOnco",
					action = "store",
					default = "15",
					help = "Mapping quality cut off to filter out reads when analyzing oncogenic IG rearrangements [default=15]") 

parser.add_argument('-mntonco', '--minNumberReadsTumorOncoIg', 
					dest = "mntonco",
					action = "store",
					default = "5",
					help = "Minimum score supporting an oncogenic IG rearrangement in order to be annotated [default=5]")

parser.add_argument('-mntoncoPass', '--minNumberReadsTumorOncoIgPass', 
					dest = "mntoncoPass",
					action = "store",
					default = "10",
					help = "Minimum score supporting an oncogenic IG rearrangement in order to be considered as high confidence [default=10]")

parser.add_argument('-mnnonco', '--maxNumberReadsNormalOncoIg', 
					dest = "mnnonco",
					action = "store",
					default = "2",
					help = "Maximum number of reads supporting an oncogenic IG rearrangement in the normal sample in order to be considered as high confidence [default=2]")

parser.add_argument('-mncPoN', '--maxNumberCountInPoN', 
					dest = "mncPoN",
					action = "store",
					default = "2",
					help = "Maximum number of count in panel of normals (PoN) in order to be considered as high confidence [default=2]")

options = parser.parse_args()

inputsFolder = options.inputsFolder
genomeVersion = options.genomeVersion
chrAnnot = options.chrAnnot
originalBamT = options.bamT
originalBamN = options.bamN
refGenome = options.refGenome
pathToSamtools = options.pathToSamtools if options.pathToSamtools == "" else options.pathToSamtools+"/"
outputPath = options.outputPath
mapq = options.mapq
baseq = options.baseq
depth = int(options.depth)
altDepth = int(options.altDepth)
vafCutoffNormal = float(options.vafCutoffNormal)
vafCutoff = options.vafCutoff.rstrip(",").split(",")
if len(vafCutoff) != 2: 
	sys.exit("IgCaller: error message... -vaf / --vafCutoff should be of the form '0.66,0.1'.")
scoreCutoff = options.scoreCutoff
tumorPurity = float(options.tumorPurity)
mntonco = int(options.mntonco)
mntoncoPass = int(options.mntoncoPass)
mnnonco = int(options.mnnonco)
mncPoN = int(options.mncPoN)
mapqOnco = options.mapqOnco
threadsForSamtools = options.threadsForSamtools
keepMiniIgBams = options.keepMiniIgBams
seq = options.seq
shortV = options.shortV

# 0) Prepare some variables and files:
print("IgCaller: checking inputs and outputs...")
## Import IgCaller defined dictionaries and functions from IgCaller_reference_files:
sys.path.insert(0, inputsFolder)
from IgCaller_functions import *

## hg19 or hg38 regions of the IG loci:
if chrAnnot == "ensembl": chrom = ""
else: chrom = "chr"
if genomeVersion == "hg19": coordsToSubset = chrom+"14:106052774-107288051 "+chrom+"22:22380000-23266000 "+chrom+"2:89131589-90274600"
elif genomeVersion == "hg38": coordsToSubset = chrom+"14:105583730-106879900 "+chrom+"22:22025700-22924000 "+chrom+"2:88832000-90235600"

## check normal BAM or reference is available
if originalBamN is None and refGenome is None:
	sys.exit("IgCaller: error message... Normal BAM file and/or reference genome must be supplied using -N and -R, resepctively.")

## Output folder and main output file:
if outputPath is None:  wkDir = originalBamT.split("/")[-1].replace(".bam", "_IgCaller")
else: wkDir = outputPath+"/"+originalBamT.split("/")[-1].replace(".bam", "_IgCaller")
if not os.path.exists(wkDir): os.mkdir(wkDir)
SUMM = open(wkDir+"/"+originalBamT.split("/")[-1].replace(".bam", "_output_filtered.tsv"), "w")
SUMM.write("Analysis\tAnnotation\tMechanism\tScore\tMQ\tV-PhaseInfo\tV-gene_pct_identity_to_germline\tLength_considered_identity\tFunctionality\tJunction_amino_acid_sequence\tSequence\n")

## Create mini IG-bam for tumor and normal (if available):
print("IgCaller: creating IG mini BAM files...")
bamT = wkDir+"/"+originalBamT.split("/")[-1]
comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -F 1024 -q "+mapq+" -b "+originalBamT+" "+coordsToSubset+" > "+bamT  # -F 1024 (read is PCR or optical duplicate)
subprocess.call(comms, shell=True)
comms = pathToSamtools+"samtools index -@ "+threadsForSamtools+" "+bamT 
subprocess.call(comms, shell=True)
bamN = None
miniBamN = None
if originalBamN is not None:
	bamN = wkDir+"/"+originalBamN.split("/")[-1]
	comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -F 1024  -q "+mapq+" -b "+originalBamN+" "+coordsToSubset+" > "+bamN
	subprocess.call(comms, shell=True)
	comms = pathToSamtools+"samtools index -@ "+threadsForSamtools+" "+bamN 
	subprocess.call(comms, shell=True)


# 1) Iterate over each IG locus:
for GENE in ["IGH", "IGK", "IGL", "CSR"]:
	
	print("IgCaller: %s..." %GENE)
	
	chromGene, coordsToSubsetLocus, bedFile, Dseqs = getGeneralInfo(GENE, chrom, genomeVersion, inputsFolder, chrAnnot)
	
	# 2) Create mini IG-locus-specific-BAM:
	miniBamT = bamT.replace(".bam", "_miniBam.bam")
	comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -b "+bamT+" "+coordsToSubsetLocus+" > "+miniBamT 
	subprocess.call(comms, shell=True)
	comms = pathToSamtools+"samtools index -@ "+threadsForSamtools+" "+miniBamT 
	subprocess.call(comms, shell=True)
	
	miniSamT = bamT.replace(".bam", "_miniSam.sam")
	comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -F 2304 "+bamT+" "+coordsToSubsetLocus+" > "+miniSamT # -F 2304 = not priamry alignment and supplementary alignment [split read/mate pair support considered from primary alignment when analyzing IG gene rearrangements]
	subprocess.call(comms, shell=True)
	
	if bamN is not None:
		miniBamN = wkDir+"/"+bamN.split("/")[-1].replace(".bam", "_miniBam.bam")
		comms = pathToSamtools+"samtools view -@ "+threadsForSamtools+" -h -b "+bamN+" -L "+bedFile+" > "+ miniBamN
		subprocess.call(comms, shell=True)
		comms = pathToSamtools+"samtools index -@ "+threadsForSamtools+" "+miniBamN 
		subprocess.call(comms, shell=True)

	# 3) Convert SAM file to annotated table and write to disc:
	## Take samfile, get columns of interest and anotate read with large insert size (insertSize) and split/soft clipped (split) reads 
	annot_table = miniSamT.replace(".sam", "_splitinsert.tsv")
	store = convertSamToAnnotatedTable(miniSamT, chromGene, GENE)
	ANNOT_TABLE = open(annot_table, "w")
	for read, data in store.items(): ANNOT_TABLE.write("%s\n" %"\t".join([str(x) for x in data]))
	ANNOT_TABLE.close()
	
	# 4) Find the J and V genes corresponding to each split/insert size position:
	JV_list = findJandVgenes(annot_table, bedFile, GENE)
	annot_table_JV = miniSamT.replace(".sam", "_splitinsert_VJ.tsv")
	ANNOT_TABLE_JV = open(annot_table_JV, "w")
	for i in JV_list: ANNOT_TABLE_JV.write(i)
	ANNOT_TABLE_JV.close()
	
	# 5) Find combinations of J-V:
	l = findCombinationsJandV(annot_table_JV, GENE)
	
	# 6) Assign positions/breaks to each J and V pairs:
	VJ_positions, data, pos = assignPositionsToJandV(l, annot_table_JV)
	
	if GENE != "CSR":
		# 7) Append to list V,J positions and number of occurrences:
		information = addPositionsAndOccurrences(pos, bedFile, shortV, data, seq)
		
		# 8) Get J and V sequences:
		information = getJandVsequences(JV_list, information, GENE, refGenome, baseq, chromGene, bamN, miniBamT, miniBamN, depth, altDepth, tumorPurity, vafCutoff, vafCutoffNormal, pathToSamtools)
		
		# 9) Get D sequences (IGH = N-D-N, IGK/IGL = N):
		information = getDsequence(information, annot_table_JV, GENE, Dseqs)
		
		# 10) Check homology and functionality (productive/unproductive):
		information = checkHomologyAndFunctionality(information, GENE)
		
		# 11) Add mapping quality and calculate score:
		information = addMapQualAndScore(information, GENE, annot_table_JV, tumorPurity)
		
		# 12) Pre-defined filter:
		trip = predefinedFilter(information, GENE, seq, scoreCutoff)	
		
		# 13) Save output:
		Vseq = open(miniSamT.replace("_miniSam.sam", "_output_"+GENE+".tsv"), "w")
		Vseq.write("Genes\tMechanisms\tN_split\tN_insertSize\tStart_J\tEnd_J\tN_split_rescued_J\tStart_V\tEnd_V\tN_split_rescued_V\tSeq_J\tSeq_D\tSeq_V\tSeq_V_normal\tSeq\tV-gene_pct_identity_to_germline\tLength_considered_identity\tFunctionality\tJunction_amino_acid_sequence\tScore\tMQ\tV-PhaseInfo\n")
		Vseq.write("\n".join(['\t'.join(map(str, (item[i] for i in list(range(0,14))+list(range(15,22))+[14]))) for item in information]))
		Vseq.close()

		if len(trip) != 0:
			SUMM.write("\n".join([GENE+"\t"+item+"\t"+'\t'.join(map(str, map(trip[item].__getitem__, [0,19,20,13,15,16,17,18,14]))) for item in trip])+"\n") ###### modified
		else:
			SUMM.write(GENE+"\tNo rearrangement found"+"\tNA"*9+"\n")
	
	elif GENE == "CSR":
		# 14) Study coverage around CSR and return info
		class_switch, class_switch_filt, reductionMeans = classSwitchAnalysis(data, bedFile, baseq, chromGene, bamT, bamN, pathToSamtools, tumorPurity)
		
		Vseq = open(miniSamT.replace("_miniSam.sam", "_output_"+GENE+".tsv"), "w")
		Vseq.write("Genes\tClass\tScore\tAdjusted_mean_pre_break\tAdjusted_mean_post_break\tPvalue\tPct_reduction_adjusted_means\n")
		if len(class_switch) > 0:
			Vseq.write("\n".join(['\t'.join(map(str, item)) for item in class_switch])+"\n")					
		Vseq.close()
		
		if len(class_switch_filt) == 0:
			SUMM.write("CSR\tIGHM\tNo CSR found"+"\tNA"*8+"\n")
		else:
			class_switch_filt = class_switch_filt[reductionMeans.index(max(reductionMeans))] # keep only the one with highest reduction
			#class_switch_filt = class_switch_filt[pvals.index(min(pvals))]
			SUMM.write("\t".join(["CSR", class_switch_filt[0], class_switch_filt[1], str(class_switch_filt[2])]+["NA"]*7)+"\n")

# 15) Genome-wide IG translocations
print("IgCaller: genome-wide IG rearrangements...")
translocationsALL, translocationsPASS = getIgTranslocations(genomeVersion, inputsFolder, pathToSamtools, threadsForSamtools, bamT, bamN, chrom, coordsToSubset, tumorPurity, mntonco, mntoncoPass, mnnonco, mapqOnco, mncPoN)

Vseq = open(miniSamT.replace("_miniSam.sam", "_output_oncogenic_IG_rearrangements.tsv"), "w")
Vseq.write("\n".join(translocationsALL))		
Vseq.close()

if len(translocationsPASS) > 0: SUMM.write("\n".join(translocationsPASS))
else: SUMM.write("Oncogenic IG rearrangement\tNo rearrangements found"+"\tNA"*9+"\n")


# 16) Clean intermediate files and close
comms = "rm "+wkDir+"/*.sam "+wkDir+"/*miniBam.bam "+wkDir+"/*miniBam.bam.bai "+wkDir+"/*miniBam_readsSpanningJV.bam "+wkDir+"/*miniBam_readsSpanningJV.bam.bai "+wkDir+"/*MutPhased.bam "+wkDir+"/*MutPhased.bam.bai "+wkDir+"/*_output_mpileup.tsv "+wkDir+"/*_readsSpanningJV.tsv "+wkDir+"/*_miniSam_splitinsert.tsv "+wkDir+"/*_miniSam_splitinsert_VJ.tsv "+wkDir+"/*MutPhased.tsv"  
subprocess.call(comms, shell=True)
if keepMiniIgBams != "yes":
	comms = "rm "+wkDir+"/*.bam "+wkDir+"/*.bam.bai"  
	subprocess.call(comms, shell=True)

SUMM.close()
print("IgCaller: done!")
