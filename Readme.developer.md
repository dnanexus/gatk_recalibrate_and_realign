GATK Indel Realigner and Quality Recalibration, Advanced Readme
===============================================================

Introduction
------------

The GATK Best Practices pipeline runs a set of programs which realign and recalibrate the quality of previously mapped reads. These steps are recommended by the Broad Institute to get the best calls out of a set of mapped reads. The app consists of three conceptual steps: marking duplicates, realigning around indels, and recalibrating quality scores. This app creates a new table of the recalibrated mappings. Unmapped reads will not be included in the new table.

### Resource Files

Certain steps of the best practices pipeline use files containing known snps and indels. GATK recommends the use of dbsnp and two known indels files. For the b37 reference genome, these can be found in the public project **Reference Genomes** in the **b37** folder.

The dbsnp file is named **dbsnp_135.b37.vcf.gz**

The two indels files are named **Mills_and_1000G_gold_standard.indels.b37.vcf.gz** and **1000G_phase1.indels.b37.vcf.gz**

In both cases, these are [bgzipped](http://samtools.sourceforge.net/tabix.shtml) files. In place of these files, you may use any dbsnp or indels file that would work with command line GATK. The app will run faster if the provided files are bgzipped.

### Marking Duplicates

The first step of the pipeline is to mark all apparent duplicate reads or read pairs. This is done with [Picard Tools](http://picard.sourceforge.net/) [MarkDuplicates](http://picard.sourceforge.net/command-line-overview.shtml#MarkDuplicates) program. For unpaired reads, this program finds all reads mapped to the same 5' position, identifies finds the read with the highest number of bases with a quality score of 15 or higher and marks all other reads as duplicates. For paired end reads, it uses the 5' mapping location of both reads in the pair, calculating the number of bases with a quality score of 15 or higher in both reads of the pair and marking as duplicates all other read pairs which share the same 5' mapping as both of the highest quality pairs. For more information, see the [Picard Tools FAQ](http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page#Q:_How_does_MarkDuplicates_work.3F)

With the default parameters, the Picard MarkDuplicates component runs the following command:

    java -jar MarkDuplicates.jar INPUT=input.bam OUTPUT=dedup.bam METRICS=metrics.txt ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

Parameters such as REMOVE_DUPLICATES can be changed by changing the app input parameters.

### Indel Realignment

In order to align reads to a reference in a reasonable amount of time, aligners such as BWA only look at a single read or read pair at a time. As a result, aligners can miss the context provided by other reads mapped to a similar location, which when looked at together would indicate that an indel is present and would provide information about the best way to align each individual read in the region.

In indel realignment, this information is extracted from the mappings and used to realign reads in the vicinity of an apparent indel. This is done with two [GATK](http://www.broadinstitute.org/gatk/) programs, [RealignerTargetCreator](http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_RealignerTargetCreator.html), which identifies regions which may benefit from realignment, and [IndelRealigner](http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_indels_IndelRealigner.html), which performs the realignment.

With the default parameters, the indel realignment component runs the following commands:

     java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R ref.fa -I input.bam -o indels.intervals

     java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ref.fa -I input.bam -targetIntervals indels.intervals \
     -o realigned.bam -known indels1.vcf.gz -known indels2.vcf.gz -known ...

### Quality Recalibration

The quality associated with sequencing runs is an attempt by the sequencing machine to provide a measure of the confidence that the assigned base at a position is correct. These quality scores have been tuned by the makers of the sequencing equipment to, on average, give the best values. However, a number of factors may mean that the quality scores over or under-estimate the quality for bases in particular contexts.

The quality recalibrator uses the fact that most true variants in a sequencing run have already been found in dbsnp to recalibrate the quality scores of the reads. This uses a number of a parameters which co-vary with the reported quality of the read and correct for biases associated with these parameters (for example, it may be that quality scores for a particular read group is generally lower or that a the quality of a particular dinucleotide is over reported).

This implementation uses the GATK programs CountCovariates, which calculates how various parameters correlate with quality, and TableRecalibrator, which performs the recalibration. These programs have been renamed and their options changed in GATK 2.0 and the legacy documentation is no longer available.

With the default parameters, the indel realignment component runs the following commands:

     java -jar GenomeAnalysisTK.jar -T CountCovariates -R ref.fa -recalFile recalibration.csv \
     -I realigned.bam -knownSites dbsnp.vcf.gz

     java -jar GenomeAnalysisTK.jar -T TableRecalibration -R ref.fa -recalFile recalibration.csv \
     -I realigned.bam -o recalibrated.bam --doNotWriteOriginalQuals

### General Options

<table>
<tr><th>Name</th><th>Label</th><th>Corresponding command line option/Notes</th><th>Class/Type</th><th>Default</th></tr>
<tr><td>mappings</td><td>Mappings Table(s)</td><td>-I (All of the mappings to realign and recalibrate. The app will act on all reads in all tables and the output will be combined)</td><td>array of gtables of type LetterMappings.</td><td>Mandatory</td>
<tr><td>reference</td><td>Reference Genome</td><td>-R </td><td>record of type:ContigSet</td><td>Mandatory</td></tr>
<tr><td>output_name</td><td>Output Name</td><td></td><td>string</td><td>Default is empty string and will name the new object  the same name of the first mapping object with "Recalibrated and Realigned" added to the end</td></tr>
<tr><td>dbsnp</td><td>dbSNP</td><td>-knownSites </td><td>a tar.gz file. If the file has been [bgzip'd](http://samtools.sourceforge.net/tabix.shtml), the program will take advantage of this</td><td>Mandatory</td></tr>
<tr><td>known_indels</td><td>Known Indels</td><td>--known (Resources which contain known, common indels in the genome)</td><td>an array of tar.gz files</td><td></td>Optional</tr>
</table>

### RealignerTargetCreator Options
<table>
<tr><th>Name</th><th>Label</th><th>Corresponding command line option/Notes</th><th>Class/Type</th><th>Default</th></tr>
<tr><td>max_interval_size</td><td>Max Interval Size</td><td>--maxIntervalSize (Intervals larger than this will not be used in indel realignment)</td><td>int</td><td>500</td>
<tr><td>min_reads_locus</td><td>Min Reads at Locus</td><td>--minReadsAtLocus (Mimimum reads at a locus to enable using entropy calculation for indel realignment)</td><td>int</td><td>4</td>
<tr><td>mismatch_fraction</td><td>Mismatch Fraction</td><td>--mismatchFraction (Fraction of base qualities needing to mismatch for a position to have high entropy for indel realignment target creator)</td><td>float</td><td>0.0</td>
<tr><td>window_size</td><td>Window Size</td><td>--windowSize (Window size for calculating entropy or SNP clusters)</td><td>int</td><td>10</td>
</table>

### IndelRealigner Options
<table>
<tr><th>Name</th><th>Label</th><th>Corresponding command line option/Notes</th><th>Class/Type</th><th>Default</th></tr>
<tr><td>consensus_model</td><td>Consensus Determination Model</td><td>--consensusDeterminationModel (Determines how to compute the possible alternate consenses in indel realignment. Must be one of the following values: [USE_READS, KNOWNS_ONLY, USE_SW])</td><td>string</td><td>Optional</td>
<tr><td>lod_threshold</td><td>LOD Cleaning Threshold</td><td>--LODThresholdForCleaning (LOD threshold above which the cleaner will clean in indel realignment. This is a measure of whether improvement is significant enough to merit realignment. Lower values are recommended in cases of low coverage or looking for indels with low allele frequency)</td><td>float</td><td>5.0</td>
<tr><td>entropy_threshold</td><td>Entropy Threshold</td><td>--entropyThreshold (Percentage of mismatches at a locus to be considered having high intropy in the indel realigner)</td><td>float</td><td>0.15</td>
<tr><td>max_consensuses</td><td>Max Consensuses</td><td>--maxConsensuses (Max alternate consensuses to try - higher numbers improve performance in deep coverage)</td><td>int</td><td>30</td>
<tr><td>max_insert_size_movement</td><td>Max Insert Movement Size</td><td>--maxIsizeForMovement (Maximum insert size of read pairs that realignment attempted for)</td><td>int</td><td>3000</td>
<tr><td>max_position_move</td><td>Max Position Move</td><td>--maxPositionalMoveAllowed (Maximum positional move in basepairs that a read can be adjusted during realignment)</td><td>int</td><td>200</td>
<tr><td>max_reads_consensus</td><td>Max Reads for Consensus</td><td>--maxReadsForConsensus (Maximum reads used for finding the alternate consensuses - higher numbers improve performance in deep coverage)</td><td>int</td><td>120</td>
<tr><td>max_reads_realignment</td><td>Max Reads for Realignment</td><td>--maxReadsForRealignment (Maximum reads allowed at an interval for realignment)</td><td>int</td><td>20000</td>
</table>

### CountCovariates and TableRecalibrator Options

Because several options are shared across CountCovariates and TableRecalibrator, their options are listed together. The covariates ReadGroup and ReportedQuality will always be used as covariates, regardless of the other covariates selected.

<table>
<tr><th>Name</th><th>Label</th><th>Corresponding command line option/Notes</th><th>Class/Type</th><th>Default</th></tr>
<tr><td>solid_recalibration_mode</td><td>SOLiD Recalibration Mode</td><td>--solid_recal_mode (Only applies to SOLiD sequencing. How to recalibrate bases in which the reference was inseted. If entered, must be one of the following options: [DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, REMOVE_REF_BIAS])</td><td>string</td><td>Optional</td>
<tr><td>solid_nocall_mode</td><td>SOLiD No-call Mode</td><td>--solid_nocall_strategy (Only applies to SOLiD sequencing. Defines behavoir when no-call encountered in color space. If entered, must be one of the following options: [DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, REMOVE_REF_BIAS])</td><td>string</td><td>Optional</td>
<tr><td>context_size</td><td>Count Covariate Context Size</td><td>--context_size (Size of the k-mer context used in count covariates)</td><td>int</td><td>Optional</td>
<tr><td>nback</td><td>Count Covariants N-back Size</td><td>--homopolymer_nback (The number of previous bases to look at in HomopolymerCovariate)</td><td>int</td><td>Optional</td>
<tr><td>cycle_covariate</td><td>Use Cycle Covariate</td><td>-cov CycleCovariate (Use cycle covariation in the quality recalibration process)</td><td>boolean</td><td>true</td>
<tr><td>dinuc_covariate</td><td>Use Dinucleotide Covariate</td><td>-cov DinucCovariate (Use dinucleotide covariation in the quality recalibration process)</td><td>boolean</td><td>true</td>
<tr><td>primer_round_covariate</td><td>Use Primer Round Covariate</td><td>-cov PrimerRoundCovariate (Use primer round covariation in the quality recalibration process)</td><td>boolean</td><td>false</td>
<tr><td>mapping_quality_covariate</td><td>Use Mapping Quality Covariate</td><td>-cov MappingQualityCovariate (Use mapping quality covariation in the quality recalibration process)</td><td>boolean</td><td>false</td>
<tr><td>gc_content_covariate</td><td>Use GC Content Covariate</td><td>-cov GCContentCovariate (Use GC content covariation in the quality recalibration process)</td><td>boolean</td><td>false</td>
<tr><td>position_covariate</td><td>Use Position Covariate</td><td>-cov PositionCovariate (Use position covariation in the quality recalibration process)</td><td>boolean</td><td>false</td>
<tr><td>minimum_nqs_covariate</td><td>Use Minimum NQS Covariate</td><td>-cov MinimumNQSCovariate (Use minimum NQS covariation in the quality recalibration process)</td><td>boolean</td><td>false</td>
<tr><td>context_covariate</td><td>Use Context Covariate</td><td>-cov ContextCovariate (Use context covariation in the quality recalibration process)</td><td>boolean</td><td>false</td>
<tr><td>preserve_qscore</td><td>Preserve Q-Scores Less Than </td><td>--preserve_qscores_less_than (Do not recalibrate quality scores below this threshold. Since many base callers use quality scores below 5 to indicate random or bad bases, it is often unsafe to recalibrate these bases)</td><td>int</td><td>5</td>
<tr><td>smoothing</td><td>Smoothing Counts </td><td>--smoothing (Number of imaginary counts to add to each bin in order to smooth out binds with few data points)</td><td>int</td><td>Optional</td>
<tr><td>max_quality</td><td>Maximum Quality Score </td><td>--max_quality_score (The value at which to cap the quality scores)</td><td>int</td><td>Optional</td>
