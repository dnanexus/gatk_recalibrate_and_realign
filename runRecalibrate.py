import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator
import time


import dxpy
import subprocess, logging
import os, sys, re, math, operator, time
from multiprocessing import Pool, cpu_count

def main():    
    ## RUN DEDUPLICATE
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'

    #mappingsTable = dxpy.open_dxgtable(job['input']['mappings'][0]['$dnanexus_link'])
    if 'output_name' in job['input']:
        outputName = job['input']['output_name']
    else:
        outputName = ''
    recalibratedTable = createNewTable(job['input']['mappings'], outputName)
    #print "Mappings Table: " + mappingsTable.get_id()
    print "Recalibrated Table: " + recalibratedTable.get_id()
    
    mappingsTable = dxpy.DXGTable(job['input']['mappings'][0]['$dnanexus_link'])
    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise dxpy.AppError("The original reference genome must be attached as a detail")
    
    
    reads = 0   
    for x in job['input']['mappings']:
        table = dxpy.DXGTable(x)
        reads += int(table.describe()['length'])
    chunks = int(reads/job['input']['reads_per_job'])+2
    #chunks = 5

    #Split the genome into chunks to parallelize
    commandList = splitGenomeLengthChromosome(originalContigSet, chunks)

    excludeInterchromosome = (chunks > 1)
    markDuplicatesJobs = []
    
    #This is a Picard Mark Duplicates job run only on interchromosomal mappings in the case that the genome is split into regions
    #This is necessary because Mark Duplicates has to look at both mates in a read pair, so interchromosomal mappings must go together
    if chunks > 1:
        mapBestPracticesInput = {
            'mappings_tables': job['input']['mappings'],
            'interchromosome_bam': '',
            'interval': commandList,
            'job_number' : -1,
            'exclude_interchromosomal': False,
            'include_interchromosomal': True
        }
        interchromosomeJobId = dxpy.new_dxjob(fn_input=mapBestPracticesInput, fn_name="mapBestPractices").get_id()
        interchromosomeJobField = { 'job': interchromosomeJobId, 'field': 'bam'}
    else:
        interchromosomeJobField = ''
    
    #This runs the Picard Mark Duplicates program to deduplicate the reads
    reduceInput = {}
    resultingFiles = []
    for i in range(len(commandList)):
        print commandList[i]
        mapBestPracticesInput = {
            'mappings_tables': job['input']['mappings'],
            'recalibrated_table_id': recalibratedTable.get_id(),
            'interchromosome_bam': interchromosomeJobField,
            'interval': commandList[i],
            'job_number' : i,
            'exclude_interchromosomal': excludeInterchromosome,
            'include_interchromosomal': False,
            'reference': job['input']['reference']['$dnanexus_link'],
            'dbsnp': job['input']['dbsnp'],
            'known_indels': job['input']['known_indels'],
            'parent_input': job['input']
        }
        mapJobId = dxpy.new_dxjob(fn_input=mapBestPracticesInput, fn_name="mapBestPractices").get_id()
        reduceInput["mapJob" + str(i)] = {'job': mapJobId, 'field': 'ok'}
    reduceInput["recalibrated_table"] = recalibratedTable.get_id()

    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceBestPractices").get_id()
    job['output'] = {'recalibrated_mappings': {'job': reduceJobId, 'field': 'recalibrated_table'}}
    #job['output'] = {'recalibrated_mappings': dxpy.dxlink(recalibratedTable.get_id())}


def reduceBestPractices():
    recalibratedTable = dxpy.DXGTable(job['input']['recalibrated_table'])
    recalibratedTable.close(block=True)
    job['output']['recalibrated_table'] = dxpy.dxlink(recalibratedTable.get_id())

def writeUnmappedReads(mappingsTable, dedupTable):
    colNames = dedupTable.get_col_names()
    col = {}
    for i in range(len(colNames)):
        col[colNames[i]] = i
    for row in mappingsTable.iterate_rows():
        entry = []
        if row[col["chr"]] != '':
            break
        if row[col["chr"]] == row[col["chr2"]]:
            dedupTable.add_rows([entry])


def mapBestPractices():
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'
    
    #mappingsTable = dxpy.open_dxgtable(job['input']['mappings_table_id'])
    
    jobNumber = job['input']['job_number']
    if job['input']['interchromosome_bam'] != '':
        interchromosomeBamList = job['input']['interchromosome_bam']
        interchromosomeBam = interchromosomeBamList[jobNumber]
        print "Interchromosome BAM: " + interchromosomeBam
        dxpy.download_dxfile(interchromosomeBam, "interchromosomeBam.bam")
    
    if job['input']['exclude_interchromosomal'] and job['input']['include_interchromosomal'] == False:        
        regionFile = open("regions.txt", 'w')
        print job['input']['interval']
        print "Include interchromosomal" + str(job['input']['include_interchromosomal'])
        print "Exclude interchromosomal" + str(job['input']['exclude_interchromosomal'])
        regionFile.write(job['input']['interval'])
        regionFile.close()

    print "Converting Table to SAM"
    for i in range(len(job['input']['mappings_tables'])):
        mappingsTable = dxpy.DXGTable(job['input']['mappings_tables'][i]['$dnanexus_link']).get_id()
        if job['input']['exclude_interchromosomal'] == False and job['input']['include_interchromosomal'] == False:
            #Commented to wait for changes to dx-toolkit to percolate through
            #command = "dx-mappings-to-sam %s --output input.sam --id_as_name" % (job['input']['mappings_table_id'])
            command = "dx_mappings_to_sam2 %s --output input.sam --id_as_name --write_row_id --read_group_platform illumina --assign_read_group %d --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)
            subprocess.check_call("java -Xmx4g net.sf.picard.sam.SortSam I=input.sam O=input.sorted.sam SORT_ORDER=coordinate", shell=True)
            subprocess.check_call("mv input.sorted.sam input."+str(i)+".sam", shell=True)
            subprocess.check_call("rm input.sam", shell=True)
        elif job['input']['exclude_interchromosomal']:
            command = "dx_mappings_to_sam2 %s --output input.sam --region_index_offset -1 --id_as_name --region_file regions.txt --no_interchromosomal_mate --write_row_id --read_group_platform illumina --assign_read_group %d --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)
            subprocess.check_call("java -Xmx4g net.sf.picard.sam.SortSam I=input.sam O=input.sorted.sam SORT_ORDER=coordinate", shell=True)
            subprocess.check_call("mv input.sorted.sam input."+str(i)+".sam", shell=True)
            subprocess.check_call("rm input.sam", shell=True)
        elif job['input']['include_interchromosomal']:
            command = "dx_mappings_to_sam2 %s --output input.sam --id_as_name --only_interchromosomal_mate --write_row_id --read_group_platform illumina --assign_read_group %d --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)
            subprocess.check_call("java -Xmx4g net.sf.picard.sam.SortSam I=input.sam O=input.sorted.sam SORT_ORDER=coordinate", shell=True)
            subprocess.check_call("mv input.sorted.sam input."+str(i)+".sam", shell=True)
            subprocess.check_call("rm input.sam", shell=True)
        #if i == 0:
        #    subprocess.check_call("mv input.0.sam input.sam", shell=True)
        #else:
        #    command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles SORT_ORDER=coordinate USE_THREADING=true INPUT=dedup.rg.bam INPUT=interchromosomeBam.bam OUTPUT=input.bam"
        #    subprocess.check_call("java -Xmx4g net.sf.picard.sam.MergeSamFiles SORT_ORDER=coordinate USE_THREADING=true INPUT=dedup.rg.bam INPUT=interchromosomeBam.bam OUTPUT=input.bam", shell=True)


    readsPresent = False
    #if len(job['input']['mappings_tables']) == 1:
    #    if checkSamContainsRead("input.0.sam"):
    #        readsPresent = True
    #        subprocess.check_call("mv input.0.sam input.sam", shell=True)
    #        subprocess.check_call("wc -l input.sam", shell=True)
    #else:
    command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=input.sam USE_THREADING=true SORT_ORDER=coordinate"
    for i in range(len(job['input']['mappings_tables'])):
        if checkSamContainsRead("input."+str(i)+".sam"):
            command += " INPUT=input."+str(i)+".sam"
            readsPresent = True
    

    if readsPresent:
        subprocess.check_call(command, shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sam O=dedup.sam METRICS_FILE=metrics.txt ASSUME_SORTED=true", shell=True)
        subprocess.check_call("samtools view -bS dedup.sam > dedup.bam", shell=True)
        
        if job['input']['job_number'] == -1:
            subprocess.check_call("samtools index dedup.bam", shell=True)
            outputFiles = []
            for i in range(len(job['input']['interval'])):
                command = "samtools view -b dedup.bam " + job['input']['interval'][i].replace(" -L ", " ") + " > chromosome.bam"
                print command
                subprocess.check_call(command, shell=True)
                #subprocess.check_call("java net.sf.picard.sam.AddOrReplaceReadGroups INPUT=chromosome.bam OUTPUT=chromosome.rg.bam RGID="+str(i)+" RGLB=library RGPL=illumina RGPU="+str(i)+" SM="+str(i), shell=True)
                outputFiles.append(dxpy.upload_local_file("chromosome.bam").get_id())
            job['output']['bam'] = outputFiles
            return
        #else:
        #    subprocess.check_call("java net.sf.picard.sam.AddOrReplaceReadGroups INPUT=dedup.bam OUTPUT=dedup.rg.bam RGID="+str(job['input']['job_number'])+" RGLB=library RGPL=illumina RGPU="+str(job['input']['job_number'])+" SM="+str(job['input']['job_number']), shell=True)
    else:
        job['output']['ok'] = True
        return
    
    
    ##THIS SEGMENT PROCEEDS AFTER MARKDUPLICATES COMPLETES
    if job['input']['interchromosome_bam'] != '':
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.MergeSamFiles SORT_ORDER=coordinate USE_THREADING=true INPUT=dedup.bam INPUT=interchromosomeBam.bam OUTPUT=input.bam", shell=True)
    else:
        subprocess.check_call("mv dedup.bam input.bam", shell=True)
    subprocess.check_call("samtools index input.bam", shell=True)

    #Download the Reference Genome
    print "Converting Contigset to Fasta"
    subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['reference']), shell=True)
    
    #RealignerTargetCreator
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T RealignerTargetCreator -R ref.fa -I input.bam -o indels.intervals "
    command += job['input']['interval']
    knownIndels = ''
    
    #Download the known indels
    knownCommand = ''
    for i in range(len(job['input']['known_indels'])):
        dxpy.download_dxfile(job['input']['known_indels'][i], "indels"+str(i)+".vcf.gz")
        knownFileName = "indels"+str(i)+".vcf.gz"
        try:
            p = subprocess.Popen("tabix -f -p vcf " + knownFileName, stderr=subprocess.PIPE, shell=True)
            if '[tabix] was bgzip' in p.communicate()[1]:
                subprocess.check_call("gzip -d " + knownFileName, shell=True)
                knownFileName = "indels"+str(i)+".vcf.gz"
        except subprocess.CalledProcessError:
            raise dxpy.AppError("An error occurred decompressing dbsnp. Expected dbsnp as a gziped file.")
        knownCommand += " -known " + knownFileName        
    command += knownIndels
    
    #Find chromosomes
    regionChromosomes = []
    for x in re.findall("-L ([^:]*):\d+-\d+", job['input']['interval']):
        regionChromosomes.append(x[0])
    
    #Commented because exporting the known indels from variants takes far too long
    #knownCommand = ''
    #for i in range(len(job['input']['known_indels'])):
    #    variantsId = job['input']['known_indels'][i]['$dnanexus_link']
    #    command = "dx_variantsToVcf2 --table_id %s --output indels%d.vcf" % (variantsId, i)
    #    for x in regionChromosomes:
    #        command += " --chr " + x
    #    subprocess.check_call(command, shell=True)
    #    knownCommand += " -known indels"+str(i)+".vcf"   
    #command += knownIndels

    #Add options for RealignerTargetCreator
    if job['input']['parent_input']['window_size'] != 10:
        command += "--windowSize " + str(job['input']['parent_input']['window_size'])
    if job['input']['parent_input']['max_interval_size'] != 500:
        command += " --maxIntervalSize "  + str(job['input']['parent_input']['max_interval_size'])
    if job['input']['parent_input']['min_reads_locus'] != 4:
        command += " --minReadsAtLocus " + str(job['input']['parent_input']['min_reads_locus'])
    if job['input']['parent_input']['mismatch_fraction'] != 0.0:
        command += " --mismatchFraction " + str(job['input']['parent_input']['mismatch_fraction'])

    print command
    subprocess.check_call(command, shell=True)
    
    #Run the IndelRealigner
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T IndelRealigner -R ref.fa -I input.bam -targetIntervals indels.intervals -o realigned.bam"
    command += job['input']['interval']
    command += knownCommand
    if "consensus_model" in job['input']['parent_input']:
        if job['input']['parent_input']['consensus_model'] != "":
            if job['input']['parent_input']['consensus_model'] == "USE_READS" or job['input']['parent_input']['consensus_model'] == "KNOWNS_ONLY" or job['input']['parent_input']['consensus_model'] == "USE_SW":
                command += " --consensusDeterminationModel " + job['input']['parent_input']['consensus_model']
            else:
                raise dxpy.AppError("The option \"Consensus Determination Model\" must be either blank or one of [\"USE_READS\", \"KNWONS_ONLY\", or \"USE_SW\"], found " + job['input']['parent_input']['consensus_model'] + " instead.")
    if job['input']['parent_input']['lod_threshold'] != 5.0:
        command += " --LODThresholdForCleaning " + str(job['input']['parent_input']['lod_threshold'])
    if job['input']['parent_input']['entropy_threshold'] != 0.15:
        command += " --entropyThreshold " + str(job['input']['parent_input']['entropy_threshold'])
    if job['input']['parent_input']['max_consensuses'] != 30:
        command += " --maxConsensuses " + str(job['input']['parent_input']['max_consensuses'])
    if job['input']['parent_input']['max_insert_size_movement'] != 3000:
        command += " --maxIsizeForMovement " + str(job['input']['parent_input']['max_insert_size_movement'])
    if job['input']['parent_input']['max_position_move'] != 200:
        command += " --maxPositionalMoveAllowed " + str(job['input']['parent_input']['max_position_move'])
    if job['input']['parent_input']['max_reads_consensus'] != 120:
        command += " --maxReadsForConsensus " + str(job['input']['parent_input']['maxReadsForRealignment'])
    if job['input']['parent_input']['max_reads_realignment'] != 20000:
        command += " --maxReadsForRealignment " + str(job['input']['parent_input']['max_reads_realignment'])
    
    print command
    subprocess.check_call(command, shell=True)
    
    #Download dbsnp
    dxpy.download_dxfile(job['input']['dbsnp'], "dbsnp.vcf.gz")
    
    dbsnpFileName = 'dbsnp.vcf.gz'
    try:
        p = subprocess.Popen("tabix -f -p vcf dbsnp.vcf.gz", stderr=subprocess.PIPE, shell=True)
        if '[tabix] was bgzip' in p.communicate()[1]:
            subprocess.check_call("gzip -d dbsnp.vcf.gz", shell=True)
            dbsnpFileName = 'dbsnp.vcf'
    except subprocess.CalledProcessError:
        raise dxpy.AppError("An error occurred decompressing dbsnp. Expected dbsnp as a gziped file.")
    
    #subprocess.check_call("gzip -d dbsnp.vcf.gz", shell=True)

    #Count Covariates
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T CountCovariates -R ref.fa -recalFile recalibration.csv -I realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --standard_covs"
    command += "-knownSites " + dbsnpFileName
    command += job['input']['interval']
    command += " --num_threads " + str(cpu_count())
    
    if "solid_recalibration_mode" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_recalibration_mode'] != "":
            if job['input']['parent_input']['solid_recalibration_mode'] == "THROW_EXCEPTION" or job['input']['parent_input']['solid_recalibration_mode'] == "LEAVE_READ_UNRECALIBRATED" or job['input']['parent_input']['solid_recalibration_mode'] == "PURGE_READ":
                command += " --solid_recal_mode " + job['input']['parent_input']['solid_recalibration_mode']
            else:
                raise dxpy.AppError("The option \"SOLiD Recalibration Mode\" must be either blank or one of [\"THROW_EXCEPTION\", \"LEAVE_READ_UNRECALIBRATED\", or \"PURGE_READ\"], found " + job['input']['parent_input']['solid_recalibration_mode'] + " instead.")
    if "solid_nocall_strategy" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_nocall_strategy'] != "":
            if job['input']['parent_input']['solid_nocall_strategy'] == "THROW_EXCEPTION" or job['input']['parent_input']['solid_nocall_strategy'] == "LEAVE_READ_UNRECALIBRATED" or job['input']['parent_input']['solid_nocall_strategy'] == "PURGE_READ":
                command += " --solid_nocall_strategy " + job['input']['parent_input']['solid_nocall_strategy']
            else:
                raise dxpy.AppError("The option \"SOLiD No-call Strategy\" must be either blank or one of [\"THROW_EXCEPTION\", \"LEAVE_READ_UNRECALIBRATED\", or \"PURGE_READ\"], found " + job['input']['parent_input']['solid_nocall_strategy'] + " instead.")
    if 'context_size' in job['input']['parent_input']:
        command += " --context_size " + str(job['input']['parent_input']['context_size'])
    if 'nback' in job['input']['parent_input']:
        command += " --homopolymer_nback " + str(job['input']['parent_input']['nback'])
    if job['input']['parent_input']['cycle_covariate']:
        command += " -cov CycleCovariate"
    if job['input']['parent_input']['dinuc_covariate']:
        command += " -cov DinucCovariate"
    if job['input']['parent_input']['primer_round_covariate']:
        command += " -cov PrimerRoundCovariate"
    if job['input']['parent_input']['mapping_quality_covariate']:
        command += " -cov MappingQualityCovariate"
    if job['input']['parent_input']['homopolymer_covariate']:
        command += " -cov HomopolymerCovariate"
    if job['input']['parent_input']['gc_content_covariate']:
        command += " -cov GCContentCovariate"
    if job['input']['parent_input']['position_covariate']:
        command += " -cov PositionCovariate"
    if job['input']['parent_input']['minimum_nqs_covariate']:
        command += " -cov MinimumNQSCovariate"
    if job['input']['parent_input']['context_covariate']:
        command += " -cov ContextCovariate"

    print command
    subprocess.check_call(command, shell=True)
    
    #Table Recalibration
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T TableRecalibration -R ref.fa -recalFile recalibration.csv -I realigned.bam -o recalibrated.bam --doNotWriteOriginalQuals"
    command += job['input']['interval']
    if "solid_recalibration_mode" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_recalibration_mode'] != "":
            command += " --solid_recal_mode " + job['input']['parent_input']['solid_recalibration_mode']
    if "solid_nocall_strategy" in job['input']['parent_input']:
        if job['input']['parent_input']['solid_nocall_strategy'] != "":
            command += " --solid_nocall_strategy " + job['input']['parent_input']['solid_nocall_strategy']
    if 'context_size' in job['input']['parent_input']:
        command += " --context_size " + str(job['input']['parent_input']['context_size'])
    if 'nback' in job['input']['parent_input']:
        command += " --homopolymer_nback " + str(job['input']['parent_input']['nback'])
    if job['input']['parent_input']['preserve_qscore'] != 5:
        command += " --preserve_qscores_less_than " + str(job['input']['parent_input']['preserve_qscore'])
    if 'smoothing' in job['input']['parent_input']:
        command += " --smoothing " + str(job['input']['parent_input']['smoothing'])
    if 'max_quality' in job['input']['parent_input']:
        command += " --max_quality_score " + str(job['input']['parent_input']['max_quality'])
    print command
    subprocess.check_call(command, shell=True)
    
    subprocess.check_call("samtools view -h -o recalibrated.sam recalibrated.bam", shell=True)
    
    result = dxpy.upload_local_file("recalibrated.bam")
    job['output']['recalibrated_bam'] = dxpy.dxlink(result.get_id())

    

    #Read SAM, extracting chr, lo, hi, qual, cigar, duplicate flag, chr2, lo2, hi2


    recalibratedTable = dxpy.DXGTable(job['input']['recalibrated_table_id'])

    default = {}

    recalibratedColNames = recalibratedTable.get_col_names()
    recalibratedCol = {}
    for i in range(len(recalibratedColNames)):
        recalibratedCol[recalibratedColNames[i]] = i
    if recalibratedCol.get("chr2") != None:
        mateLocations = {}

    for x in recalibratedTable.describe()['columns']:
        if "int" in x["type"]:
            default[x["name"]] = sys.maxint
        elif x["type"] == "float":
            default[x["name"]] = float(sys.maxint)
        elif x["type"] == "boolean":
            default[x["name"]] = False
        else:
            default[x["name"]] = ""

    print "Recalibrated SAM: " + dxpy.upload_local_file("recalibrated.sam").get_id()

    result = {}
    for line in open("recalibrated.sam", 'r'):
        if line[0] != "@":
            rowId = int(re.findall("ZD:Z:(\d+)", line)[0])
            tabSplit = line.split("\t")
            lo = int(tabSplit[3])-1
            qual = tabSplit[10]
            alignLength = 0
            cigar = re.split('(\d+)', tabSplit[5])
            duplicate = (int(tabSplit[1]) & 0x400 == True)
            if int(tabSplit[1]) & 0x10:
                qual = qual[::-1]
            for p in range(len(cigar)):
                c = cigar[p]
                if c == 'M' or c == 'D' or c == 'N' or c == 'X' or c == 'P' or c == '=':
                    alignLength += int(cigar[p-1])
            cigar = tabSplit[5]
            hi = lo + alignLength
            if recalibratedCol.get("chr2") != None:
                templateId = tabSplit[0]
                if int(tabSplit[1]) & 0x1 & 0x40:
                    mateLocations[templateId] = {0: {"lo":lo, "hi":hi}}
                elif int(tabSplit[1]) & 0x1 & 0x80:
                    mateLocations[templateId] = {1: {"lo":lo, "hi":hi}}
            result[rowId]= {"lo":lo, "hi":hi, "cigar":cigar, "qual":qual, "duplicate":duplicate}
            

    #Lookup start row for each chromosome
    for i in range(len(job['input']['mappings_tables'])):
        mappingsTable = dxpy.DXGTable(job['input']['mappings_tables'][i]['$dnanexus_link'])
        mappingsCol = {}
        mappingsColNames = mappingsTable.get_col_names()
        for i in range(len(mappingsColNames)):
            mappingsCol[mappingsColNames[i]] = i
        #Iterate through the table, if id encountered add to recalibrated table
        for x in regionChromosomes:
            query = mappingsTable.genomic_range_query(x[0],0,sys.maxint, index='gri')
            for row in mappingsTable.get_rows(query=query, limit=1)['data']:
                startRow =  row[0]
                for row in mappingsTable.iterate_rows(start=startRow):
                    if row[5] != x[0]:
                        break
                    rowId = row[0]
                    if result.get(rowId) != None:
                        entry = []
                        rowEntries = row[1:]
                        for j in range(len(recalibratedColNames)):
                            c = recalibratedColNames[j]
                            if mappingsCol.get(c) != None:
                                entry.append(rowEntries[mappingsCol[c]])
                            else:
                                entry.append(default[c])
                        entry[recalibratedCol["lo"]] = result[rowId]["lo"]
                        entry[recalibratedCol["hi"]] = result[rowId]["hi"]
                        entry[recalibratedCol["quality"]] = result[rowId]["qual"]
                        entry[recalibratedCol["duplicate"]] = result[rowId]["duplicate"]
                        entry[recalibratedCol["cigar"]] = result[rowId]["cigar"]
                        if recalibratedCol.get("chr2") != None:
                            location = mateLocations.get(str(i)+"_"+str(entry[recalibratedCol["template_id"]]))
                            if location != None:
                                mate = location.get(int(entry[recalibratedCol["mate_id"]]))
                                if mate != None:
                                    entry[recalibratedCol["lo2"]] = mate["lo"]
                                    entry[recalibratedCol["hi2"]] = mate["hi"]
                        recalibratedTable.add_rows([entry])
        recalibratedTable.flush()

    job['output']['ok'] = True


def splitGenomeLengthChromosome(contig_set, chunks):
    details = dxpy.DXRecord(contig_set).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']
    
    commandList = []
    chunkSizes = []
    
    totalSize = sum(sizes)
    readsPerChunk = int(float(totalSize)/chunks)
    
    
    for i in range(chunks):
        commandList.append('')
        chunkSizes.append(0)

    chromosome = 0
    position = 0
    while chromosome < len(names):
        chunkSizes[position] += sizes[chromosome]
        commandList[position] += " -L %s:%d-%d" % (names[chromosome], 1, sizes[chromosome])
        chromosome += 1
        if chunkSizes[position] > readsPerChunk and position < chunks-1:
            position += 1

    return commandList


def checkIntervalRange(includeList, chromosome, lo, hi):
    included = False
    command = ''
    if len(includeList) == 0:
        return " -L %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
            min = lo
            max = hi
            if (lo >= x[0] and lo <= x[1]) or (hi <= x[1] and hi >= x[0]):
                if lo >= x[0] and lo <= x[1]:
                    min = lo
                elif lo <= x[0]:
                    min = x[0]
                if hi <= x[1] and hi >= x[0]:
                    max = hi
                elif hi >= x[1]:
                    max = x[1]
                command += " -L %s:%d-%d" % (chromosome, min, max)
    return command

def checkSamContainsRead(samFileName):
    for line in open(samFileName, 'r'):
        if line[0] != "@":
            return True
    return False


def createNewTable(mappingsArray, recalibratedName):
    
    columns = []
    tags = []
    indices = []
    types = []
    for i in range(len(mappingsArray)):
        oldTable = dxpy.DXGTable(mappingsArray[i]['$dnanexus_link'])        
        for x in oldTable.describe()['columns']:
            if x not in columns:
                columns.append(x)
        for x in oldTable.describe()['indices']:
            if x not in indices:
                indices.append(x)
        for x in oldTable.describe()['tags']:
            if x not in tags:
                tags.append(x)
        for x in oldTable.describe()['types']:
            if x not in types:
                types.append(x)
    oldTable = dxpy.DXGTable(mappingsArray[0]['$dnanexus_link'])
    if recalibratedName == '':
        recalibratedName = oldTable.describe()['name'] + " Realigned and Recalibrated"        
    details = oldTable.get_details()    
    newTable = dxpy.new_dxgtable(columns=columns, indices=indices)
    newTable.add_tags(tags)
    newTable.set_details(details)
    newTable.add_types(types)
    newTable.rename(recalibratedName)
    
    return newTable

