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

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    mappingSchema = mappingsTable.describe()['columns']

    #chunks = int(mappingsTable.describe()['length']/100000000)+1
    chunks = 5

    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise dxpy.AppError("The original reference genome must be attached as a detail")
        
    #Split the genome into chunks to parallelize
    commandList = splitGenomeLengthChromosome(originalContigSet, chunks)

    reduceInput = {}
    excludeInterchromosome = (chunks > 1)
    markDuplicatesJobs = []
    reduceInput = {}
    
    #This runs the Picard Mark Duplicates program to deduplicate the reads
    for i in range(len(commandList)):
        print commandList[i]
        if len(commandList[i]) > 0:
            mapMarkDuplicatesInput = {
                'mappings_table_id': mappingsTable.get_id(),
                #'dedup_table_id': dedupTable.get_id(),
                'interval': commandList[i],
                'discard_duplicates': job['input']['discard_duplicates'],
                'discard_unmapped': job['input']['discard_unmapped'],
                'job_number' : i,
                'exclude_interchromosomal': excludeInterchromosome,
                'include_interchromosomal': False
            }
            # Run a "map" job for each chunk
            markDuplicatesJobs.append(dxpy.new_dxjob(fn_input=mapMarkDuplicatesInput, fn_name="mapMarkDuplicates").get_id())
    
    #This is a Picard Mark Duplicates job run only on interchromosomal mappings in the case that the genome is split into regions
    #This is necessary because Mark Duplicates has to look at both mates in a read pair, so interchromosomal mappings must go together
    if chunks > 1:
        mapMarkDuplicatesInput = {
            'mappings_table_id': mappingsTable.get_id(),
            #'dedup_table_id': dedupTable.get_id(),
            'interval': commandList,
            'discard_duplicates': job['input']['discard_duplicates'],
            'discard_unmapped': job['input']['discard_unmapped'],
            'job_number' : i,
            'exclude_interchromosomal': False,
            'include_interchromosomal': True
        }
        interchromosomeJob = dxpy.new_dxjob(fn_input=mapMarkDuplicatesInput, fn_name="mapMarkDuplicates").get_id()
    else:
        interchromosomeJob = ''
    
    #knownIndels = ['indels/Mills_and_1000G_gold_standard.indels.b37.vcf.gz', 'indels/1000G_phase1.indels.b37.vcf.gz']
    for i in range(len(commandList)):
        mapRealignInput = {
            'region_bam' : { 'job': markDuplicatesJobs[i], 'field': 'bam'},
            'interchromosome_bam': {'job': interchromosomeJob, 'field': 'bam'},
            'region': commandList[i],
            'reference': contigSetId,
            'dbsnp': job['input']['dbsnp'],
            'known_indels': job['input']['known_indels'],
            'job_number': i,
            }
        realignJobId = dxpy.new_dxjob(fn_input=mapRealignInput, fn_name="mapRealigner").get_id()
        reduceInput["realignJob" + str(chunks)] = {'job': realignJobId, 'field': 'id'}
        #if job['input']['discard_unmapped'] == False:
            #writeUnmappedReads(mappingsTable, dedupTable)

    #reduceInput['dedupTableId'] = dedupTable.get_id()
    reduceJob = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceRealigner")
    reduceJobId = reduceJob.get_id()
    job['output'] = {'recalibrated_bam': {'job': reduceJobId, 'field': 'recalibrated_bam'}}

def mapRealigner():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar'
    fileNum = 0
    print job
    region = job['input']['region']
    regionBamList = job['input']['region_bam']
    if len(regionBamList) > 0:
        regionBam = regionBamList[0]
        print "region BAM: "+regionBam
        jobNumber = job['input']['job_number']
        print "Converting Contigset to Fasta"
        subprocess.check_call("contigset2fasta %s ref.fa" % (job['input']['reference']), shell=True)
    
        if job['input']['interchromosome_bam'] != '':
            interchromosomeBamList = job['input']['interchromosome_bam']
            if len(interchromosomeBamList) < job['input']['job_number']:
                interchromosomeBam = interchromosomeBamList[jobNumber]
                print "Interchromosome BAM: " + interchromosomeBAM
                dxpy.download_dxfile(regionBam, "regionBam.bam")
                dxpy.download_dxfile(interchromosomeBam, "interchromosomeBam.bam")
                subprocess.check_call("java -Xmx4g net.sf.picard.sam.MergeSamFiles SORT_ORDER=coordinate USE_THREADING=true INPUT=regionBam.bam INPUT=interchromosomeBam.bam OUTPUT=input.bam", shell=True)
            else:
                dxpy.download_dxfile(regionBam, "input.bam")
        else:
            dxpy.download_dxfile(regionBam, "input.bam")
            
        subprocess.check_call("samtools index input.bam", shell=True)
        
        #Run the RealignerTargetCreator
        command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T RealignerTargetCreator -R ref.fa -I input.bam -o indels.intervals "
        command += job['input']['region']
        knownIndels = ''
        
        #Download the known indels
        knownCommand = ''
        for i in range(len(job['input']['known_indels'])):
            dxpy.download_dxfile(job['input']['known_indels'][i], "indels"+str(i)+".vcf.gz")
            subprocess.check_call("gzip -d indels"+str(i)+".vcf.gz", shell=True)
            knownCommand += " -known indels"+str(i)+".vcf"        
        command += knownIndels
        print command
        subprocess.check_call(command, shell=True)
        
        #Run the Realigner
        command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T IndelRealigner -R ref.fa -I input.bam -targetIntervals indels.intervals -o realigned.bam"
        command += job['input']['region']
        command += knownCommand
        print command
        subprocess.check_call(command, shell=True)
        
        #Download dbsnp
        dxpy.download_dxfile(job['input']['dbsnp'], "dbsnp.vcf.gz")
        
        #Count Covariates
        command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T CountCovariates -R ref.fa -recalFile recalibration.csv -knownSites dbSNP.vcf.gz -I realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariante --standard_covs"
        command += " --num_threads " + str(cpu_count())
        print command
        subprocess.check_call(command, shell=True)
        
        #Table Recalibration
        command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T TableRecalibration -R ref.fa -recalFile recalibration.csv -I realigned.bam -o recalibrated.bam"
        command += " --num_threads " + str(cpu_count())
        print command
        subprocess.check_call(command, shell=True)
        
        subprocess.check_call("samtools view -h -o recalibrated.sam recalibrated.bam", shell=True)
        
        result = dxpy.upload_local_file("recalibrated.sam")
        job['output'] = {'recalibrated_bam': [dxpy.dxlink(result.get_id())]}
    else:
        job['output'] = {'recalibrated_bam': []}

def reduceRealigner():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/MergeSamFiles.jar'
    fileNum = 0
    command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=input.sorted.bam SORT_ORDER=coordinate USE_THREADING=true"
    print job
    originalInput = dxpy.DXJob(job['id']).describe()['originalInput']
    print originalInput
    for k, v in originalInput.iteritems():
        mapJob = dxpy.DXJob(v['job'])
        fileList = mapJob.describe()['output']['recalibrated_bam']
        if len(fileList) > 0:
            fileName = fileList[0]
            dxpy.download_dxfile(fileName, "sam"+str(fileNum)+".bam")
            command += " INPUT=sam"+str(fileNum)+".bam"
            fileNum += 1
    
    subprocess.check_call(command, shell=True)
    result = dxpy.upload_local_file("input.sorted.bam")
    job['output'] = {'recalibrated_bam': dxpy.dxlink(result.get_id())}
    
    #t = dxpy.open_dxgtable(job['input']['dedupTableId'])
    #print "Closing Table"
    #t.close()
    #job['output']['deduplicated mappings'] = dxpy.dxlink(t.get_id())

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


def mapMarkDuplicates():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar'
    
    mappingsTable = dxpy.open_dxgtable(job['input']['mappings_table_id'])
    
    if job['input']['exclude_interchromosomal'] and job['input']['include_interchromosomal'] == False:        
        regionFile = open("regions.txt", 'w')
        print job['input']['interval']
        print "Include interchromosomal" + str(job['input']['include_interchromosomal'])
        print "Exclude interchromosomal" + str(job['input']['exclude_interchromosomal'])
        regionFile.write(job['input']['interval'])
        regionFile.close()

    print "Converting Table to SAM"
    if job['input']['exclude_interchromosomal'] == False and job['input']['include_interchromosomal'] == False:
        command = "dx-mappings-to-sam %s --output input.sam --id_as_name" % (job['input']['mappings_table_id'])
        if job['input']['discard_unmapped']:
            command += " --discard_unmapped"
        print command
        subprocess.check_call(command, shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.SortSam I=input.sam O=input.sorted.sam SORT_ORDER=coordinate", shell=True)
        subprocess.check_call("mv input.sorted.sam input.sam", shell=True)
    elif job['input']['exclude_interchromosomal']:
        command = "dx-mappings-to-sam %s --output input.sam --region_index_offset -1 --id_as_name --region_file regions.txt --no_interchromosomal_mate" % (job['input']['mappings_table_id'])
        print command
        subprocess.check_call(command, shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.SortSam I=input.sam O=input.sorted.sam SORT_ORDER=coordinate", shell=True)
        subprocess.check_call("mv input.sorted.sam input.sam", shell=True)
    elif job['input']['include_interchromosomal']:
        command = "dx-mappings-to-sam %s --output input.sam --id_as_name --only_interchromosomal_mate" % (job['input']['mappings_table_id'])
        print command
        subprocess.check_call(command, shell=True)
        subprocess.check_call("java -Xmx4g net.sf.picard.sam.SortSam I=input.sam O=input.sorted.sam SORT_ORDER=coordinate", shell=True)
        subprocess.check_call("mv input.sorted.sam input.sam", shell=True)

    if checkSamContainsRead("input.sam"):
        subprocess.call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sam O=dedup.sam METRICS_FILE=metrics.txt ASSUME_SORTED=true", shell=True)
        subprocess.check_call("samtools view -bS dedup.sam > dedup.bam", shell=True)
        
        if job['input']['include_interchromosomal']:
            subprocess.check_call("samtools index dedup.bam", shell=True)
            outputFiles = []
            for i in range(len(job['input']['interval'])):
                command = "samtools view -b dedup.bam " + job['input']['interval'][i].replace(" -L ", " ") + " > chromosome.bam"
                print command
                subprocess.check_call(command, shell=True)
                outputFiles.append(dxpy.upload_local_file("chromosome.bam").get_id())
            job['output']['bam'] = outputFiles
        else:
            outputFile = dxpy.upload_local_file("dedup.bam")
            job['output']['bam'] = [outputFile.get_id()]
    else:
        job['output']['bam'] = []

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


def buildGatkCommand(job):
    
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T UnifiedGenotyper -R ref.fa -I input.sorted.bam -o output.vcf "
    if job['input']['output_mode'] != "EMIT_VARIANTS_ONLY":
        command += " -out_mode " + (job['input']['output_mode'])
    if job['input']['call_confidence'] != 30.0:
        command += " -stand_call_conf " +str(job['input']['call_confidence'])
    if job['input']['emit_confidence'] != 30.0:
        command += " -stand_emit_conf " +str(job['input']['emit_confidence'])
    if job['input']['pcr_error_rate'] != 0.0001:
        command += " -pcr_error " +str(job['input']['pcr_error_rate'])
    if job['input']['heterozygosity'] != 0.001:
        command += " -hets " + str(job['input']['heterozygosity'])
    if job['input']['indel_heterozygosity'] != 0.000125:
        command += " -indelHeterozygosity " + str(job['input']['indel_heterozygosity'])
    if job['input']['genotype_likelihood_model'] != "BOTH":
        command += " -glm " + job['input']['genotype_likelihood_model']
    if job['input']['minimum_base_quality'] != 17:
        command += " -mbq " + str(job['input']['minimum_base_quality'])
    if job['input']['max_alternate_alleles'] != 3:
        command += " -maxAlleles " + str(job['input']['max_alternate_alleles'])
    if job['input']['max_deletion_fraction'] != 0.05:
        command += " -deletions " + str(job['input']['max_deletion_fraction'])
    if job['input']['min_indel_count'] != 5:
        command += " -minIndelCnt " + str(job['input']['min_indel_count'])
    if job['input']['non_reference_probability_model'] != "EXACT":
        if job['input']['non_reference_probability_model'] != "GRID_SEARCH":
            raise AppError("Option \"Probability Model\" must be either \"EXACT\" or \"GRID_SEARCH\". Found " + job['input']['non_reference_probability_model'] + " instead")
        command += " -pnrm " + str(job['input']['non_reference_probability_model'])
    
    command += " --num_threads " + str(cpu_count())
    command += " -L regions.interval_list"

    if job['input']['downsample_to_coverage'] != 250:
        command += " -dcov " + str(job['input']['downsample_to_coverage'])
    elif job['input']['downsample_to_fraction'] != 1.0:
        command += " -dfrac " + str(job['input']['downsample_to_fraction'])

    if job['input']['nondeterministic']:
        command += " -ndrs "

    if job['input']['calculate_BAQ'] != "OFF":
        if job['input']['calculate_BAQ'] != "CALCULATE_AS_NECESSARY" and job['input']['calculate_BAQ'] != "RECALCULATE":
            raise AppError("Option \"Calculate BAQ\" must be either \"OFF\" or or \"CALCULATE_AS_NECESSARY\" \"RECALCULATE\". Found " + job['input']['calculate_BAQ'] + " instead")
        command += "-baq " + job['input']['calculate_BAQ']
        if job['input']['BAQ_gap_open_penalty'] != 40.0:
            command += "-baqGOP " + str(job['input']['BAQ_gap_open_penalty'])
    if job['input']['no_output_SLOD']:
        command += "-nosl"
    
    #print command
    return command

def extractHeader(vcfFileName, elevatedTags):
    result = {'columns': '', 'tags' : {'format' : {}, 'info' : {} }, 'filters' : {}}
    for line in open(vcfFileName):
        tag = re.findall("ID=(\w+),", line)
        if len(tag) > 0:
          tagType = ''
          if line.count("FORMAT") > 0:
            tagType = 'format'
          elif line.count("INFO") > 0:
            tagType = 'info'
          elif line.count("FILTER") > 0:
            result['filters'][re.findall("ID=(\w+),")[0]] = re.findall('Description="(.*)"')[0]
      
          typ = re.findall("Type=(\w+),", line)
          if tagType != '':
            number = re.findall("Number=(\w+)", line)
            description = re.findall('Description="(.*)"', line)
            if len(number) == 0:
              number = ['.']
            if len(description) == 0:
              description = ['']
            if "format_"+tag[0] not in elevatedTags:
                result['tags'][tagType][tag[0]] = {'type':typ[0], 'description' : description[0], 'number' : number[0]}
        if line[0] == "#" and line[1] != "#":
          result['columns'] = line.strip()
        if line == '' or line[0] != "#":
            break
    return result

#
#
#def main():
#    
#    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'
#
#    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
#    numberMappedReads = mappingsTable.describe()['size']
#    mappingSchema = mappingsTable.describe()['columns']
#    
#    covariateTable = dxpy.new_dxgtable(mappingSchema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
#    recalibratedTable = dxpy.new_dxgtable(mappingSchema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
#    
#    covariateFile = dxpy.new_dxfile()
#
#    try:
#        ##NEED TO CHANGE TO NEW CONTIGSET
#        #contigSetId = mappingsTable.get_details()['originalContigSet']['$dnanexus_link']
#        #originalContigSet = mappingsTable.get_details()['original_contigset']
#        contigSetId = mappingsTable.get_details()['originalContigSet']['$dnanexus_link']
#        originalContigSet = mappingsTable.get_details()['originalContigSet']
#    except:
#        raise Exception("The original reference genome must be attached as a detail")
#
#    print "Unpacking Reference"
#    subprocess.check_call("contigset2fasta %s ref.fa" % (contigSetId), shell=True)
#    referenceSequence = dxpy.dxlink(dxpy.upload_local_file("ref.fa"))
#    print "Indexing Dictionary"
#    subprocess.check_call("java net.sf.picard.sam.CreateSequenceDictionary R=ref.fa O=ref.dict", shell=True)
#    referenceDictionary = dxpy.dxlink(dxpy.upload_local_file("ref.dict"))
#
#    commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['mark_duplicates_chunks'])
#    print commandList
#    
#    ##Add known sites vcf
#    print referenceSequence
#    print referenceDictionary
#    print covariateTable.get_id()
#    print covariateFile.get_id()
#    countCovariateInput = {
#        'reference_sequence': referenceSequence,
#        'reference_dict': referenceDictionary,
#        'covariate_mappings_table_id': covariateTable.get_id(),
#        'covariate_file' : covariateFile.get_id(),
#    }
#    covariateJobId = dxpy.new_dxjob(fn_input=countCovariateInput, fn_name="countCovariates").get_id()
#    
#    
#    
#    #Run Realigner target creator
#    realignerTargetFile = dxpy.new_dxfile()
#    realignerTargetCreatorInput = {
#        'recalibrated_table_id': recalibratedTable.get_id(),
#        'reference_sequence': referenceSequence,
#        'reference_dict': referenceDictionary,
#        'known_indels': job['input']['known_indels'],
#        'target_file': realignerTargetFile.get_id()
#    }
#    realignerTargetCreatorJobId = dxpy.new_dxjob(fn_input=realignerTargetCreatorInput, fn_name="realignerTargetCreator").get_id()
#    
#    ##Add known sites and indels
#    reduceInput = {}
#    for i in range(len(commandList)):
#        print commandList[i]
#        if len(commandList[i]) > 0:
#            print i
#            print mappingsTable.get_id()
#            print covariateTable.get_id()
#            print covariateFile.get_id()
#            print recalibratedTable.get_id()
#            print commandList[i]
#            #print covariateJobId
#            mapMarkDuplicatesInput = {
#                'mappings_table_id': mappingsTable.get_id(),
#                'covariate_mappings_table_id': covariateTable.get_id(),
#                'covariate_file': covariateFile.get_id(),
#                'covariate_job': covariateJobId,
#                'target_creator_file': realignerTargetFile.get_id(),
#                'target_creator_job': realignerTargetCreatorJobId,
#                'recalibrated_table_id': recalibratedTable.get_id(),
#                'reference_sequence': referenceSequence,
#                'reference_dict': referenceDictionary,
#                'known_dbsnp': job['input']['known_dbsnp'],
#                'known_indels': job['input']['known_indels'],
#                'interval': commandList[i],
#                'job_number' : i
#            }
#            # Run a "map" job for each chunk
#            mapJobId = dxpy.new_dxjob(fn_input=mapMarkDuplicatesInput, fn_name="mapMarkDuplicates").get_id()
#            reduceInput["markDuplicatesJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}
#    reduceInput['dedupTableId'] = recalibratedTable.get_id()
#    reduceJob = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceMarkDuplicates")
#    reduceJobId = reduceJob.get_id()
#    reduceJob.wait_on_done()
#    #print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})
#    #job['output'] = {'dedup': {'job': reduceJobId, 'field': 'dedup'}}
#
#    
#def mapMarkDuplicates():
#    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar'
#
#    dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
#    dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
#
#    print job['input']['mappings_table_id'] 
#    print job['input']['interval']
#    
#    for i in range(len(job['input']['known_indels'])):
#        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output indels%d.vcf --region_index_offset -1 %s" % (job['input']['known_indels'][i], i, job['input']['interval']), shell=True)
#
#    print "Converting Table to SAM"
#    subprocess.check_call("dx-mappings-to-sam --table_id %s --output input.sam --region_index_offset -1 %s" % (job['input']['mappings_table_id'], job['input']['interval']), shell=True)
#    subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.bam O=input.sorted.sam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1 SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT", shell=True)
#    subprocess.call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sorted.sam O=dedup.sam REMOVE_DUPLICATES=true METRICS_FILE=metrics.txt ASSUME_SORTED=true", shell=True)
#    metricsFile = open("metrics.txt", 'r')
#    print metricsFile.read()
#
#    ##Get covariate information from the covariates job
#    job['input']['covariate_job'].wait_on_done()
#    dxpy.download_dxfile(job['input']['covariate_file'], 'recal.csv')
#
#    #Convert from SAM to BAM for GATK steps, Run TableRecalibration
#    print "Converting to BAM"
#    subprocess.call("java -Xmx4g net.sf.picard.sam.SamFormatConverter INPUT=dedup.sam OUTPUT=dedup.bam", shell=True)
#    
#    subprocess.check_call("java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T TableRecalibration -R ref.fa -I dedup.rg.bam -o  recal.bam -recalFile recal.csv", shell=True)
#    
#    ##Get intervals from the target creator job
#    job['input']['target_creator_job'].wait_on_done()
#    dxpy.download_dxfile(job['input']['target_creator_file'], 'indels.intervals')
#    
#    ##Run Indel Realigner
#    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T IndelRealigner -targetIntervals indels.intervals -R ref.fa -I recal.bam -o realign.bam"
#    for i in range(len(job['input']['known_indels'])):
#        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output indels%d.vcf" % (job['input']['known_indels'][i], i), shell=True)
#        command += " -known indels%d.vcf" % (i)
#    subprocess.check_call(command, shell=True)
#    subprocess.call("java -Xmx4g net.sf.picard.sam.SamFormatConverter INPUT=realign.bam OUTPUT=realign.sam", shell=True)
#    
#    
#    job["output"]["ok"] = True
#    
#
#def countCovariates():
#    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar'
#    
#    covariateFile = dxpy.open_dxfile(job['input']['covariate_file'])
#    dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
#    dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
#    
# 
#    for i in range(len(job['input']['known_dbsnp'])):
#        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output dbsnp%d.vcf" % (job['input']['known_dbsnp'][i], i), shell=True)
#        
#    
#    
#    startTime = time.clock()
#    while 1:
#        print time.clock() - startTime
#        if time.clock() - startTime >= 60*60*24*7:
#            break
#        elif job['input']['covariate_mappings_table_id'] == True:
#            ##Do count covariate stuff
#            subprocess.check_call("dx-mappings-to-sam --table_id %s --output input.sam" % (job['input']['covariate_mappings_table_id']), shell=True)
#            subprocess.call("java -Xmx4g net.sf.picard.sam.SamFormatConverter INPUT=input.sam OUTPUT=input.bam", shell=True)
#            subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.bam O=input.rg.bam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1", shell=True)
#            command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T CountCovariates -R ref.fa -I input.rg.bam -recalFile recal.csv -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --standard_covs"
#            for i in range(len(job['input']['known_dbsnp'])):
#                command += " -knownSites dpsnp%d.vcf" % (i)
#            subprocess.check_call(command, shell=True)
#            print "Finished Count Covariates"
#            break
#        time.sleep(1)
#    job['output']['ok'] = True
#        
#def realignerTargetCreator():
#    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'
#    
#    dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
#    dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
#    
#    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T RealignerTargetCreator -R ref.fa -o indels.intervals"
#    
#    for i in range(len(job['input']['known_indels'])):
#        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output indels%d.vcf" % (job['input']['known_indels'][i], i), shell=True)
#        command += " --known indels%d.vcf" % (i)
#    subprocess.check_call(command, shell=True)
#    
#    results = open('indels.intervals', 'r').read()
#    job['input']['target_file'].write(results)
#    job['input']['tartget_file'].close()
#    job['output']['ok'] = True
#    
#    
#
#
#
#def reduceMarkDuplicates():
#    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'
#    
#
#def reduceGatk():
#    #t = dxpy.open_dxgtable(job['input']['dedupTableId'])
#    #t.close(block=True)
#    #print "Closing Table"
#    #job['output']['dedup'] = dxpy.dxlink(t.get_id())
#    job["output"]["ok"] = True
#
#    
#def splitGenomeLengthLargePieces(contig_set, chunks):
#    details = dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
#    sizes = details['contigs']['sizes']
#    names = details['contigs']['names']
#    offsets = details['contigs']['offsets']
#
#    commandList = []
#    for i in range(chunks):
#        commandList.append('')
#    position = 0
#    chromosome = 0
#    chunkSize = sum(sizes)/chunks
#    currentChunk = 0
#    currentLength = 0
#    
#    while chromosome < len(names):
#        if position + (chunkSize - currentLength) >= sizes[chromosome]:
#            print chromosome
#            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, sizes[chromosome])
#            currentLength += sizes[chromosome] - position
#            chromosome += 1
#            position = 0
#        else:
#            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, position+chunkSize+1)
#            position += chunkSize + 1
#            if currentChunk < chunks-1:
#                currentChunk += 1
#            currentLength = 0
#    return commandList
#
#def checkIntervalRange(includeList, chromosome, lo, hi):
#    included = False
#    command = ''
#    if len(includeList) == 0:
#        return " -L %s:%d-%d" % (chromosome, lo, hi)
#    if includeList.get(chromosome) != None:
#        for x in includeList[chromosome]:
#            print "List"
#            print x
#            min = lo
#            max = hi
#            if (lo >= x[0] and lo <= x[1]) or (hi <= x[1] and hi >= x[0]):
#                if lo >= x[0] and lo <= x[1]:
#                    min = lo
#                elif lo <= x[0]:
#                    min = x[0]
#                if hi <= x[1] and hi >= x[0]:
#                    max = hi
#                elif hi >= x[1]:
#                    max = x[1]
#                command += " -L %s:%d-%d" % (chromosome, min, max)
#    return command
#
#        
