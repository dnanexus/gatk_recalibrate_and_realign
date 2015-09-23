# Copyright (C) 2013 DNAnexus, Inc.
#
# This file is part of gatk_recalibrate_and_realign (DNAnexus platform app).
#
# (The MIT Expat License)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the "Software"),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator
import time
import string

import dxpy
import subprocess, logging
import os, sys, re, math, operator, time
from multiprocessing import Pool, cpu_count

@dxpy.entry_point('main')
def main(**job_inputs):
    ## RUN DEDUPLICATE
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'

    #mappingsTable = dxpy.open_dxgtable(job_inputs['mappings'][0]['$dnanexus_link'])
    if 'output_name' in job_inputs:
        outputName = job_inputs['output_name']
    else:
        outputName = ''
    recalibratedTable = createNewTable(job_inputs['mappings'], outputName)

    #print "Mappings Table: " + mappingsTable.get_id()
    print "Recalibrated Table: " + recalibratedTable.get_id()

    mappingsTable = dxpy.DXGTable(job_inputs['mappings'][0]['$dnanexus_link'])
    for x in job_inputs['mappings']:
        if 'quality' not in dxpy.DXGTable(x).get_col_names():
            if len(job_inputs['mappings']) > 1:
                raise dxpy.AppError("One of the provided mappings did not have quality scores, for example %s. GATK can't recalibrate mappings without quality scores" % dxpy.DXGTable(x).describe()['name'])
            else:
                raise dxpy.AppError("The provided mappings did not have quality scores. GATK can't recalibrate mappings without quality scores")

    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise dxpy.AppError("The original reference genome must be attached as a detail")

    if contigSetId != job_inputs['reference']['$dnanexus_link']:
        raise dxpy.AppError("The reference genome of the mappings does not match the provided reference genome")

    reads = 0
    for x in job_inputs['mappings']:
        table = dxpy.DXGTable(x)
        reads += int(table.describe()['length'])
    chunks = int(reads/job_inputs['reads_per_job'])+1

    #Split the genome into chunks to parallelize
    commandList = splitGenomeLengthChromosome(originalContigSet, chunks)
    chunks = len(commandList)
    if chunks == 1:
        job_inputs['deduplicate_interchromosomal_pairs'] = False

    excludeInterchromosome = (chunks > 1)
    markDuplicatesJobs = []

    #This is a Picard Mark Duplicates job run only on interchromosomal mappings in the case that the genome is split into regions
    #This is necessary because Mark Duplicates has to look at both mates in a read pair, so interchromosomal mappings must go together
    reduceInterchromosomeInput = {}
    bamFiles = []
    if job_inputs['deduplicate_interchromosomal_pairs']:
        for i in xrange(-1, chunks):
            bamFiles.append(dxpy.new_dxfile().get_id())
            mapInterchromosomeInput = {
            'mappings_tables': job_inputs['mappings'],
            'interval': commandList[i],
            'job_number' : i,
            'separate_read_groups' : job_inputs['separate_read_groups']
            }
            interchromosomeJobId = dxpy.new_dxjob(fn_input=mapInterchromosomeInput, fn_name="mapInterchromosome").get_id()
            reduceInterchromosomeInput["mapJob"+str(i)] = {'job': interchromosomeJobId, 'field': 'file_id'}
        #interchromosomeJobField = { 'job': interchromosomeJobId, 'field': 'bam'}

        #Make a reduce job for the interchromosome component
        reduceInterchromosomeInput["file_list"] =  bamFiles
        reduceInterchromosomeInput["interval"] = commandList
        reduceInterchromosomeInput["discard_duplicates"] = job_inputs['discard_duplicates']
        reduceJobId = dxpy.new_dxjob(fn_input=reduceInterchromosomeInput, fn_name="reduceInterchromosome").get_id()
        deduplicateInterchromosome = True
    else:
        interchromosomeJobField = ''
        deduplicateInterchromosome = False

    #This runs the Picard Mark Duplicates program to deduplicate the reads
    reduceInput = {}
    for i in range(len(commandList)):
        print commandList[i]
        mapBestPracticesInput = {
            'mappings_tables': job_inputs['mappings'],
            'recalibrated_table_id': recalibratedTable.get_id(),
            'file_list': bamFiles,
            'interval': commandList[i],
            'job_number' : i,
            'reference': job_inputs['reference']['$dnanexus_link'],
            'dbsnp': job_inputs['dbsnp'],
            'separate_read_groups' : job_inputs['separate_read_groups'],
            'discard_duplicates': job_inputs['discard_duplicates'],
            'parent_input': job_inputs,
            'intervals_to_include': job_inputs.get('intervals_to_process'),
            'intervals_to_exclude': job_inputs.get('intervals_to_exclude'),
            'intervals_merging': job_inputs['intervals_merging'],
            'deduplicate_interchromosome': deduplicateInterchromosome
        }
        if 'known_indels' in job_inputs:
            mapBestPracticesInput['known_indels'] = job_inputs['known_indels']


        mapJobId = dxpy.new_dxjob(fn_input=mapBestPracticesInput, fn_name="mapBestPractices").get_id()
        reduceInput["mapJob" + str(i)] = {'job': mapJobId, 'field': 'ok'}
    reduceInput["recalibrated_table"] = recalibratedTable.get_id()

    reduceJobId = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceBestPractices").get_id()
    output = {'recalibrated_mappings': {'job': reduceJobId, 'field': 'recalibrated_table'}}
    #output = {'recalibrated_mappings': dxpy.dxlink(recalibratedTable.get_id())}

    return output

def runAndCatchGATKError(command, shell=True):
    # Added to capture any errors outputted by GATK
    try:
        subprocess.check_output(command, stderr=subprocess.STDOUT, shell=shell)
    except subprocess.CalledProcessError, e:
        print e 
        error = '\n'.join([l for l in e.output.splitlines() if l.startswith('##### ERROR MESSAGE:')])
        if error: 
            raise dxpy.AppError("App failed with GATK error. Please see logs for more information: {err}".format(err=error))
        else: 
            raise dxpy.AppInternalError("App failed with error. Please see logs for more information: {err}".format(err=e))         

@dxpy.entry_point('reduceBestPractices')
def reduceBestPractices(**job_inputs):
    startTime = time.time()
    recalibratedTable = dxpy.DXGTable(job_inputs['recalibrated_table'])
    recalibratedTable.close()
    print "Table closing completed in " + str(int((time.time()-startTime)/60)) + " minutes"
    output = {'recalibrated_table': dxpy.dxlink(recalibratedTable.get_id())}

    return output

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

@dxpy.entry_point('mapInterchromosome')
def mapInterchromosome(**job_inputs):
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'

    regionFile = open("regions.txt", 'w')
    regionFile.write(job_inputs['interval'])
    regionFile.close()

    jobNumber = job_inputs['job_number']
    if jobNumber == -1:
        for i in range(len(job_inputs['mappings_tables'])):
            mappingsTable = dxpy.DXGTable(job_inputs['mappings_tables'][i]['$dnanexus_link']).get_id()
            command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --id_as_name --only_unmapped --read_group_platform illumina --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)
    else:
        for i in range(len(job_inputs['mappings_tables'])):
            mappingsTable = dxpy.DXGTable(job_inputs['mappings_tables'][i]['$dnanexus_link']).get_id()
            command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --id_as_name --region_index_offset -1 --region_file regions.txt --read_group_platform illumina --only_interchromosomal_mate --id_prepend %d_" % (mappingsTable, i, i)
            print command
            subprocess.check_call(command, shell=True)

    output = {}
    readsPresent = False
    if len(job_inputs['mappings_tables']) == 1:
        if checkSamContainsRead("input.0.sam"):
            readsPresent = True
            subprocess.check_call("mv input.0.sam output.sam", shell=True)
    else:
        command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=output.sam USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        for i in range(len(job_inputs['mappings_tables'])):
            if checkSamContainsRead("input."+str(i)+".sam"):
                command += " INPUT=input."+str(i)+".sam"
                readsPresent = True
        runAndCatchGATKError(command, shell=True)
    if readsPresent:
        subprocess.check_call("samtools view -bS output.sam > output.bam", shell=True)
        fileId = dxpy.upload_local_file("output.bam").get_id()
        output['file_id'] = fileId
    else:
        output['file_id'] = ''

    return output

@dxpy.entry_point('reduceInterchromosome')
def reduceInterchromosome(**job_inputs):
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'
    output = {}

    i = -1
    command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=merge.bam USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
    noFiles = True
    while "mapJob"+str(i) in job_inputs:
        if job_inputs["mapJob"+str(i)] != '':
            command += " INPUT=input."+str(i+1)+".bam"
            dxpy.DXFile(job_inputs["mapJob"+str(i)]).wait_on_close()
            dxpy.download_dxfile(job_inputs["mapJob"+str(i)], "input."+str(i+1)+".bam")
            noFiles = False
        i+=1

    if noFiles:
        for x in job_inputs['file_list']:
            dxpy.DXFile(x).close()
    else:
        runAndCatchGATKError(command, shell=True)
        runAndCatchGATKError("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=merge.bam O=dedup.bam METRICS_FILE=metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=%s" % job_inputs['discard_duplicates'], shell=True)
        #subprocess.check_call("samtools view -bS dedup.sam > dedup.bam", shell=True)
        subprocess.check_call("samtools index dedup.bam", shell=True)
        for i in range(len(job_inputs['interval'])):
            try:
                startTime = time.time()
                chromosomeFile = open("chromosome.bam", 'w')
                chromosomeFile.close()
                fileId = dxpy.DXFile(job_inputs['file_list'][i])
                command = "samtools view -b dedup.bam " + job_inputs['interval'][i].replace(" -L ", " ") + " > chromosome.bam"
                print command
                subprocess.check_call(command, shell=True)
                if os.stat('chromosome.bam').st_size > 0:
                    dxpy.upload_local_file("chromosome.bam", use_existing_dxfile=fileId)
                    print "Result upload completed in " + str(int((time.time()-startTime)/60)) + " minutes"
                else:
                    print "This interval was empty of reads"
                    dxpy.DXFile(job_inputs['file_list'][i]).close()
            except:
                print "Could not write this interval"
                dxpy.DXFile(job_inputs['file_list'][i]).close()
    output['ok'] = True
    return output

@dxpy.entry_point('mapBestPractices')
def mapBestPractices(**job_inputs):
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar:/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MergeSamFiles.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/AddOrReplaceReadGroups.jar'
    output = {}

    #mappingsTable = dxpy.open_dxgtable(job_inputs['mappings_table_id'])

    jobNumber = job_inputs['job_number']

    if job_inputs['intervals_merging'] != "INTERSECTION" and ("intervals_to_include" in job_inputs) and job_inputs.get("intervals_to_include") != "":
        job_inputs['interval'] = splitUserInputRegions(job_inputs['interval'], job_inputs['intervals_to_include'], "-L")
        if job_inputs['interval'] == '':
            output['ok'] = True
            return output

    gatkRegionList = ""
    if job_inputs['intervals_merging'] == "INTERSECTION":
        if job_inputs.get('intervals_to_process') != None:
            gatkRegionList += " " + job_inputs['intervals_to_process']
    if job_inputs.get('intervals_to_exclude') != None:
        gatkRegionList += " " + job_inputs['intervals_to_exclude']


    regionFile = open("regions.txt", 'w')
    print job_inputs['interval']
    regionFile.write(job_inputs['interval'])
    regionFile.close()

    readGroups = 0
    print "Converting Table to SAM"
    for i in range(len(job_inputs['mappings_tables'])):
        mappingsTable = dxpy.DXGTable(job_inputs['mappings_tables'][i]['$dnanexus_link']).get_id()
        if job_inputs['deduplicate_interchromosome']:
            command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --region_index_offset -1 --id_as_name --region_file regions.txt --no_interchromosomal_mate --write_row_id --read_group_platform illumina --id_prepend %d_" % (mappingsTable, i, i)
        else:
            command = "pypy /usr/bin/dx_mappings_to_sam2 %s --output input.%d.sam --region_index_offset -1 --id_as_name --region_file regions.txt --write_row_id --read_group_platform illumina --id_prepend %d_" % (mappingsTable, i, i)
        if job_inputs['separate_read_groups']:
            command += " --add_to_read_group " + str(readGroups)
            readGroups += len(dxpy.DXGTable(job_inputs['mappings_tables'][i]['$dnanexus_link']).get_details()['read_groups'])
        print command
        startTime = time.time()
        subprocess.check_call(command, shell=True)
        print "Download mappings completed in " + str(int((time.time()-startTime)/60)) + " minutes"

    readsPresent = False

    if len(job_inputs['mappings_tables']) == 1:
        if checkSamContainsRead("input.0.sam"):
            readsPresent = True
            subprocess.check_call("mv input.0.sam input.sam", shell=True)
    else:
        command = "java -Xmx4g net.sf.picard.sam.MergeSamFiles OUTPUT=input.sam USE_THREADING=true SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
        for i in range(len(job_inputs['mappings_tables'])):
            if checkSamContainsRead("input."+str(i)+".sam"):
                command += " INPUT=input."+str(i)+".sam"
                readsPresent = True

    if readsPresent:
        runAndCatchGATKError(command, shell=True)
        runAndCatchGATKError("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sam O=dedup.sam METRICS_FILE=metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=%s" % job_inputs["discard_duplicates"], shell=True)
        startTime = time.time()
        subprocess.check_call("samtools view -bS dedup.sam > dedup.bam", shell=True)
        print "Conversion to BAM completed in " + str(int((time.time()-startTime)/60)) + " minutes"
        #else:
        #    subprocess.check_call("java net.sf.picard.sam.AddOrReplaceReadGroups INPUT=dedup.bam OUTPUT=dedup.rg.bam RGID="+str(job_inputs['job_number'])+" RGLB=library RGPL=illumina RGPU="+str(job_inputs['job_number'])+" SM="+str(job_inputs['job_number']), shell=True)
    else:
        output['ok'] = True
        return output

    ##THIS SEGMENT PROCEEDS AFTER MARKDUPLICATES COMPLETES
    sleepTime = 10
    sleepCounter = 0
    startTime = time.time()
    if jobNumber >  -1 and job_inputs['deduplicate_interchromosome']:
        if job_inputs['file_list'][jobNumber] != '':
            interchromosomeBamId = job_inputs['file_list'][jobNumber]
            interchromosomeBam = dxpy.DXFile(interchromosomeBamId)
            while 1:
                if interchromosomeBam.describe()['state'] != 'closed':
                    time.sleep(sleepTime)
                    sleepCounter += sleepTime
                else:
                    print "Waited " + str(sleepTime/60) + " minutes for interchromosome job"
                    dxpy.download_dxfile(interchromosomeBamId, "interchromosomeBam.bam")
                    print "Interchromosome BAM: " + interchromosomeBam.get_id()
                    print "Download interchromosome BAM in " + str(int((time.time()-startTime)/60)) + " minutes"
                    break



    print "dedup file:" + dxpy.upload_local_file("dedup.bam").get_id()

    if job_inputs['file_list'] == []:
        subprocess.check_call("mv dedup.bam input.bam", shell=True)
    elif job_inputs['file_list'][jobNumber] != '':
        runAndCatchGATKError("java -Xmx4g net.sf.picard.sam.MergeSamFiles SORT_ORDER=coordinate USE_THREADING=true INPUT=dedup.bam INPUT=interchromosomeBam.bam OUTPUT=input.bam VALIDATION_STRINGENCY=SILENT", shell=True)
    else:
        subprocess.check_call("mv dedup.bam input.bam", shell=True)
    startTime = time.time()
    subprocess.check_call("samtools index input.bam", shell=True)
    print "Index BAM completed in " + str(int((time.time()-startTime)/60)) + " minutes"

    #Download the Reference Genome
    print "Converting Contigset to Fasta"
    subprocess.check_call("dx-contigset-to-fasta %s ref.fa" % (job_inputs['reference']), shell=True)

    #RealignerTargetCreator
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T RealignerTargetCreator -R ref.fa -I input.bam -o indels.intervals %s -rf BadCigar" % gatkRegionList
    command += job_inputs['interval']

    # Download the known indels
    knownCommand = ''
    if 'known_indels' in job_inputs:
        for i in range(len(job_inputs['known_indels'])):
            dxpy.download_dxfile(job_inputs['known_indels'][i], "indels"+str(i)+".vcf.gz")
            knownFileName = "indels"+str(i)+".vcf.gz"
            try:
                p = subprocess.Popen("tabix -f -p vcf " + knownFileName, stderr=subprocess.PIPE, shell=True)
                if '[tabix] was bgzip' in p.communicate()[1]:
                    subprocess.check_call("zcat -f " + knownFileName + " | bgzip -c >temp.vcf.gz && mv -f temp.vcf.gz " + knownFileName + " && tabix -p vcf " + knownFileName, shell=True)
            except subprocess.CalledProcessError:
                raise dxpy.AppError("An error occurred while trying to index the provided known indels with tabix. Please make sure the provided known indels are valid VCF files.")
            knownCommand += " -known " + knownFileName
        command += knownCommand

    #Find chromosomes
    regionChromosomes = []
    for x in re.findall("-L ([^:]*):\d+-\d+", job_inputs['interval']):
        regionChromosomes.append(x)

    #Add options for RealignerTargetCreator
    if job_inputs['parent_input']['window_size'] != 10:
        command += "--windowSize " + str(job_inputs['parent_input']['window_size'])
    if job_inputs['parent_input']['max_interval_size'] != 500:
        command += " --maxIntervalSize "  + str(job_inputs['parent_input']['max_interval_size'])
    if job_inputs['parent_input']['min_reads_locus'] != 4:
        command += " --minReadsAtLocus " + str(job_inputs['parent_input']['min_reads_locus'])
    if job_inputs['parent_input']['mismatch_fraction'] != 0.0:
        command += " --mismatchFraction " + str(job_inputs['parent_input']['mismatch_fraction'])

    print command
    runAndCatchGATKError(command, shell=True)
    
    #Run the IndelRealigner
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T IndelRealigner -R ref.fa -I input.bam -targetIntervals indels.intervals -o realigned.bam %s -rf BadCigar" % gatkRegionList
    command += job_inputs['interval']
    command += knownCommand
    if "consensus_model" in job_inputs['parent_input']:
        if job_inputs['parent_input']['consensus_model'] != "":
            if job_inputs['parent_input']['consensus_model'] == "USE_READS" or job_inputs['parent_input']['consensus_model'] == "KNOWNS_ONLY" or job_inputs['parent_input']['consensus_model'] == "USE_SW":
                command += " --consensusDeterminationModel " + job_inputs['parent_input']['consensus_model']
            else:
                raise dxpy.AppError("The option \"Consensus Determination Model\" must be either blank or one of [\"USE_READS\", \"KNWONS_ONLY\", or \"USE_SW\"], found " + job_inputs['parent_input']['consensus_model'] + " instead.")
    if job_inputs['parent_input']['lod_threshold'] != 5.0:
        command += " --LODThresholdForCleaning " + str(job_inputs['parent_input']['lod_threshold'])
    if job_inputs['parent_input']['entropy_threshold'] != 0.15:
        command += " --entropyThreshold " + str(job_inputs['parent_input']['entropy_threshold'])
    if job_inputs['parent_input']['max_consensuses'] != 30:
        command += " --maxConsensuses " + str(job_inputs['parent_input']['max_consensuses'])
    if job_inputs['parent_input']['max_insert_size_movement'] != 3000:
        command += " --maxIsizeForMovement " + str(job_inputs['parent_input']['max_insert_size_movement'])
    if job_inputs['parent_input']['max_position_move'] != 200:
        command += " --maxPositionalMoveAllowed " + str(job_inputs['parent_input']['max_position_move'])
    if job_inputs['parent_input']['max_reads_consensus'] != 120:
        command += " --maxReadsForConsensus " + str(job_inputs['parent_input']['maxReadsForRealignment'])
    if job_inputs['parent_input']['max_reads_realignment'] != 20000:
        command += " --maxReadsForRealignment " + str(job_inputs['parent_input']['max_reads_realignment'])

    print command
    runAndCatchGATKError(command, shell=True)

    # Download dbSNP
    startTime = time.time()
    dxpy.download_dxfile(job_inputs['dbsnp'], "dbsnp.vcf.gz")
    print "Download dbsnp completed in " + str(int((time.time()-startTime)/60)) + " minutes"

    dbsnpFileName = 'dbsnp.vcf.gz'
    try:
        p = subprocess.Popen("tabix -f -p vcf dbsnp.vcf.gz", stderr=subprocess.PIPE, shell=True)
        if '[tabix] was bgzip' in p.communicate()[1]:
            subprocess.check_call("zcat -f dbsnp.vcf.gz | bgzip -c >temp.vcf.gz && mv -f temp.vcf.gz dbsnp.vcf.gz && tabix -p vcf dbsnp.vcf.gz", shell=True)
    except subprocess.CalledProcessError:
        raise dxpy.AppError("An error occurred while trying to index the provided dbSNP file with tabix. Please make sure the provided dbSNP file is a valid VCF file.")

    #subprocess.check_call("gzip -d dbsnp.vcf.gz", shell=True)

    #Count Covariates
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T CountCovariates -R ref.fa -recalFile recalibration.csv -I realigned.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --standard_covs %s -rf BadCigar" % gatkRegionList
    command += " -knownSites " + dbsnpFileName
    command += job_inputs['interval']
    if job_inputs['parent_input'].get('single_threaded') != True:
        command += " --num_threads " + str(cpu_count())

    if "solid_recalibration_mode" in job_inputs['parent_input']:
        command += " --solid_recal_mode " + job_inputs['parent_input']['solid_recalibration_mode']
    if "solid_nocall_strategy" in job_inputs['parent_input']:
        command += " --solid_nocall_strategy " + job_inputs['parent_input']['solid_nocall_strategy']
    if 'context_size' in job_inputs['parent_input']:
        command += " --context_size " + str(job_inputs['parent_input']['context_size'])
    if 'nback' in job_inputs['parent_input']:
        command += " --homopolymer_nback " + str(job_inputs['parent_input']['nback'])
    if job_inputs['parent_input']['cycle_covariate']:
        command += " -cov CycleCovariate"
    if job_inputs['parent_input']['dinuc_covariate']:
        command += " -cov DinucCovariate"
    if job_inputs['parent_input']['primer_round_covariate']:
        command += " -cov PrimerRoundCovariate"
    if job_inputs['parent_input']['mapping_quality_covariate']:
        command += " -cov MappingQualityCovariate"
    if job_inputs['parent_input']['homopolymer_covariate']:
        command += " -cov HomopolymerCovariate"
    if job_inputs['parent_input']['gc_content_covariate']:
        command += " -cov GCContentCovariate"
    if job_inputs['parent_input']['position_covariate']:
        command += " -cov PositionCovariate"
    if job_inputs['parent_input']['minimum_nqs_covariate']:
        command += " -cov MinimumNQSCovariate"
    if job_inputs['parent_input']['context_covariate']:
        command += " -cov ContextCovariate"

    print command
    runAndCatchGATKError(command, shell=True)

    #Table Recalibration
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T TableRecalibration -R ref.fa -recalFile recalibration.csv -I realigned.bam -o recalibrated.bam --doNotWriteOriginalQuals %s -rf BadCigar" % gatkRegionList
    command += job_inputs['interval']
    command += " --solid_recal_mode " + job_inputs['parent_input']['solid_recalibration_mode']
    command += " --solid_nocall_strategy " + job_inputs['parent_input']['solid_nocall_strategy']
    if 'context_size' in job_inputs['parent_input']:
        command += " --context_size " + str(job_inputs['parent_input']['context_size'])
    if 'nback' in job_inputs['parent_input']:
        command += " --homopolymer_nback " + str(job_inputs['parent_input']['nback'])
    if job_inputs['parent_input']['preserve_qscore'] != 5:
        command += " --preserve_qscores_less_than " + str(job_inputs['parent_input']['preserve_qscore'])
    if 'smoothing' in job_inputs['parent_input']:
        command += " --smoothing " + str(job_inputs['parent_input']['smoothing'])
    if 'max_quality' in job_inputs['parent_input']:
        command += " --max_quality_score " + str(job_inputs['parent_input']['max_quality'])
    print command
    runAndCatchGATKError(command, shell=True)

    subprocess.check_call("samtools view -h -o recalibrated.sam recalibrated.bam", shell=True)

    result = dxpy.upload_local_file("recalibrated.bam")
    output['recalibrated_bam'] = dxpy.dxlink(result.get_id())
    print "Recalibrated file: " + result.get_id()



    #Read SAM, extracting chr, lo, hi, qual, cigar, duplicate flag, chr2, lo2, hi2


    recalibratedTable = dxpy.DXGTable(job_inputs['recalibrated_table_id'])

    default = {}

    recalibratedColNames = recalibratedTable.get_col_names()
    recalibratedCol = {}
    for i in range(len(recalibratedColNames)):
        recalibratedCol[recalibratedColNames[i]] = i

    for x in recalibratedTable.describe()['columns']:
        if "int" in x["type"]:
            default[x["name"]] = sys.maxint
        elif x["type"] == "float":
            default[x["name"]] = float(sys.maxint)
        elif x["type"] == "boolean":
            default[x["name"]] = False
        else:
            default[x["name"]] = ""

    print "Writing mate pair information for lookup"
    startTime = time.time()
    if recalibratedCol.get("chr2") != None:
        mateLocations = {}
        for line in open("recalibrated.sam", 'r'):
            if line[0] != "@":
                tabSplit = line.split("\t")
                chr = tabSplit[2]
                lo = int(tabSplit[3])-1
                templateId = tabSplit[0]
                cigar = re.split('(\d+)', tabSplit[5])
                alignLength = 0
                for p in range(len(cigar)):
                    c = cigar[p]
                    if c == 'M' or c == 'D' or c == 'N' or c == 'X' or c == 'P' or c == '=':
                        alignLength += int(cigar[p-1])
                hi = lo + alignLength

                recalibrationTags = re.findall("zd:Z:([^\t]*)[\t\n]", line)[0].split("##&##")
                reportedLo = int(recalibrationTags[1])
                reportedHi = int(recalibrationTags[2])
                if lo != reportedLo or hi != reportedHi:
                    if int(tabSplit[1]) & 0x1 & 0x40:
                        mateLocations[templateId] = {0: {"lo":lo, "hi":hi, "chr":chr}}
                    elif int(tabSplit[1]) & 0x1 & 0x80:
                        mateLocations[templateId] = {1: {"lo":lo, "hi":hi, "chr":chr}}
        print str(len(mateLocations)) + " Interchromosomal reads changed lo or hi"
    print "Interchromosome changes to lo and hi recorded in " + str(int((time.time()-startTime)/60)) + " minutes"

    complement_table = string.maketrans("ATGCatgc", "TACGtacg")
    rowsWritten = 0

    startTime = time.time()
    for line in open("recalibrated.sam", 'r'):
        if line[0] != "@":
            tabSplit = line.split("\t")
            templateId = int(tabSplit[0].split("_")[1])
            flag = int(tabSplit[1])
            chr = tabSplit[2]
            lo = int(tabSplit[3])-1
            qual = tabSplit[10]
            alignLength = 0
            mapq = int(tabSplit[4])
            cigar = re.split('(\d+)', tabSplit[5])
            duplicate = (flag & 0x400 == True)
            sequence = tabSplit[9]
            recalibrationTags = re.findall("zd:Z:([^\t]*)[\t\n]", line)[0].split("##&##")

            name = recalibrationTags[0]
            readGroup = int(re.findall("RG:Z:(\d+)", line)[0])

            if flag & 0x4:
                status = "UNMAPPED"
            elif flag & 0x100:
                status = "SECONDARY"
            else:
                status = "PRIMARY"
            if flag & 0x200:
                qcFail = True
            else:
                qcFail = False
            if flag & 0x10 == 0 or status == "UNMAPPED":
                negativeStrand = False
            else:
                negativeStrand = True
                sequence = sequence.translate(complement_table)[::-1]
                qual = qual[::-1]

            for p in range(len(cigar)):
                c = cigar[p]
                if c == 'M' or c == 'D' or c == 'N' or c == 'X' or c == 'P' or c == '=':
                    alignLength += int(cigar[p-1])
            cigar = tabSplit[5]
            hi = lo + alignLength

            if recalibratedCol.get("chr2") != None:
                properPair=False
                if (flag & 0x1) and (flag & 0x2):
                    properPair = True

                reportedChr2 = recalibrationTags[3]
                reportedLo2 = int(recalibrationTags[4])
                reportedHi2 = int(recalibrationTags[5])
                status2 = recalibrationTags[6]

                if not flag & 0x1:
                    negativeStrand2 = False
                elif flag & 0x20 and status2 != "UNMAPPED":
                    negativeStrand2 = True
                else:
                    negativeStrand2 = False

                if flag & 0x1:
                    if flag & 0x40:
                        mateId = 0
                    elif flag & 0x80:
                        mateId = 1
                else:
                    mateId = -1

                try:
                    if mateId == 1:
                        chr2 = mateLocations[tabSplit[0]][0]["chr2"]
                        lo2 = mateLocations[tabSplit[0]][0]["lo"]
                        hi2 = mateLocations[tabSplit[0]][0]["hi"]
                    elif mateId == 0:
                        chr2 = mateLocations[tabSplit[0]][1]["chr2"]
                        lo2 = mateLocations[tabSplit[0]][1]["lo"]
                        hi2 = mateLocations[tabSplit[0]][1]["hi"]
                    else:
                        chr2 = reportedChr2
                        lo2 = reportedLo2
                        hi2 = reportedHi2
                except:
                    chr2 = reportedChr2
                    lo2 = reportedLo2
                    hi2 = reportedHi2
                recalibratedTable.add_rows([[sequence, name, qual, status, chr, lo, hi, negativeStrand, mapq, qcFail, duplicate, cigar, templateId, readGroup, mateId, status2, chr2, lo2, hi2, negativeStrand2, properPair]])
            else:
                recalibratedTable.add_rows([[sequence, name, qual, status, chr, lo, hi, negativeStrand, mapq, qcFail, duplicate, cigar, templateId, readGroup]])
            rowsWritten += 1
            if rowsWritten%100000 == 0:
                print "Imported " + str(rowsWritten) + " rows. Time taken: " + str(int((time.time()-startTime)/60)) + " minutes"
                recalibratedTable.flush()
    recalibratedTable.flush()
    output['ok'] = True

    return output


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
    totalSizes = []
    print chunkSizes
    for x in range(chunks):
        totalSizes.append(0)
    while chromosome < len(names):
        try:
            minimum = min([x for x in chunkSizes if x > 0])
            position = chunkSizes.index(minimum)
        except:
            position = 0
        if chunkSizes[position] + sizes[chromosome] > readsPerChunk:
            try:
                position = chunkSizes.index(0)
            except:
                position = chunkSizes.index(min([x for x in chunkSizes if x > 0]))
        commandList[position] += " -L %s:%d-%d" % (names[chromosome], 1, sizes[chromosome])
        chunkSizes[position] += sizes[chromosome]
        totalSizes[position] += sizes[chromosome]
        chromosome += 1
    commandList = filter(None, commandList)
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
    read_groups = []
    for i in range(len(mappingsArray)):
        oldTable = dxpy.DXGTable(mappingsArray[i]['$dnanexus_link'])
        if oldTable.get_details()['read_groups'] != None:
            read_groups.extend(oldTable.get_details()['read_groups'])
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

    print "combined columns"
    print columns

    schema = [{"name": "sequence", "type": "string"}]
    schema.append({"name": "name", "type": "string"})
    schema.append({"name": "quality", "type": "string"})
    schema.extend([{"name": "status", "type": "string"},
                          {"name": "chr", "type": "string"},
                          {"name": "lo", "type": "int32"},
                          {"name": "hi", "type": "int32"},
                          {"name": "negative_strand", "type": "boolean"},
                          {"name": "error_probability", "type": "uint8"},
                          {"name": "qc_fail", "type": "boolean"},
                          {"name": "duplicate", "type": "boolean"},
                          {"name": "cigar", "type": "string"},
                          {"name": "template_id", "type": "int64"},
                          {"name": "read_group", "type": "int32"}])
    if {"type":"string", "name":"chr2"} in columns:
        schema.extend([{"name": "mate_id", "type": "int32"}, # TODO: int8
                              {"name": "status2", "type": "string"},
                              {"name": "chr2", "type": "string"},
                              {"name": "lo2", "type": "int32"},
                              {"name": "hi2", "type": "int32"},
                              {"name": "negative_strand2", "type": "boolean"},
                              {"name": "proper_pair", "type": "boolean"}])

    oldTable = dxpy.DXGTable(mappingsArray[0]['$dnanexus_link'])
    if recalibratedName == '':
        recalibratedName = oldTable.describe()['name'] + " Realigned and Recalibrated"
    details = oldTable.get_details()
    details['read_groups'] = read_groups
    newTable = dxpy.new_dxgtable(columns=schema, indices=indices)
    newTable.add_tags(tags)
    newTable.set_details(details)
    newTable.add_types(types)
    newTable.rename(recalibratedName)

    return newTable

def splitUserInputRegions(jobRegions, inputRegions, prefix):
    jobList = re.findall("-L ([^:]*):(\d+)-(\d+)", jobRegions)
    inputList = re.findall("-L ([^:]*):(\d+)-(\d+)", inputRegions)
    
    result = ""
    for x in inputList:
        for y in jobList:
            if(x[0] == y[0]):
                lo = max(int(x[1]), int(y[1]))
                hi = min(int(x[2]), int(y[2]))
                if hi > lo:
                    result += " %s %s:%d-%d" % (prefix, x[0], lo, hi)
                    
    return result
