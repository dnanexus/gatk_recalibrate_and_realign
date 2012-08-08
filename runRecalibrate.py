import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator
import time

def main():
    
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    numberMappedReads = mappingsTable.describe()['size']
    mappingSchema = mappingsTable.describe()['columns']
    
    covariateTable = dxpy.new_dxgtable(mappingSchema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    recalibratedTable = dxpy.new_dxgtable(mappingSchema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])
    
    covariateFile = dxpy.new_dxfile()

    try:
        ##NEED TO CHANGE TO NEW CONTIGSET
        #contigSetId = mappingsTable.get_details()['originalContigSet']['$dnanexus_link']
        #originalContigSet = mappingsTable.get_details()['original_contigset']
        contigSetId = mappingsTable.get_details()['originalContigSet']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['originalContigSet']
    except:
        raise Exception("The original reference genome must be attached as a detail")

    print "Unpacking Reference"
    subprocess.check_call("contigset2fasta %s ref.fa" % (contigSetId), shell=True)
    referenceSequence = dxpy.dxlink(dxpy.upload_local_file("ref.fa"))
    print "Indexing Dictionary"
    subprocess.check_call("java net.sf.picard.sam.CreateSequenceDictionary R=ref.fa O=ref.dict", shell=True)
    referenceDictionary = dxpy.dxlink(dxpy.upload_local_file("ref.dict"))

    commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['mark_duplicates_chunks'])
    print commandList
    
    ##Add known sites vcf
    print referenceSequence
    print referenceDictionary
    print covariateTable.get_id()
    print covariateFile.get_id()
    countCovariateInput = {
        'reference_sequence': referenceSequence,
        'reference_dict': referenceDictionary,
        'covariate_mappings_table_id': covariateTable.get_id(),
        'covariate_file' : covariateFile.get_id(),
    }
    covariateJobId = dxpy.new_dxjob(fn_input=countCovariateInput, fn_name="countCovariates").get_id()
    
    
    
    #Run Realigner target creator
    realignerTargetFile = dxpy.new_dxfile()
    realignerTargetCreatorInput = {
        'recalibrated_table_id': recalibratedTable.get_id(),
        'reference_sequence': referenceSequence,
        'reference_dict': referenceDictionary,
        'known_indels': job['input']['known_indels'],
        'target_file': realignerTargetFile.get_id()
    }
    realignerTargetCreatorJobId = dxpy.new_dxjob(fn_input=realignerTargetCreatorInput, fn_name="realignerTargetCreator").get_id()
    
    ##Add known sites and indels
    reduceInput = {}
    for i in range(len(commandList)):
        print commandList[i]
        if len(commandList[i]) > 0:
            print i
            print mappingsTable.get_id()
            print covariateTable.get_id()
            print covariateFile.get_id()
            print recalibratedTable.get_id()
            print commandList[i]
            #print covariateJobId
            mapMarkDuplicatesInput = {
                'mappings_table_id': mappingsTable.get_id(),
                'covariate_mappings_table_id': covariateTable.get_id(),
                'covariate_file': covariateFile.get_id(),
                'covariate_job': covariateJobId,
                'target_creator_file': realignerTargetFile.get_id(),
                'target_creator_job': realignerTargetCreatorJobId,
                'recalibrated_table_id': recalibratedTable.get_id(),
                'reference_sequence': referenceSequence,
                'reference_dict': referenceDictionary,
                'known_dbsnp': job['input']['known_dbsnp'],
                'known_indels': job['input']['known_indels'],
                'interval': commandList[i],
                'job_number' : i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapMarkDuplicatesInput, fn_name="mapMarkDuplicates").get_id()
            reduceInput["markDuplicatesJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}
    reduceInput['dedupTableId'] = recalibratedTable.get_id()
    reduceJob = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceMarkDuplicates")
    reduceJobId = reduceJob.get_id()
    reduceJob.wait_on_done()
    #print "SimpleVar table" + json.dumps({'table_id':simpleVar.get_id()})
    #job['output'] = {'dedup': {'job': reduceJobId, 'field': 'dedup'}}

    
def mapMarkDuplicates():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar'

    dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")

    print job['input']['mappings_table_id'] 
    print job['input']['interval']
    
    for i in range(len(job['input']['known_indels'])):
        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output indels%d.vcf --region_index_offset -1 %s" % (job['input']['known_indels'][i], i, job['input']['interval']), shell=True)

    print "Converting Table to SAM"
    subprocess.check_call("dx-mappings-to-sam --table_id %s --output input.sam --region_index_offset -1 %s" % (job['input']['mappings_table_id'], job['input']['interval']), shell=True)
    subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.bam O=input.sorted.sam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1 SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT", shell=True)
    subprocess.call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sorted.sam O=dedup.sam REMOVE_DUPLICATES=true METRICS_FILE=metrics.txt ASSUME_SORTED=true", shell=True)
    metricsFile = open("metrics.txt", 'r')
    print metricsFile.read()

    ##Get covariate information from the covariates job
    job['input']['covariate_job'].wait_on_done()
    dxpy.download_dxfile(job['input']['covariate_file'], 'recal.csv')

    #Convert from SAM to BAM for GATK steps, Run TableRecalibration
    print "Converting to BAM"
    subprocess.call("java -Xmx4g net.sf.picard.sam.SamFormatConverter INPUT=dedup.sam OUTPUT=dedup.bam", shell=True)
    
    subprocess.check_call("java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T TableRecalibration -R ref.fa -I dedup.rg.bam -o  recal.bam -recalFile recal.csv", shell=True)
    
    ##Get intervals from the target creator job
    job['input']['target_creator_job'].wait_on_done()
    dxpy.download_dxfile(job['input']['target_creator_file'], 'indels.intervals')
    
    ##Run Indel Realigner
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T IndelRealigner -targetIntervals indels.intervals -R ref.fa -I recal.bam -o realign.bam"
    for i in range(len(job['input']['known_indels'])):
        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output indels%d.vcf" % (job['input']['known_indels'][i], i), shell=True)
        command += " -known indels%d.vcf" % (i)
    subprocess.check_call(command, shell=True)
    subprocess.call("java -Xmx4g net.sf.picard.sam.SamFormatConverter INPUT=realign.bam OUTPUT=realign.sam", shell=True)
    
    
    job["output"]["ok"] = True
    

def countCovariates():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar'
    
    covariateFile = dxpy.open_dxfile(job['input']['covariate_file'])
    dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
    
 
    for i in range(len(job['input']['known_dbsnp'])):
        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output dbsnp%d.vcf" % (job['input']['known_dbsnp'][i], i), shell=True)
        
    
    
    startTime = time.clock()
    while 1:
        print time.clock() - startTime
        if time.clock() - startTime >= 60*60*24*7:
            break
        elif job['input']['covariate_mappings_table_id'] == True:
            ##Do count covariate stuff
            subprocess.check_call("dx-mappings-to-sam --table_id %s --output input.sam" % (job['input']['covariate_mappings_table_id']), shell=True)
            subprocess.call("java -Xmx4g net.sf.picard.sam.SamFormatConverter INPUT=input.sam OUTPUT=input.bam", shell=True)
            subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.bam O=input.rg.bam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1", shell=True)
            command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T CountCovariates -R ref.fa -I input.rg.bam -recalFile recal.csv -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate --standard_covs"
            for i in range(len(job['input']['known_dbsnp'])):
                command += " -knownSites dpsnp%d.vcf" % (i)
            subprocess.check_call(command, shell=True)
            print "Finished Count Covariates"
            break
        time.sleep(1)
    job['output']['ok'] = True
        
def realignerTargetCreator():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'
    
    dxpy.download_dxfile(job['input']['reference_sequence'], "ref.fa")
    dxpy.download_dxfile(job['input']['reference_dict'], "ref.dict")
    
    command = "java -Xmx4g org.broadinstitute.sting.gatk.CommandLineGATK -T RealignerTargetCreator -R ref.fa -o indels.intervals"
    
    for i in range(len(job['input']['known_indels'])):
        subprocess.check_call("dx-simplevar-to-vcf --table_id %s --output indels%d.vcf" % (job['input']['known_indels'][i], i), shell=True)
        command += " --known indels%d.vcf" % (i)
    subprocess.check_call(command, shell=True)
    
    results = open('indels.intervals', 'r').read()
    job['input']['target_file'].write(results)
    job['input']['tartget_file'].close()
    job['output']['ok'] = True
    
    



def reduceMarkDuplicates():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/GenomeAnalysisTK.jar:/opt/jar/CreateSequenceDictionary.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'
    

def reduceGatk():
    #t = dxpy.open_dxgtable(job['input']['dedupTableId'])
    #t.close(block=True)
    #print "Closing Table"
    #job['output']['dedup'] = dxpy.dxlink(t.get_id())
    job["output"]["ok"] = True

    
def splitGenomeLengthLargePieces(contig_set, chunks):
    details = dxpy.DXRecord(contig_set['$dnanexus_link']).get_details()
    sizes = details['contigs']['sizes']
    names = details['contigs']['names']
    offsets = details['contigs']['offsets']

    commandList = []
    for i in range(chunks):
        commandList.append('')
    position = 0
    chromosome = 0
    chunkSize = sum(sizes)/chunks
    currentChunk = 0
    currentLength = 0
    
    while chromosome < len(names):
        if position + (chunkSize - currentLength) >= sizes[chromosome]:
            print chromosome
            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, sizes[chromosome])
            currentLength += sizes[chromosome] - position
            chromosome += 1
            position = 0
        else:
            commandList[currentChunk] += checkIntervalRange({}, names[chromosome], position+1, position+chunkSize+1)
            position += chunkSize + 1
            if currentChunk < chunks-1:
                currentChunk += 1
            currentLength = 0
    return commandList

def checkIntervalRange(includeList, chromosome, lo, hi):
    included = False
    command = ''
    if len(includeList) == 0:
        return " -L %s:%d-%d" % (chromosome, lo, hi)
    if includeList.get(chromosome) != None:
        for x in includeList[chromosome]:
            print "List"
            print x
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

        
