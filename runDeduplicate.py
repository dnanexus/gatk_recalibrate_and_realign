import dxpy
import subprocess, logging
import os
import sys
import re
import math
import operator
import time

def main():
    
    os.environ['CLASSPATH'] = '/opt/jar/MarkDuplicates.jar:/opt/jar/MergeSamFiles.jar'

    mappingsTable = dxpy.open_dxgtable(job['input']['mappings']['$dnanexus_link'])
    numberMappedReads = mappingsTable.describe()['size']
    mappingSchema = mappingsTable.describe()['columns']
    
    print mappingSchema
    
    dedupTable = dxpy.new_dxgtable(mappingSchema, indices=[dxpy.DXGTable.genomic_range_index("chr","lo","hi", 'gri')])

    try:
        contigSetId = mappingsTable.get_details()['original_contigset']['$dnanexus_link']
        originalContigSet = mappingsTable.get_details()['original_contigset']
    except:
        raise Exception("The original reference genome must be attached as a detail")

    commandList = splitGenomeLengthLargePieces(originalContigSet, job['input']['mark_duplicates_chunks'])
    print commandList
    
    reduceInput = {}
    for i in range(len(commandList)):
        print commandList[i]
        if len(commandList[i]) > 0:
            print i
            print mappingsTable.get_id()
            print dedupTable.get_id()
            print commandList[i]
            #print covariateJobId
            mapMarkDuplicatesInput = {
                'mappings_table_id': mappingsTable.get_id(),
                'dedup_table_id': dedupTable.get_id(),
                'interval': commandList[i],
                'discard_duplicates': job['input']['discard_duplicates'],
                'job_number' : i
            }
            # Run a "map" job for each chunk
            mapJobId = dxpy.new_dxjob(fn_input=mapMarkDuplicatesInput, fn_name="mapMarkDuplicates").get_id()
            reduceInput["markDuplicatesJob" + str(i) + "TableId"] = {'job': mapJobId, 'field': 'id'}
    reduceInput['dedupTableId'] = dedupTable.get_id()
    reduceJob = dxpy.new_dxjob(fn_input=reduceInput, fn_name="reduceMarkDuplicates")
    reduceJobId = reduceJob.get_id()
    reduceJob.wait_on_done()
    print "SimpleVar table" + json.dumps({'table_id':dedupTable.get_id()})
    job['output'] = {'dedup': {'job': reduceJobId, 'field': 'dedup'}}

    
def mapMarkDuplicates():
    os.environ['CLASSPATH'] = '/opt/jar/AddOrReplaceReadGroups.jar:/opt/jar/MarkDuplicates.jar:/opt/jar/SamFormatConverter.jar:/opt/jar/SortSam.jar'

    print job['input']['mappings_table_id']
    print job['input']['dedup_table_id']
    print job['input']['interval']
    
    intervalMatch = re.findall("(\w+):(\d+)-(\d+)", job['input']['interval'])
    
    dedupTable = dxpy.open_dxgtable(job['input']['dedup_table_id'])
    mappingsTable = dxpy.open_dxgtable(job['input']['mappings_table_id'])
    
    
    
    #minId = -1
    #for x in re.findall("(\w+):(\d+)-(\d+)", job['input']['interval']):
    #    query = mappingsTable.genomic_range_query(chr=x[0], lo=int(x[1]), hi=int(x[2]))
    #    for row in mappingsTable.iterate_query_rows(query=query):
    #        minId = row[0]
    #    if minId > -1:
    #        break
    #
    #print minId
    #if minId > -1:
    
    print "Converting Table to SAM"
    subprocess.check_call("dx_mappingsTableToSam --table_id %s --output input.sam --region_index_offset -1 --output_ids%s" % (job['input']['mappings_table_id'], job['input']['interval']), shell=True)
    #subprocess.call("java -Xmx4g net.sf.picard.sam.AddOrReplaceReadGroups I=input.bam O=input.sorted.sam RGPL=illumina RGID=1 RGSM=1 RGLB=1 RGPU=1 SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT", shell=True)
    subprocess.call("java -Xmx4g net.sf.picard.sam.MarkDuplicates I=input.sam O=dedup.sam METRICS_FILE=metrics.txt ASSUME_SORTED=true", shell=True)
    metricsFile = open("metrics.txt", 'r')
    print metricsFile.read()
    
    colNames = dedupTable.get_col_names()
    col = {}
    for i in range(len(colNames)):
        col[colNames[i]] = i
        

    dedupFile = open("dedup.sam", 'r')
    for line in dedupFile:
        if line[0] != "@":
            tabSplit = line.split("\t")
            rowId = re.findall("ZD:Z:(\d+)", line)[0]
            if (int(tabSplit[1]) & 0x400) == False:
                entry = []
                row = mappingsTable.get_rows(starting=int(rowId), limit=1)['data'][0]
                for i in range(1,len(row)):
                    entry.append(row[i])
                dedupTable.add_rows([entry])
            elif job['input']['discard_duplicates'] == False:
                entry = []
                row = mappingsTable.get_rows(starting=int(rowId), limit=1)['data'][0]
                for i in range(1,len(row)):
                    entry.append(row[i])
                    entry[col["qc"]] = "PCR or optical duplicate"
                dedupTable.add_rows([entry])
                print entry
                print "Marked Duplicate"
                
def reduceMarkDuplicates():
    t = dxpy.open_dxgtable(job['input']['dedupTableId'])
    t.close(block=True)
    print "Closing Table"
    job['output']['dedup'] = dxpy.dxlink(t.get_id())


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

        
