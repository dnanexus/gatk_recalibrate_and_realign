#!/usr/bin/env python
import os, sys, unittest, json, subprocess

import dxpy, dxpy.app_builder
from dxpy.exceptions import *

src_dir = os.path.join(os.path.dirname(__file__), "..")
test_resources_dir = os.path.join(src_dir, "test", "resources")

from optparse import OptionParser

def makeInputsBwa():
    try:
        contigset_importer = dxpy.DXApplet(dxpy.find_data_objects(classname="applet", properties={"name": "fasta_contigset_importer"}).next()['id'])

    except StopIteration:
        raise Exception("fasta_contigset_importer or Letter Space FASTQ importer not found, please upload them")
    
    genome_archive = dxpy.upload_local_file(os.path.join(test_resources_dir, "hg19_chr22.fa.xz"), wait_on_close=True)
    contigset_importer_input = {"name": "hg19_chr22", "sequence_file": dxpy.dxlink(genome_archive)}
    print "Running fasta_contigset_importer with", contigset_importer_input
    job = contigset_importer.run(contigset_importer_input)
    job.wait_on_done()
    contig_set = job.describe()["output"]["contig_set"]
    
    return {"reference": contig_set}


class TestMyApp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
  
        bundled_resources = dxpy.app_builder.upload_resources(src_dir)
        program_id = dxpy.app_builder.upload_program(src_dir, bundled_resources, overwrite=True)
        cls.program = dxpy.DXApplet(program_id)

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_full_pipeline(self):
        #mappings = {"$dnanexus_link":mappingsId}
        #print {'mappings':mappings}
        job = self.program.run({'mappings':{"$dnanexus_link":"gtable-9yjpgx00000FvFXPg98Q0003"}, 'known_dbsnp':[{"$dnanexus_link":"gtable-9ykJkz800009pbgJgV8Q003G"}], 'known_indels':[{"$dnanexus_link":"gtable-9ykJx3Q00008k68jj58Q006V"}, {"$dnanexus_link":"gtable-9ykK16j00008k68jj58Q006j"}]})            
        job.wait_on_done()
        print "Recalibration output:"
        print json.dumps(job.describe()["output"])

if __name__ == '__main__':
    #If a mappings_id is provided, test suite will assume you want to resume to test at GATK with the output from bwa
    #parser = OptionParser("Usage: % mappings_id")
    #parser.add_option("--mappings_id", dest="mappings_id", default=False, help="gtable with the bwa output of a previous run")
    #(opts, args) = parser.parse_args()
    #mappingsId = opts.mappings_id
    #del sys.argv[1:]
    unittest.main()
