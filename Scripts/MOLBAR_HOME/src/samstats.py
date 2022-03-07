'''
Statistics from SAM/BAM files.
This was used to comare bowtie, hisat, star aligners.
Inputs do not need to be sorted or indexed but must be SAM or BAM format.
By default, print summary for every *.bam *.sam in current directory.
Summary shows filename, num_aligns, num_no_indel, num_spliced.
Optionally, print detail for one file given on command line.
Detail shows target, num_aligns, num_no_indel, num_spliced.
Lack of indels is assumed if num cigar tuples is 1.
Lack of splice is detected by CICAR does not contain N.
Write to STDOUT.
'''
import argparse
import traceback
import sys
import os
import glob
import fnmatch
from pathlib import Path
import pysam

class SamStats (object):

    def __init__(self, sf, debug):
        self.samfile = sf
        self.debug = debug
        # pysam stores bam header in this data structure.
        # pysam rejects sam files missing a header.
        self.wrapper = pysam.AlignmentFile(self.samfile, 'r')

    def __del__(self):
        try:
            self.wrapper.close()
        except:
            pass

    def get_file(self,file):
        if os.path.isfile(file):
            return file
        else:
            print('File not found: ' +file)
            raise Exception

    def is_indel_free(self,alm):
        # If the number of cigar tuples is 1,
        # then the 1 tuple must be like 151M for 151 matches,
        # and there are no indels.
        if len(alm.cigartuples)==1:
            return True
        return False

    def is_spliced(self,alm):
        # If the cigar string contains an N,
        # then it skips over some number of bases denoted after the N,
        # and the alignment is considered spliced.
        # Note some aligners output small gaps like N5
        # that are clearly not splice sites in the biological sense.
        if 'N' in alm.cigarstring:
            return True
        return False

    def count_alignments(self):
        # cnt = self.wrapper.count()   # requires indexed bam
        self.allcount = 0
        self.noindel = 0
        self.spliced = 0
        for alm in self.wrapper:
            self.allcount += 1
            #print (str(type(alm)))
            if alm.cigarstring:
                if self.is_indel_free(alm):
                    self.noindel += 1
                if self.is_spliced(alm):
                    self.spliced += 1

    def load_counts_per_target(self):
        self.cnt_per_target={}
        self.noindel_per_target={}
        self.spliced_per_target={}
        cpt = self.cnt_per_target
        npt = self.noindel_per_target
        spt = self.spliced_per_target
        for alm in self.wrapper:
            rn = alm.reference_name
            if not rn in cpt:
                # Assume all three dicts are always in sync.
                cpt[rn] = 0
                npt[rn] = 0
                spt[rn] = 0
            cpt[rn] += 1
            if self.is_indel_free(alm):
                npt[rn] += 1
            if self.is_spliced(alm):
                spt[rn] += 1

    def print_summary(self):
        self.count_alignments()
        print("File=%s, Pairs=%d, NoIndel=%d, Spliced=%d" %
              (self.samfile,self.allcount,self.noindel,self.spliced))

    def print_details(self):
        self.load_counts_per_target()
        cpt = self.cnt_per_target
        npt = self.noindel_per_target
        spt = self.spliced_per_target
        print("%s\t%s\t%s\t%s\t%s" % 
              ("gene","allele","pairs","indel-free","spliced"))
        for xx in sorted(cpt):
            (gene,allele) = xx.split("_")
            print("%s\t%s\t%d\t%d\t%d" % 
                  (gene,allele,cpt[xx],npt[xx],spt[xx]))

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description = 'Generate statistics from SAM/BAM file.')
    parser.add_argument('--file', help='One file to process in detail.',
        type=str)
    parser.add_argument('--debug', help='print traceback',
        action='store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    try:
        args_parse()
        if args.file:
            instance = SamStats(args.file,args.debug)
            sys.stderr.write('Generate detail stats for %s\n' % args.file)
            sys.stderr.flush()
            instance.print_details()
        else:
            for filename in glob.glob('*.[bs]am'):
                instance = SamStats(filename,args.debug)
                sys.stderr.write('Generate summary stats for %s\n' % filename)
                sys.stderr.flush()
                instance.print_summary()
    except Exception as e:
        if args.debug:
            print(traceback.format_exc())
        else:
            print('ERROR, run with --debug for traceback')
        exit(1)
