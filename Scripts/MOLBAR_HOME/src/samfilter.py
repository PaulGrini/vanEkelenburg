'''
Read a SAM file, write a BAM file.
Retain only primary alignments (i.e. best alignment per pair).
Retain only unspliced alignments.
This was used to compare bowtie, hisat, star.
'''
import argparse
import traceback
import os
import glob
import fnmatch
from pathlib import Path
import pysam

class SamFilter (object):

    def __init__(self, infile, outfile):
        self.infile=infile
        self.outfile=outfile
        # Use the read-only flag on the input file.
        # No need for the binary flag on the input file.
        # The pysam object will auto-detect sam/bam.
        self.reader = pysam.AlignmentFile(self.infile, 'r')
        # Use the read-write binary flags on the output file.
        # The b flag tells the pysam object to write bam not sam.
        self.writer = pysam.AlignmentFile(self.outfile, 'wb',
            template=self.reader)
        # The template allows re-production of the sam header.

    def __del__(self):
        if self.reader:
            self.reader.close()
        if self.writer:
            self.writer.close()

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

    def is_secondary(self,alm):
        # 'samtools view -f 0x100' shows only secondary alignments.
        # 'samtools view -F 0x100' shows only the primary alignment.
        secondary = False
        if alm.flag & 0x0100:
            secondary = True

    def main(self):
        aligns_in = 0
        aligns_out = 0
        for alm in self.reader:
            # How to print all the field names in class alm:
            # print(dir(alm))
            aligns_in += 1
            if not args.noseconds or not self.is_secondary(alm):
                if not args.noindel or self.is_indel_free(alm):
                    if not args.nosplice or not self.is_spliced(alm):
                        self.writer.write(alm)
                        aligns_out +=1
        print(aligns_in," alignments input")
        print(aligns_out," alignments output")

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description = 'Read SAM | filter | write BAM.')
    parser.add_argument('infile', type=str,
        help='Input SAM file (ex: in.sam or in.bam)')
    parser.add_argument('outfile', type=str,
        help='Output BAM file (ex: out.bam)')
    parser.add_argument('--noseconds', help='only primary alignments',
        action='store_true')
    parser.add_argument('--noindel', help='only indel-free alignments',
        action='store_true')
    parser.add_argument('--nosplice', help='only unspliced alignments',
        action='store_true')
    parser.add_argument('--debug', help='print traceback',
        action='store_true')
    args = parser.parse_args()

if __name__ == '__main__':
    try:
        args_parse()
        instance=SamFilter(args.infile,args.outfile)
        instance.main()
    except Exception as e:
        if args.debug:
            print(traceback.format_exc())
        else:
            print('ERROR, run with --debug for traceback')
        exit(1)
