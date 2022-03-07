import argparse
import traceback
import os
import errno
import pathlib
import logging
import science
import numpy
import sys

class NormalizeCounts(object):

    def __init__(self, infile, outfile, debug):
        filename = os.path.basename(__file__)
        f = filename.split(".")
        prefix_filename = f[0]
        self.infile = infile
        self.output = outfile
        self.debug = debug
        self.DELIMITER=','  # TO DO: make this settable
        self.PSEUDOCOUNT=1   # TO DO: make this settable

    def assert_file_exists (self, filename):
        if not os.path.isfile(filename):
            self.logger.error('File not found: '+filename)
            raise FileNotFoundError(filename)

    def normalize(self):
        array_counts = self.file_to_array(self.infile)
        normalized = self.apply_normal(array_counts)
        if (self.output != None):
            numpy.savetxt(self.output, normalized, fmt='%s')
            self.assert_file_exists(self.output)
        else:
            numpy.savetxt(sys.stdout, normalized, fmt='%s')

    def file_to_array(self, input):
        """Converts a file to a 1D array with each line of the file a different value in the array"""
        self.assert_file_exists(input)
        a = []
        infile = open(input, 'r')
        for line in infile:
            line = line.strip() #removes \n for new line from the values
            a.append(line)
        infile.close()
        return a

    def parse_line(self,dataline):
        # Expect: gene,1,2,3,4,5,6
        parsed=[]
        fields = dataline.split(self.DELIMITER)
        f=0
        while f <= 6:   # process fields 1 through 6
            if f==0:
                parsed.append(fields[f])  # gene name
            else:
                val = int(fields[f])    # read count
                val += self.PSEUDOCOUNT
                parsed.append(val)
            f += 1
        return parsed

    def float_to_count(self,floatingpoint):
        s = str(int(round(floatingpoint)))
        return s

    def apply_normal(self, array_counts):
        normalized_array = []
        base = 0  # the normal case
        if self.PSEUDOCOUNT==0:   
            base = 1  # avoid divide by zero
        norm1_4=base
        norm2_5=base
        norm3_6=base
        for r in range(len(array_counts)):
            line = array_counts[r]
            fields = self.parse_line(line)
            norm1_4 += fields[1]+fields[4]
            norm2_5 += fields[2]+fields[5]
            norm3_6 += fields[3]+fields[6]
        total_reads =  norm1_4 + norm2_5 +  norm3_6
        average_reads = total_reads/3
        print("total_reads,%d"%total_reads)
        print("average_reads,%f"%average_reads)
        print("normalize_BR1,%f"%(average_reads/norm1_4))
        print("normalize_BR2,%f"%(average_reads/norm2_5))
        print("normalize_BR3,%f"%(average_reads/norm3_6))
        for r in range(len(array_counts)):
            line = array_counts[r]
            fields = self.parse_line(line)
            fields[1] = self.float_to_count(fields[1]*average_reads/norm1_4)
            fields[2] = self.float_to_count(fields[2]*average_reads/norm2_5)
            fields[3] = self.float_to_count(fields[3]*average_reads/norm3_6)
            fields[4] = self.float_to_count(fields[4]*average_reads/norm1_4)
            fields[5] = self.float_to_count(fields[5]*average_reads/norm2_5)
            fields[6] = self.float_to_count(fields[6]*average_reads/norm3_6)
            revised = self.DELIMITER.join(fields)
            normalized_array.append(revised)
        return(normalized_array)


def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description = 'Normalize counts for each of three replicates.')
    parser.add_argument('infile', help='Text file of gene + six counts per line.',
        type=str)
    parser.add_argument('outfile', help='Text file of gene + six counts per line.',
        type=str)
    parser.add_argument('--debug', action='store_true', help='print traceback')
    args = parser.parse_args()

if __name__ =="__main__":
    try:
        args_parse()
        instance = NormalizeCounts(args.infile, args.outfile, args.debug)
        instance.normalize()
    except Exception as e:
        print('ERROR!')
        if args.debug:
            print(traceback.format_exc())
        else:
            print('ERROR, run with --debug for traceback')
        exit(1)
