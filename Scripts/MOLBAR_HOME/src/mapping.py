'''
Maps reads to reference/target file
Assumes end to end pairing
'''
import argparse
import traceback
import os
import logging
from run_bowtie import BowtieRunner
import science
import fnmatch
from pathlib import Path

class RunMapping(object):

    def __init__(self, prog, hetero, target, r1, r2, output, debug):
        self.program = prog
        self.debug = debug
        self.log_filename = 'mapping.log'
        self.logger = self.get_logger(self.log_filename)
        self.logger.info('------------mapping.py----------------')
        self.logger.info(self.program)
        # TO DO: make it possible to change the options via a method call from outside this class.
        self.options = " --no-unal --no-mixed --no-discordant --sensitive --end-to-end --threads 4 "
        # In heterozygous mode, retain best 2 maps per read.
        # Else (homozygous mode), retain best 1 map per read.
        if hetero:
            self.options = self.options + ' -k 2 '
        else:
            self.options = self.options + ' -k 1 '
        self.logger.info(self.options)
        self.r1 = self.get_file(r1)
        self.r2 = self.get_file(r2)
        self.output = output
        self.bowtie = BowtieRunner(self.program,self.options,debug)
        self.ref_index = None
        try:
            self.ref_index = self.bowtie.get_index(target)
        except Exception as e:
            self.logger.error('Failed to find or create bowtie index: ' + target)
        self.logger.info(self.ref_index)

    def __del__(self):
        logging.shutdown()
        if self.logger and self.logger.handlers:
            self.logger.handlers.clear()

    def get_file(self,file):
        if os.path.isfile(file):
            return file
        else:
            self.logger.error('File not found: ' +file)
            raise Exception

    def get_logger(self, filename):
        logger = logging.getLogger(filename)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter(
            '%(asctime)s: %(levelname)s: %(message)s')
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return logger

    def mapIt(self):
        self.logger.info('Call bowtie_align(')
        self.logger.info('ref_index='+self.ref_index)
        self.logger.info(',r1='+self.r1)
        self.logger.info(',r2='+self.r2)
        self.logger.info(',output='+self.output+')')
        self.bowtie.bowtie_align(self.ref_index,
            self.r1, self.r2, self.output)
        # TO DO: detect whether the mapper actually finished.
        self.logger.info('Mapping finished')

    def main(self):
        self.mapIt()

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description = 'Align paired reads to reference.')
    parser.add_argument(
        '--program', help='hisat2 or bowtie2 (default)',
        type=str)
    parser.add_argument(
        '--heterozygous',
        help='Output two best alignments per read pair.', action='store_true')
    parser.add_argument('target', help='Reference (*.fasta)',
        type=str)
    parser.add_argument('read1', help='Read 1 (*.fastq).',
        type=str)
    parser.add_argument('read2', help='Read 2 (*.fastq)',
        type=str)
    parser.add_argument('output', help='Output (*.sam).',
        type=str)
    parser.add_argument('--debug', action='store_true', help='print traceback')
    args = parser.parse_args()

if __name__ == '__main__':
    try:
        args_parse()
        prog = args.program
        if (not prog):
            prog = "bowtie2"  # default
        instance = RunMapping(prog,
            args.heterozygous, args.target,
            args.read1, args.read2, args.output, args.debug)
        instance.main()
    except Exception as e:
        instance.logger.error(e)
        print('ERROR! check mapping.log')
        if args.debug:
            print(traceback.format_exc())
        else:
            print('ERROR, run with --debug for traceback')
        exit(1)
