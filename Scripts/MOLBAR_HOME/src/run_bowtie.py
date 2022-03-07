'''
Run the Bowtie2 program one time.
'''
import argparse
import traceback
import os
import glob
import errno
import logging
from pathlib import Path
import science
import config

class BowtieRunner:

    def __init__(self, prog, options, debug, work_dir=''):
        self.work_dir = work_dir
        self.debug = debug
        self.program = prog
        self.options = options
        self.log_filename = 'run_bowtie.log'
        self.logger = self.get_logger(self.log_filename)
        self.logger.info('-------------------------------------------')
        self.logger.info('Running run_bowtie.py')
        self.logger.info('Program option is '+prog)
        self.logger.info('Options setting is: '+options)
        self.logger.info('Debug setting is: '+ str(debug))
        #self.conf = config.Config()
        #self.alignment = self.conf.get_value_by_key('bowtie2', 'options')

    def __del__(self):
        logging.shutdown()
        self.logger.handlers.clear()


    def get_index(self,target):
        pattern='*.bt2'  # default, bowtie2
        if (self.program=="hisat2"):
            pattern='*.ht2'
        fileList = glob.glob(target+pattern)
        if len(fileList) > 0:
            self.logger.info('Found the bowtie index: '+target)
            return target
        self.logger.info('Bowtie index not found. Build it now.')
        self.build_index(target)
        return target

    def build_index(self, target):
        fastafilename = target+'.fasta'
        self.logger.error('Building an index of: '+fastafilename)
        self.assert_file_exists(fastafilename)
        # Bowtie writes lots of stats to stdout and stderr.
        # We use bash to redirect both (2>&1) to append to our log (>>).
        self.build_command = 'bowtie2-build '
        if (self.program=="hisat2"):
            self.build_command = 'hisat2-build '
        self.build_command = self.build_command \
            +fastafilename+' ' \
            +target \
            +' >>'+self.log_filename \
            +' 2>&1 '
        self.logger.info(self.build_command)
        retval=os.system(self.build_command)
        if (retval==0):
            self.logger.info('Built index of: '+fastafilename)
        else:
            self.logger.warning('Possibly failed to build index!')
            self.logger.warning('Return value: ' + str(retval))

    def bowtie_align(self, index, r1, r2, base):
        # It is done in mapping.py: take the options from the config.db file.
        # TO DO: put these options in the config.db file.
        # DONE: make sure the option --sensitive is included.
        # DONE: make sure the option --end-to-end is included.
        # TO DO: either way, test the job success and log it.
        self.index = index
        self.r1 = r1
        self.r2 = r2
        # set output filename
        self.base = self.work_dir + base
        self.assert_file_exists(self.r1)
        self.assert_file_exists(self.r2)
        self.logger.info("Begin align.")
        # Bowtie reports mapping percentages to stderr.
        # Use bash to redirect stderr to the log (2>) in append mode (>>).
        self.align_command = 'bowtie2 '
        if (self.program=="hisat2"):
            self.align_command = 'hisat2 '
        self.align_command = self.align_command + self.options \
            + ' -x ' + index \
            + ' -1 ' + r1 + ' -2 '+ r2 \
            + ' -S ' + base \
            + ' 2>> ' + self.log_filename
        self.logger.info(self.align_command)
        #Check if bowtie returns a non zero
        if(os.system(self.align_command)):
            #print("ERROR! alignment failed, check "+self.log_filename)
            #print('Run with --debug for stack trace')
            self.logger.error("Aligner has returned an error")
            if self.debug:
                print(traceback.format_exc())
            return
        self.logger.info("Alignment finished")
        self.logger.info("Logging information to science.db")
        self.write_science_log()

    def get_logger(self,filename):
        logger = logging.getLogger(filename)
        logger.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s: %(levelname)s: %(message)s')
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        return logger

    def write_science_log(self):
        dbms = science.ScienceDB()
        in1 = self.index
        in2 = self.r1
        in3 = self.r2
        out1 = self.work_dir + self.base #Output file

        infile1 = dbms.add_data_file(in1,"Index","","Index")
        infile2 = dbms.add_data_file_auto_date(in2,"Reads 1 of pair","fastq")
        infile3 = dbms.add_data_file_auto_date(in3,"Reads 2 of pair","fastq")
        outfile1 = dbms.add_data_file_auto_date(out1,"SAM","Read alignments")

        infiles=[[in1,infile1],[in2,infile2],[in3,infile3]]
        outfiles=[[out1,outfile1]]

        dbms.add_process("align",''.join(self.align_command),infiles,outfiles)
        self.logger.info("Successfully logged information to science.db")

    def assert_file_exists (self, filename):
        if not os.path.isfile(filename):
            self.logger.error('File not found: '+filename)
            raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), filename)

if __name__ == '__main__':
    try:
        print("This program is a wrapper for bowtie2 and hisat2.")
        print("This program is part of a pipeline.")
        print("It does not support direct invocation from the command line.")
        # TO DO: Come on, support a command-line invocation!
    except Exception as e:
        print("\nERROR! run_bowtie.py has crashed\n")
        if args.debug:
            print(traceback.format_exc())
        else:
            print ('Run with --debug for traceback.')

