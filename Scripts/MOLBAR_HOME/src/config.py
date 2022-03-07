import traceback
import argparse
import csv
import os

FILENAME = "config.csv"

#class Configuration:
#    def __init__(self,resource,attribute,value):
#        self.resource=resource
#        self.attribute=attribute
#        self.value=value

class Config:
    """
    Manage the list of environment settings.
    """

    def __del__(self):
        pass

    def __init__(self, dbName=FILENAME):
        self.filename=dbName

    def get_value(self,res,att):
        """
        Return the value from the first matching line.
        """
        val=None
        if not os.path.exists(self.filename):
            # TO DO: decide on exception or error logging
            print("WARNING: Missing file "+self.filename)
            return val
        with open(self.filename,"r") as ins:
            reader=csv.DictReader(ins,delimiter=',')
            for line in reader:
                if line['resource']==res:
                    if line['attribute']==att:
                        val=line['value']
        return val

    def add_value(self,res,att,val):
        """
        Append a line to the file. No check for duplicates.
        """
        line={'resource':res,'attribute':att,'value':val}
        fn=('resource','attribute','value')
        newfile = not os.path.exists(self.filename)
        with open(self.filename,"a+") as outs:
            writer=csv.DictWriter(outs,delimiter=',',fieldnames=fn)
            if newfile:
                writer.writeheader()
            writer.writerow(line)

    def update_value(self, res, att, up_val):
        """
        Updates the value should the user see fit.
        """
        match_found = False
        csv_file = open(self.filename, "r")
        lines = csv_file.readlines()
        line_counter = 0
        for line in lines:
            columns = line.split(',')
            if columns[0] == res and columns[1] == att:
                match_found = True
            else:
                line_counter += 1
            if match_found:
                lines[line_counter] = "{},{},{}\n".format(res, att, up_val)
        csv_file.close()
        csv_file = open(self.filename, "w")
        csv_file.writelines(lines)
        csv_file.close()

def args_parse():
    global args
    parser = argparse.ArgumentParser(
        description='Manage the config.db file.')
    parser.add_argument('resource', help='trimmomatic', type=str)
    parser.add_argument('attribute', help='jarpath', type=str)
    parser.add_argument('--add', help='/home', type=str)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--update', help='update value', type=str)
    args = parser.parse_args()  # on error, program exits

if __name__ == '__main__':
    args_parse()
    c = Config()
    try:
        if (args.add is not None):
            c.add_value(args.resource,args.attribute,args.add)
        elif (args.update is not None):
            c.update_value(args.resource, args.attribute, args.update)
        else:
            print(c.get_value(args.resource, args.attribute))
    except Exception as e:
        print("\nThere was an error.")
        if args.debug:
            print(traceback.format_exc())
        else:
            print('Run with --debug for traceback.')
