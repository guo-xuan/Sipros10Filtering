'''
Created on Jun 24, 2016

@author: xgo
'''
import sys
import getopt

# # Help message
help_message = '''
Usage:
    python Ft2ToMs2.py [options]

Inputs:
    FT2 file
    output ms2 filename

Options:
    -h/--help
    -v/--version
    -i/--input-file ./path    
    -o/--output-file ./path

Outputs:
    output file
'''

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:o:",
                                    ["help",
                                     "version",
                                     "input-file",
                                     "output-file"])

    # Default working dir and config file
    input_file = "./"
    output_file = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "Ft2ToMs2.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input-file"):
            input_file = value
        if option in ("-o", "--output-file"):
            output_file = value
            
    if input_file == "./" or output_file == "":
        print help_message
        sys.exit(0)

    return (input_file, output_file)

def extract(input_file, output_file):
    outputFile = open(output_file, "w")
    sLine = ""
    lSplit = []
    with open(input_file) as inputFile:
        for sLine in inputFile:
            if sLine[0].isdigit():
                lSplit = sLine.split()
                outputFile.write(lSplit[0])
                outputFile.write("\t")
                outputFile.write(lSplit[1])
                outputFile.write("\n")
            else:
                outputFile.write(sLine)  
    
    outputFile.close()
    print "DONE"

def main(argv=None):
    if argv is None:
        argv = sys.argv
    # parse options
    (input_file, output_file) = parse_options(argv)
    # process file
    extract(input_file, output_file)


if __name__ == '__main__':
    sys.exit(main())
