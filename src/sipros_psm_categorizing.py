'''
Created on Jul 20, 2016

@author: xgo
'''

import getopt, sys, os

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_psm_tabulating.py [options]

Inputs:
    tab file
    output directory

Options:
    -h/--help
    -v/--version
    -i/--input-folder ./path    
    -o/--output-folder ./path

Outputs:
    output PSM table
'''

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:o:",
                                    ["help",
                                     "version",
                                     "input-folder",
                                     "output-folder"])

    # Default working dir and config file
    input_file = ""
    output_folder = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "sipros_peptides_filtering.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input-folder"):
            input_file = value
        if option in ("-o", "--output-folder"):
            output_folder = value

    if input_file == "" or output_folder == "":
        print help_message
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, output_folder)

def categorize_pro_pep(val):
    if val == 3:
        return '3'
    elif val == 2:
        return '2'
    else:
        return '1'

def categorize_pep_psm(val):
    if val == 0:
        return '0'
    else:
        return '1'

def categorize_score_agreement(val):
    if val == 3:
        return '3'
    elif val == 2:
        return '2'
    else:
        return '1'

def categorize_parent_charge(val):
    if val == 1:
        return '1'
    elif val == 2:
        return '2'
    elif val >= 3:
        return '3'

def get_pro_pep(sLine, sDelimiter='\t', iFirstDelimiter=15):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    if iPos == -1:
        iPos = len(sLine)
    iProPep = int(sLine[iBegin:iPos])
    return iProPep

def get_pep_psm(sLine, sDelimiter='\t', iFirstDelimiter=14):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    if iPos == -1:
        iPos = len(sLine)
    iPepPsm = int(sLine[iBegin:iPos])
    return iPepPsm

# # get the score agreement count, this number is in the fourteenth column
def get_score_agreement(sLine, sDelimiter='\t', iFirstDelimiter=13):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    if iPos == -1:
        iPos = len(sLine)
    iScoreAgreement = int(sLine[iBegin:iPos])
    return iScoreAgreement

# # get the parent charge from the line, parent change is in the third column
def get_parent_charge(sLine, sDelimiter='\t', iFirstDelimiter=2):
    iPos = -1
    while iFirstDelimiter > 0:
        iPos = sLine.find(sDelimiter, iPos + 1)
        iFirstDelimiter -= 1
    iBegin = iPos + 1
    iPos = sLine.find(sDelimiter, iBegin)
    if iPos == -1:
        iPos = len(sLine)
    iCharge = int(sLine[iBegin:iPos])
    return iCharge

def categorize(input_file, output_folder):    
    
    with open(input_file, 'r') as f:
        dictGroup = {}
        # skip header
        sHeader = f.readline()
        if sHeader.count('\t') == 15:
            # read PSM
            while True:
                sLine = f.readline()
                if not sLine:
                    break
                sParentCharge = categorize_parent_charge(get_parent_charge(sLine))
                sScoreAgreement = categorize_score_agreement(get_score_agreement(sLine))
                sPepPsm = categorize_pep_psm(get_pep_psm(sLine))
                sProPep = categorize_pro_pep(get_pro_pep(sLine))
                sKey = sParentCharge + '_' + sScoreAgreement + '_' + sPepPsm + '_' + sProPep
                if not dictGroup.has_key(sKey):
                    fWriter = open(output_folder + os.path.splitext(os.path.basename(input_file))[0] + '_' + sKey + '.tab', 'w')
                    dictGroup[sKey] = fWriter
                    fWriter.write(sHeader)
                    fWriter.write(sLine)
                else:
                    fWriter = dictGroup[sKey]
                    fWriter.write(sLine)
        else:
            # read PSM
            while True:
                sLine = f.readline()
                if not sLine:
                    break
                sParentCharge = categorize_parent_charge(get_parent_charge(sLine))
                sScoreAgreement = categorize_score_agreement(get_score_agreement(sLine))
                sKey = sParentCharge + '_' + sScoreAgreement
                if not dictGroup.has_key(sKey):
                    fWriter = open(output_folder + os.path.splitext(os.path.basename(input_file))[0] + '_' + sKey + '.tab', 'w')
                    dictGroup[sKey] = fWriter
                    fWriter.write(sHeader)
                    fWriter.write(sLine)
                else:
                    fWriter = dictGroup[sKey]
                    fWriter.write(sLine)

    # close all the open files
    for _key, f in dictGroup.iteritems():
        f.close()

def clean_folder(output_folder):
    for the_file in os.listdir(output_folder):
        file_path = os.path.join(output_folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            # elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse options
    (input_file, output_folder) = parse_options(argv)

    # clean the folder
    clean_folder(output_folder)

    # categorize PSM table
    categorize(input_file, output_folder)

    #
    print "Done."

if __name__ == '__main__':
    sys.exit(main())
