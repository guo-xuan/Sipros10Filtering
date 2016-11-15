'''
Created on Jun 6, 2016

@author: xgo
'''

import sys
import getopt
from random import shuffle
import re
# # Import Sipros package modules
import sipros_post_module
import parseconfig

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_shuffle_database.py [options]

Inputs:
    database file
    sipros config file
    output filename

Options:
    -h/--help
    -v/--version
    -i/--input-file ./path    
    -c/--configuration-file SiprosConfig.cfg    # SIPROS config file
    -o/--output-file ./path

Outputs:
    output file
'''

# # Global variables
pep_iden_str = '[Peptide_Identification]'
cleave_after_str = 'Cleave_After_Residues'
cleave_before_str = 'Cleave_Before_Residues'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
aa_after_cleave_str = ""
aa_before_cleave_str = ""
decoy_str = ""

mark_beg_str = "(((("
mark_end_str = "))))"

# # check_file_exist
check_file_exist = sipros_post_module.check_file_exist
# # Import classes and definitions in the Sipros module
# # Class Usage
Usage = sipros_post_module.Usage


# # Parse options
def parse_options(argv):

    try:
        opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:",
                                    ["help",
                                     "version",
                                     "input-file",
                                     "configuration-file",
                                     "output-file"])

    # Error handling of options
    except getopt.error, msg:
        raise Usage(msg)

    # Default working dir and config file
    input_file = "./"
    config_file = "SiprosConfig.cfg"
    output_file = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            raise Usage(help_message)
        if option in ("-v", "-V", "--version"):
            print "sipros_peptides_filtering.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input-file"):
            input_file = value
        if option in ("-c", "--configuration-file"):
            config_file = value
        if option in ("-o", "--output-file"):
            output_file = value

    return (input_file, config_file, output_file)

def ShufflePep(seq_str):
    start = 0
    end = 1
    seq_new_str = ""
    pep_str = ""
    i = 0
    for i in range(0, len(seq_str) - 1):
        if seq_str[i] in aa_before_cleave_str and seq_str[i + 1] in aa_after_cleave_str:
            end = i
            pep_str = list(seq_str[start:end])
            shuffle(pep_str)
            seq_new_str += "".join(pep_str)
            seq_new_str += "".join(seq_str[i])
            start = end + 1
    # if seq_str[i+1] not in aa_after_cleave_str or seq_str[i] not in aa_before_cleave_str:
    if start <= len(seq_str) - 1:
        pep_str = list(seq_str[start:])
        shuffle(pep_str)
        seq_new_str += "".join(pep_str)
    return seq_new_str

def shuffle_pep_mutation(protein_id_str, protein_seq_str):
    # get the mutation information
    # # original position, mutation sequence, original index, new position
    all_mutations_list = []
    pos_list = []
    all_links_list = []
    if mark_beg_str in protein_id_str and mark_end_str in protein_id_str:
        split_list = protein_id_str[protein_id_str.find(mark_beg_str)+len(mark_beg_str):protein_id_str.find(mark_end_str)].split(',')
        mutation_list = re.split('([a-zA-Z]+)', split_list[0])
        # the last element is empty
        pos = 0
        idx = 0
        for i in range((len(mutation_list)-1)/2):
            pos += int(mutation_list[i*2])
            all_mutations_list.append([pos, mutation_list[i*2+1], idx])
            idx += 1
        # handle the link info
        if len(split_list) > 1:
            for i in range(1, len(split_list)):
                link_list = split_list[i].strip().split('-')
                all_links_list.append(link_list)
    # normal way to digest a protein            
    protein_seq_list = []
    cleavage_site_list = [0]
    for i in range(0, len(protein_seq_str) - 1):
        if protein_seq_str[i] in aa_before_cleave_str and protein_seq_str[i + 1] in aa_after_cleave_str:
            cleavage_site_list.append(i + 1);
    cleavage_site_list.append(len(protein_seq_str))
    # the last peptide needs full shuffle
    for i in range(1, len(cleavage_site_list) - 1):
        x = [j for j in range(cleavage_site_list[i - 1], cleavage_site_list[i] - 1)]
        shuffle(x)
        x.append(cleavage_site_list[i] - 1)
        pos_list.extend(x)
        pep_new = [protein_seq_str[j] for j in x]
        protein_seq_list.extend(pep_new)
    # handle the last peptide
    x = [i for i in range(cleavage_site_list[-2],  cleavage_site_list[-1])]
    shuffle(x)
    pos_list.extend(x)
    pep_new = [protein_seq_str[j] for j in x]
    protein_seq_list.extend(pep_new)
    protein_seq_str = ''.join(protein_seq_list)
    # if mutation is there, get the new position, new index
    if all_mutations_list:
        for a_list in all_mutations_list:
            a_list.append(pos_list.index(a_list[0]))
        new_id_list = sorted(range(len(all_mutations_list)),key=lambda x:all_mutations_list[x][3])
        for a_list in all_links_list:
            for i in range(len(a_list)):
                a_list[i] = str(new_id_list.index(int(a_list[i])))
        all_mutations_list = sorted(all_mutations_list, key = lambda a_list:a_list[3])
        mutation_new_list = [protein_id_str[0:protein_id_str.find(mark_beg_str)], mark_beg_str]
        pos = all_mutations_list[0][3]
        mutation_new_list.append(str(all_mutations_list[0][3]))
        mutation_new_list.append(str(all_mutations_list[0][1]))
        for i in range(1, len(all_mutations_list)):
            mutation_new_list.append(str(all_mutations_list[i][3] - pos))
            mutation_new_list.append(str(all_mutations_list[i][1]))
            pos = all_mutations_list[i][3]
        if all_links_list:
            for a_list in all_links_list:
                mutation_new_list.append(',')
                mutation_new_list.append('-'.join(a_list))
        mutation_new_list.append(mark_end_str)
        protein_id_str = ''.join(mutation_new_list)
    
    return (protein_id_str, protein_seq_str)

def reverse_seq_mutation(inputFileName, outputFileName, config_dict):
    outputFile = open(outputFileName, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    inputFile = open(inputFileName)
    line_str = ""
    iCount = 0
    for line_str in inputFile:
        if line_str[0] == '>':
            if not id_str.startswith(decoy_str) and seq_str != "":
                seq_str = seq_str.rstrip()
                outputFile.write(id_str)
                outputFile.write(seq_str)
                # generate shuffled sequences
                (id_new_str, seq_new_str) = shuffle_pep_mutation(id_str[1:].rstrip(), seq_str)
                outputFile.write("\n>Sh1_")
                outputFile.write(id_new_str)
                outputFile.write("\n")
                outputFile.write(seq_new_str)
                outputFile.write("\n")
                (id_new_str, seq_new_str) = shuffle_pep_mutation(id_str[1:].rstrip(), seq_str)
                outputFile.write('>Sh2_')
                outputFile.write(id_new_str)
                outputFile.write("\n")
                outputFile.write(seq_new_str)
                outputFile.write("\n")
                iCount += 1
                if iCount % 1000 == 0:
                    print "processed %d reads" % iCount
            id_str = line_str
            seq_str = ""
        else:
            seq_str += line_str.strip()
    if not id_str.startswith(decoy_str) and seq_str != "":
        seq_str = seq_str.rstrip()
        outputFile.write(id_str)
        outputFile.write(seq_str)
        # generate shuffled sequences
        (id_new_str, seq_new_str) = shuffle_pep_mutation(id_str[1:].rstrip(), seq_str)
        outputFile.write("\n>Sh1_")
        outputFile.write(id_new_str)
        outputFile.write("\n")
        outputFile.write(seq_new_str)
        outputFile.write("\n")
        (id_new_str, seq_new_str) = shuffle_pep_mutation(id_str[1:].rstrip(), seq_str)
        outputFile.write('>Sh2_')
        outputFile.write(id_new_str)
        outputFile.write("\n")
        outputFile.write(seq_new_str)
        outputFile.write("\n")
    inputFile.close()
    outputFile.close()

def ReverseSeq(inputFileName, outputFileName, config_dict) :
    outputFile = open(outputFileName, "w")
    id_str = ""
    seq_str = ""
    seq_new_str = ""
    inputFile = open(inputFileName)
    line_str = ""
    for line_str in inputFile:
        if line_str[0] == '>':
            if not id_str.startswith(decoy_str) and seq_str != "":
                seq_str = seq_str.rstrip()
                seq_new_str = ShufflePep(seq_str)
                outputFile.write(id_str)
                outputFile.write(seq_str)
                outputFile.write("\n>Sh1_")
                outputFile.write(id_str[1:])
                outputFile.write(seq_new_str)
                outputFile.write("\n")
                seq_new_str = ShufflePep(seq_str)
                outputFile.write('>Sh2_')
                outputFile.write(id_str[1:])
                outputFile.write(seq_new_str)
                outputFile.write("\n")
            '''    
            else:
                if id_str.startswith(decoy_str) and seq_str != "":
                    outputFile.write(id_str)
                    outputFile.write(seq_str)
                    outputFile.write("\n")
            '''
            id_str = line_str
            seq_str = ""
        else:
            seq_str += line_str.strip()
    if not id_str.startswith(decoy_str) and seq_str != "":
        seq_str = seq_str.rstrip()
        seq_new_str = ShufflePep(seq_str)
        outputFile.write(id_str)
        outputFile.write(seq_str)
        outputFile.write("\n>Sh1_")
        outputFile.write(id_str[1:])
        outputFile.write(seq_new_str)
        outputFile.write("\n")
        seq_new_str = ShufflePep(seq_str)
        outputFile.write('>Sh2_')
        outputFile.write(id_str[1:])
        outputFile.write(seq_new_str)
        outputFile.write("\n")
    '''
    else:
        if id_str.startswith(decoy_str) and seq_str != "":
            outputFile.write(id_str)
            outputFile.write(seq_str)
            outputFile.write("\n")
    '''
    inputFile.close()
    outputFile.close()

# # Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}  # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    # Save all config values to dictionary
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)

    # valiables were defined in global
    # pep_iden_str      = '[Protein_Identification]'
    # cleave_after_str  = 'Decoy_Prefix'
    # cleave_before_str = 'FDR_Filtering'
    # FDR_threshold_str = 'FDR_Threshold'

    # only save protein_identification config info to config_dict
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + cleave_after_str):
            config_dict[cleave_after_str] = value
        elif key == (pep_iden_str + cleave_before_str):
            config_dict[cleave_before_str] = value
        elif key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        else:
            continue

    # return config dictionary
    return config_dict

def main(argv=None):

    # try to get arguments and error handling
    try:
        if argv is None:
            argv = sys.argv
        try:
            # parse options
            (input_file, config_filename, output_file) = parse_options(argv)
            # Call parse_config to open and read config file
            config_dict = parse_config(config_filename)
            global aa_after_cleave_str
            aa_after_cleave_str = config_dict[cleave_before_str]
            global aa_before_cleave_str
            aa_before_cleave_str = config_dict[cleave_after_str]
            global decoy_str
            decoy_str = ">" + config_dict[decoy_prefix_str]
            # random shuffle database
            # ReverseSeq(input_file, output_file, config_dict)
            # new way shuffle sequence
            reverse_seq_mutation(input_file, output_file, config_dict)
            print "Done."

        # Error handling
        except Usage, err:
            print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
            return 2

    # Error handling
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "for help use -h/--help"
        return 2

if __name__ == '__main__':
    sys.exit(main())
