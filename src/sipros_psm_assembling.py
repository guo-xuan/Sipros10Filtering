'''
Created on Jul 26, 2016

@author: xgo
'''

import sys, getopt, os
import sipros_post_module
from _collections import defaultdict
import sipros_peptides_assembling

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_psm_assembling.py [options]

Inputs:
    input directory containing
        pep.txt file(s)
        psm.txt file(s)
    output directory

Options:
    -h/--help
    -v/--version
    -i/--input-folder
    -o/--output-folder

Outputs:
    output PSM table
'''


# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:f:p:r",
                                    ["help",
                                     "version",
                                     "input-folder",
                                     "configuration",
                                     "output-folder",
                                     "psm_tab_file,"
                                     "psm_tab_folder"])

    # Default working dir and config file
    input_folder = ""
    output_folder = ""
    sConfig = ""
    psm_tab_folder = ""
    psm_tab_file = ""
    FirstRound = True
    
    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "sipros_psm_assembling.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input-folder"):
            input_folder = value
        if option in ("-o", "--output-folder"):
            output_folder = value
        if option in ("-c", "--configuration"):
            sConfig = value
        if option in ("-p", "--psm_tab_folder"):
            psm_tab_folder = value
        if option in ("-f", "--psm_tab_file"):
            psm_tab_file = value
        if option in ("-r"):
            FirstRound = False

    if input_folder == "" or output_folder == "":
        print help_message
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')
    input_folder = os.path.join(input_folder, '')
    psm_tab_folder = os.path.join(psm_tab_folder, '')
    
    return (input_folder, sConfig, output_folder, psm_tab_file, psm_tab_folder, FirstRound)

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'

import parseconfig

## Check file exist
check_file_exist = sipros_post_module.check_file_exist

# defaul value
decoy_prefix = 'Rev_'
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'

## Parse config file
def parse_config(config_filename):

    # Save all config values to dictionary
    all_config_dict = {}    # initialize dictionay
    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    all_config_dict = parseconfig.parseConfigKeyValues(config_filename)

    # only save protein_identification config info to config_dict
    config_dict[decoy_prefix_str] = decoy_prefix
    config_dict[min_peptide_per_protein_str] = min_peptide_per_protein
    config_dict[min_unique_peptide_per_protein_str] = min_unique_peptide_per_protein
    config_dict[remove_decoy_identification_str] = remove_decoy_identification
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + fasta_database_str):
            config_dict[fasta_database_str] = value
        elif key == (pro_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pro_iden_str + min_peptide_per_protein_str):
            config_dict[min_peptide_per_protein_str] = value
        elif key == (pro_iden_str + min_unique_peptide_per_protein_str):
            config_dict[min_unique_peptide_per_protein_str] = value
        elif key == (pro_iden_str + remove_decoy_identification_str):
            config_dict[remove_decoy_identification_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


get_file_list_with_ext = sipros_post_module.get_file_list_with_ext
get_base_out = sipros_post_module.get_base_out

protein_type = sipros_post_module.protein_type

def merge_psm_files(input_folder, output_folder):
    # Get base_out for output
    base_out_default = 'Sipros_searches'
    sipros_file_list = get_file_list_with_ext(input_folder, 'psm.txt')
    if len(sipros_file_list) == 1:
        return sipros_file_list[0]
    base_out = get_base_out(sipros_file_list, base_out_default, input_folder)
    base_out_filename = base_out.split('/')[-1]
    base_out = output_folder + base_out_filename

    iNum = [0, 0, 0]

    with open(base_out + ".psm.txt", 'w') as fw:
        # for psm out
        psm_out_list = ['Filename',  # 0
                    'ScanNumber',  # 1
                    'ParentCharge',  # 2
                    'MeasuredParentMass',  # 3
                    'CalculatedParentMass',  # 4
                    'MassErrorDa',  # 5 CalculatedParentMass - MeasuredParentMass
                    'MassErrorPPM',  # 6 MassErrorDa / CalculatedParentMass
                    'ScanType',  # 7
                    'SearchName',  # 8
                    'ScoringFunction',  # 9
                    'Score',  # 10
                    'DeltaZ',  # 11 the difference score between the rank 1 and 2
                    'DeltaP',  # 12
                    'IdentifiedPeptide',  # 13
                    'OriginalPeptide',  # 14
                    'ProteinNames',  # 15
                    'ProteinCount',  # 16
                    'TargetMatch']  # 17
        fw.write('\t'.join(psm_out_list) + '\n')
        for input_file in sipros_file_list:
            with open(input_file, 'r') as fr:
                for sLine in fr:
                    fw.write(sLine)
                    # debug
                    asWords = sLine.split('\t')
                    iType = protein_type(asWords[15])
                    iNum[iType - 1] += 1
    print iNum
    return base_out + ".psm.txt"

class Peptide:
    
    def __init__(self):
        self.IdentifiedPeptide = ''
        self.ParentCharge = ''
        self.OriginalPeptide = ''
        self.ProteinNames = []
        self.ProteinCount = 0
        self.SpectralCount = 0
        self.BestScore = 0.0
        self.PSMs = []
        self.ScanType = []
        self.SearchName = []
        
    def add(self, asWords):
        for sProtein in self.splitWords(asWords[3]):
            if sProtein not in self.ProteinNames:
                self.ProteinNames.append(sProtein)
        self.SpectralCount += int(asWords[6])
        if self.BestScore < float(asWords[7]):
            self.BestScore = float(asWords[7])
        self.PSMs.append(asWords[8][1:-1])
        self.ScanType.append(asWords[9][1:-1])
        self.SearchName.append(asWords[10][1:-1])
        
        
    def set(self, asWords):
        self.IdentifiedPeptide = asWords[0]
        self.ParentCharge = asWords[1]
        self.OriginalPeptide = asWords[2]
        for sProtein in self.splitWords(asWords[3]):
            if sProtein not in self.ProteinNames:
                self.ProteinNames.append(sProtein)
        self.SpectralCount = int(asWords[6])
        if self.BestScore < float(asWords[7]):
            self.BestScore = float(asWords[7])
        self.PSMs.append(asWords[8][1:-1])
        self.ScanType.append(asWords[9][1:-1])
        self.SearchName.append(asWords[10][1:-1])
    
    def check_decoy_match(self, decoy_prefix):
        for sProtein in self.ProteinNames:
            if not sProtein.startswith(decoy_prefix):
                self.TargetMatch = 'T'
                return 
        self.TargetMatch = 'F'
        
    def splitWords(self, sWord):
        sWord = sWord.replace('{','')
        sWord = sWord.replace('}','')
        asWord = sWord.split(',')
        return asWord
        
    def __repr__(self):
        l = [self.IdentifiedPeptide,
             self.ParentCharge,
             self.OriginalPeptide,
             ('{'+','.join(self.ProteinNames)+'}'),
             str(len(self.ProteinNames)),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{'+','.join(self.PSMs)+'}'),
             ('{'+','.join(self.ScanType)+'}'),
             ('{'+','.join(self.SearchName)+'}')]
        
        return '\t'.join(l) 

def merge_pep_files(input_folder, decoy_prefix, output_folder):
    # Get base_out for output
    base_out_default = 'Sipros_searches'
    sipros_file_list = get_file_list_with_ext(input_folder, 'pep.txt')
    if len(sipros_file_list) == 1:
        return sipros_file_list[0]
    base_out = get_base_out(sipros_file_list, base_out_default, input_folder)
    base_out_filename = base_out.split('/')[-1]
    base_out = output_folder + base_out_filename
    
    # pep_sub_dict for preparing pep_out
    pep_sub_dict = defaultdict(list)    # initialize dict of list
    
    for input_file in sipros_file_list:
        with open(input_file, 'r') as fr:
            for sLine in fr:
                asWords = sLine.split('\t')
                pep_ID = asWords[0] + '_+_' + asWords[1]
                if pep_ID in pep_sub_dict:
                    pep_sub_dict[pep_ID].add(asWords)
                else:
                    oPeptide = Peptide()
                    oPeptide.set(asWords)
                    pep_sub_dict[pep_ID] = oPeptide
    
    with open(base_out + ".pep.txt", 'w') as fw:
        # for pep out
        pep_out_list = ['IdentifiedPeptide',    #0
                    'ParentCharge',         #1
                    'OriginalPeptide',      #2
                    'ProteinNames',         #3
                    'ProteinCount',         #4
                    'TargetMatch',          #5
                    'SpectralCount',        #6 number of PSMs matched to this peptide
                    'BestScore',            #7 the highest score of those PSMs
                    'PSMs',                 #8 a list of PSMs matched to this peptide. Use{Filename[ScanNumber],Filename[ScanNumber]} format
                    'ScanType',             #9
                    'SearchName']           #10
        fw.write('\t'.join(pep_out_list) + '\n')
        for _pep_id, oPeptide in pep_sub_dict.iteritems():
            oPeptide.check_decoy_match(decoy_prefix)
            fw.write(repr(oPeptide))
            fw.write('\n')
        
        
    return base_out + ".pep.txt"

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'

## check pep and psm files pair set, and save run#
def get_run_num(pep_file_list, psm_file_list):

    # dictionary of Run# for each pep_file
    run_num_dict = {}
    psm_run_num_dict = {}

    # read multiple pep files
    for pep_file_idx, pep_file in enumerate(pep_file_list):
    
        # check file exist and open
        check_file_exist(pep_file)

        # If base common prefix ends with '.pep.txt', then remove '.pep.txt'
        base_pep_file = pep_file.replace(pep_file_ext, "")

        # make a psm filename using pep filename
        psm_file = base_pep_file + psm_file_ext

        # check psm file exist
        check_file_exist(psm_file)

        # for Run#
        run_tag = 'Run' + str(pep_file_idx + 1)

        # run_num_dict
        run_num_dict[pep_file] = run_tag
        psm_run_num_dict[psm_file] = run_tag

    return (run_num_dict, psm_run_num_dict)

read_run_files = sipros_peptides_assembling.read_run_files
greedy_alg = sipros_peptides_assembling.greedy_alg
check_decoy_match = sipros_peptides_assembling.check_decoy_match
divide = sipros_post_module.divide
PrettyFloat = sipros_post_module.PrettyFloat
set_float_digit = sipros_post_module.set_float_digit
PepOutFields = sipros_post_module.PepOutFields
report_output = sipros_peptides_assembling.report_output
read_fasta_file = sipros_peptides_assembling.read_fasta_file

## Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile
import csv

# # need pep.txt file and pro.txt file, and big PSM table file
# # Spectral Count for original/identified peptide
# # Unique peptide counts and total peptide counts for protein
def feature_update(run_num_dict, pro_file, psm_tab, psm_tab_new, FirstRound):
    
    # save the pep file and pro file data to the defaultdict
    id_pep_data_dict = {}
    or_pep_data_dict = {}
    pro_data_dict = defaultdict(list)
    
    # key = pep_file, val = run_num , sorted by Run# index
    for pep_file, _run_num in sorted(run_num_dict.items(), key=lambda x: x[1][-1]):
        f = open(pep_file, 'rb')
        # read line with csv
        pep_reader = csv.reader(CommentedFile(f),
                                   delimiter='\t')
        # skip header
        _headline = pep_reader.next()
        
        # get data
        for pep_line in pep_reader:
            # adapt class PepOutFields
            pep_obj = PepOutFields._make(pep_line)
            identified_peptide = pep_obj.IdentifiedPeptide
            identified_peptide = identified_peptide.strip()
            ParentCharge = pep_obj.ParentCharge
            identified_peptide = identified_peptide +'_' + ParentCharge
            original_peptide = pep_obj.OriginalPeptide
            original_peptide = original_peptide.strip()
            
            iSpectralCount = int(pep_obj.SpectralCount.strip())
            
            if identified_peptide in id_pep_data_dict:
                id_pep_data_dict[identified_peptide] += iSpectralCount
            else:
                id_pep_data_dict[identified_peptide] = iSpectralCount
            
            if original_peptide in or_pep_data_dict:
                or_pep_data_dict[original_peptide] += iSpectralCount
            else:
                or_pep_data_dict[original_peptide] = iSpectralCount
        f.close()
    
    f = open(pro_file, 'rb')
    pro_reader = csv.reader(CommentedFile(f),
                                   delimiter='\t')
    # skip header
    headline = pro_reader.next()
    iNumRuns = (len(headline) - 2)/6
    iUniquePeptideCounts = 0
    iTotalPeptideCounts = 0
    for asWords in pro_reader:
        identified_protain = asWords[0]
        iUniquePeptideCounts = 0
        iTotalPeptideCounts = 0
        for i in range(iNumRuns):
            iUniquePeptideCounts += int(asWords[i*6+1])
            iTotalPeptideCounts += int(asWords[i*6+2])
        pro_data_dict[identified_protain] = [iUniquePeptideCounts, iTotalPeptideCounts]
    
    f.close()
    
    fr = open(psm_tab, 'rb')
    psm_reader = csv.reader(CommentedFile(fr),
                                   delimiter='\t')
    
    # skip header
    headline = psm_reader.next()
    
    with open(psm_tab_new, 'w') as f:
        # header
        f.write('FileName\t')
        f.write('ScanNumber\t')
        f.write('ParentCharge\t')
        f.write('MeasuredParentMass\t')
        f.write('ScanType\t')
        f.write('SearchName\t')
        f.write('IdentifiedPeptide\t')
        f.write('OriginalPeptide\t')
        f.write('CalculatedParentMass\t')
        f.write('MVH\t')
        f.write('Xcorr\t')
        f.write('WDP\t')
        f.write('ProteinNames\t')
        f.write('ScoreAgreement\t')
        f.write('pep_psm\t')
        f.write('pro_pep')
        f.write('\n')
        
        if FirstRound:
            for asWords in psm_reader:
                
                original_peptide = asWords[7]
                if original_peptide in or_pep_data_dict:
                    iNumPepPsm = or_pep_data_dict[original_peptide]
                    iNumPepPsm -= 1
                    if iNumPepPsm < 0:
                        print 'Error'
                        exit(1)
                else:
                    iNumPepPsm = 0
                '''
                identified_peptide = asWords[6]
                ParentCharge = asWords[2]
                identified_peptide = identified_peptide +'_' + ParentCharge
                if identified_peptide in id_pep_data_dict:
                    iNumPepPsm = id_pep_data_dict[identified_peptide]
                    iNumPepPsm -= 1
                    if iNumPepPsm < 0:
                        print 'Error'
                        exit(1)
                else:
                    iNumPepPsm = 0
                '''
                sProteins = asWords[12]
                asProteins = (sProteins[1:-1]).split()
                iNumTotalPep = 0
                iNumUniqPep = 0
                for sProtein in asProteins:
                    if sProtein in pro_data_dict:
                        l = pro_data_dict[sProtein]
                        iNumTotalPep += l[1]
                        iNumUniqPep += l[0]
                f.write('\t'.join(asWords))
                f.write('\t')
                f.write(str(iNumPepPsm))
                f.write('\t')
                if iNumUniqPep > 1:
                    f.write('3')
                elif iNumTotalPep > 1:
                    f.write('2')
                else:
                    f.write('1')
                f.write('\n')
        else:
            for asWords in psm_reader:
                original_peptide = asWords[7]
                if original_peptide in or_pep_data_dict:
                    iNumPepPsm = or_pep_data_dict[original_peptide]
                else:
                    iNumPepPsm = 0
                sProteins = asWords[12]
                asProteins = (sProteins[1:-1]).split()
                iNumTotalPep = 0
                iNumUniqPep = 0
                for sProtein in asProteins:
                    if sProtein in pro_data_dict:
                        l = pro_data_dict[sProtein]
                        iNumTotalPep += l[1]
                        iNumUniqPep += l[0]
                asWords = asWords[:-2]
                f.write('\t'.join(asWords))
                f.write('\t')
                f.write(str(iNumPepPsm))
                f.write('\t')
                if iNumUniqPep > 1:
                    f.write('3')
                elif iNumTotalPep > 1:
                    f.write('2')
                else:
                    f.write('1')
                f.write('\n')
    fr.close()
    
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
    
    # Parse options and get config file
    sys.stderr.write('[Step 1] Parse options and read fasta file: Running -> ')
    # Call parse_config to open and read config file
    # parse options
    (input_folder, config_filename, output_folder, psm_tab_file, psm_tab_folder, FirstRound) = parse_options(argv)
    
    # clean the folder
    clean_folder(output_folder)
    
    config_dict = parse_config(config_filename)
    merge_psm_files(input_folder, output_folder)
    merge_pep_files(input_folder, config_dict[decoy_prefix_str], output_folder)
    # Read fasta file and retrieve protein ID and description
    fasta_ID_dict = read_fasta_file(output_folder, config_dict)
    
    # Get .pep.txt output file(s) in working directory
    pep_file_list = get_file_list_with_ext(output_folder, pep_file_ext)
    # Get .psm.txt output file(s) in working directory
    psm_file_list = get_file_list_with_ext(output_folder, psm_file_ext)
    # check pep and psm files pair set, and save run#
    (run_num_dict, psm_run_num_dict) = get_run_num(pep_file_list, psm_file_list)
    # Get base_out for output
    base_out_default = 'Sipros_searches'
    base_out = get_base_out(pep_file_list, base_out_default, output_folder)
    sys.stderr.write('Done!\n')
    
    # Read and load pep and psm files
    sys.stderr.write('[Step 2] Load %s file(s): Running -> ' % (pep_file_ext))
    (pep_data_dict, psm_data_dict, pro_pep_dict, pep_pro_dict) = read_run_files(run_num_dict)
    sys.stderr.write('Done!\n')

    # Merge indistinguishable proteins that have an identical set of peptides
    sys.stderr.write('[Step 3] Merge indistinguishable proteins: Running -> Done!\n')

    # extract proteins that have >2 peptides and at least one of those is unique
    # then iteratively extract a protein at a time that covers the most peptides
    sys.stderr.write('[Step 4] Greedy algorithm for identifying a list of proteins: Running -> ')
    (pro_greedy_list) = greedy_alg(config_dict, pro_pep_dict, pep_pro_dict)
    sys.stderr.write('Done!\n')
    
    # Report output
    sys.stderr.write('[Step 5] Report output: Running -> ')

    # Report output files
    report_output(config_dict,
                  run_num_dict,
                  psm_run_num_dict,
                  pep_data_dict,
                  psm_data_dict,
                  pro_pep_dict,
                  pep_pro_dict,
                  pro_greedy_list,
                  fasta_ID_dict,
                  base_out)
    sys.stderr.write('Done!\n')
    
    # Feature Update
    sys.stderr.write('[Step 6] Feature update: Running -> ')
    base_out_filename = psm_tab_file.split('/')[-1]
    psm_tab_new = psm_tab_folder + base_out_filename
    feature_update(run_num_dict, base_out+'.pro.txt', psm_tab_file, psm_tab_new, FirstRound)
    
    sys.stderr.write('Done!\n')
    
    


if __name__ == '__main__':
    sys.exit(main())
