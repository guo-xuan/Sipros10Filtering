'''
Created on Jul 21, 2016

@author: Xuan Guo
'''

import getopt, sys, os, math
import sipros_post_module
import parseconfig
from _collections import defaultdict
import numpy as np
from scipy.stats import rankdata

bDebug = True
ISCOREAGREE = 1

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_psm_filtering.py [options]

Inputs:
    tab file
    sipros configuration file
    output directory

Options:
    -h/--help
    -v/--version
    -i/--input-file
    -c/--configuration
    -o/--output-folder

Outputs:
    output PSM table
'''


# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:",
                                    ["help",
                                     "version",
                                     "input-file",
                                     "configuration",
                                     "output-folder"])

    # Default working dir and config file
    input_file = ""
    sConfig = ""
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
        if option in ("-c", "--configuration"):
            sConfig = value
        if option in ("-o", "--output-folder"):
            output_folder = value

    if input_file == "" or sConfig == "" or output_folder == "":
        print help_message
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, sConfig, output_folder)

# # check_file_exist
check_file_exist = sipros_post_module.check_file_exist

pep_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
FDR_filtering_str = 'FDR_Filtering'
FDR_threshold_str = 'FDR_Threshold'

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
    # decoy_prefix_str  = 'Decoy_Prefix'
    # FDR_filtering_str = 'FDR_Filtering'
    # FDR_threshold_str = 'FDR_Threshold'

    # only save protein_identification config info to config_dict
    for key, value in all_config_dict.items():
        if key == (pep_iden_str + decoy_prefix_str):
            config_dict[decoy_prefix_str] = value
        elif key == (pep_iden_str + FDR_filtering_str):
            config_dict[FDR_filtering_str] = value
        elif key == (pep_iden_str + FDR_threshold_str):
            config_dict[FDR_threshold_str] = value
        else:
            continue

    # return config dictionary
    return config_dict

protein_type = sipros_post_module.protein_type

class PSM:

    iNumScores = 3

    def __init__(self):
        self.lfScores = []
        self.bIsTarget = False
        self.lRanks = []

    # 9th MVH, 10th Xcorr, 11th Wdp, 12th Proteins
    # FileName    ScanNumber    ParentCharge    MeasuredParentMass    ScanType
    # SearchName    IdentifiedPeptide    OriginalPeptide    CalculatedParentMass
    # MVH    Xcorr    WDP    ProteinNames    ScoreAgreement
    def readInfo(self, sLine="", decoy_prefix="Rev_"):
        asWords = sLine.strip().split('\t')
        self.FileName = asWords[0]
        self.ScanNumber = asWords[1]
        self.ParentCharge = asWords[2]
        self.MeasuredParentMass = asWords[3]
        self.ScanType = asWords[4]
        self.SearchName = asWords[5]
        self.IdentifiedPeptide = asWords[6]
        self.OriginalPeptide = asWords[7]
        self.CalculatedParentMass = asWords[8]
        self.lProtein = []
        for i in range(9, 12):
            self.lfScores.append(float(asWords[i]))
        self.check_decoy_match(asWords[12], decoy_prefix)
        self.fRankProduct = 0.0
        if len(asWords) >= 15:
            self.fRankProduct = float(asWords[14])

    def check_decoy_match(self, protein_sequence, decoy_prefix):
        sProteins = protein_sequence.replace('{', '')
        sProteins = sProteins.replace('}', '')
        asProteins = sProteins.split(',')
        self.lProtein.extend(asProteins)
        for sProtein in asProteins:
            if not sProtein.startswith(decoy_prefix):
                self.bIsTarget = True
                return
        self.bIsTarget = False

    def __repr__(self):
        if self.bIsTarget:
            sT = 'T'
        else:
            sT = 'F'
        l = [self.FileName,
             self.ScanNumber,
             self.ParentCharge,
             self.MeasuredParentMass,
             self.CalculatedParentMass,
             str(float(self.MeasuredParentMass) - float(self.CalculatedParentMass)),
             'NA',
             self.ScanType,
             self.SearchName,
             'Sipros10',
             'NA',
             'NA',
             'NA',
             self.IdentifiedPeptide,
             self.OriginalPeptide,
             '{' + ','.join(self.lProtein) + '}',
             str(len(self.lProtein)),
             sT]

        return '\t'.join(l)

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
        self.PassThreshold = False

    def add(self, oPsm):
        for sProtein in oPsm.lProtein:
            if sProtein not in self.ProteinNames:
                self.ProteinNames.append(sProtein)
        self.PSMs.append(oPsm.FileName + '[' + oPsm.ScanNumber + ']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if self.BestScore < oPsm.lfScores[-1]:
            self.BestScore = oPsm.lfScores[-1]
        self.SpectralCount += 1

    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide

    def check_decoy_match(self, decoy_prefix):
        for sProtein in self.ProteinNames:
            if not sProtein.startswith(decoy_prefix):
                self.TargetMatch = 'T'
                return
        self.TargetMatch = 'F'

    def __repr__(self):
        l = [self.IdentifiedPeptide,
             self.ParentCharge,
             self.OriginalPeptide,
             ('{' + ','.join(self.ProteinNames) + '}'),
             str(len(self.ProteinNames)),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{' + ','.join(self.PSMs) + '}'),
             ('{' + ','.join(self.ScanType) + '}'),
             ('{' + ','.join(self.SearchName) + '}')]

        return '\t'.join(l)

def read_tab_file(input_file, decoy_prefix):
    lPsm = []

    # pep_sub_dict for preparing pep_out
    pep_sub_dict = defaultdict(list)  # initialize dict of list

    with open(input_file, 'r') as f:
        # skip the header
        _sLine = f.readline()
        # read the PSM
        for sLine in f:
            oPsm = PSM()
            oPsm.readInfo(sLine, decoy_prefix)
            lPsm.append(oPsm)
            # pep ID is unique with IdentifiedPeptide and ParentCharge
            pep_ID = oPsm.IdentifiedPeptide + '_+_' + oPsm.ParentCharge
            if pep_ID in pep_sub_dict:
                pep_sub_dict[pep_ID].add(oPsm)
            else:
                oPeptide = Peptide()
                oPeptide.set(oPsm)
                oPeptide.add(oPsm)
                pep_sub_dict[pep_ID] = oPeptide

    return (lPsm, pep_sub_dict)

def get_max_min_avg(lPsm, iScoreIndex):

    if len(lPsm) == 0:
        return (0.0, 0.0, 0.0)

    max_score = lPsm[0].lfScores[iScoreIndex]
    min_score = lPsm[0].lfScores[iScoreIndex]
    avg_score = 0.0

    for oPsm in lPsm:
        fScore = oPsm.lfScores[iScoreIndex]
        if fScore > max_score:
            max_score = fScore
        if fScore < min_score:
            min_score = fScore
        avg_score += fScore
    avg_score /= len(lPsm)
    return (max_score, min_score, avg_score)

def get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex):
    F_num = 0.0
    T_num = 0.0
    for oPsm in lPsm:
        fScore = oPsm.lfScores[iScoreIndex]
        if fScore >= cutoff_score:
            if oPsm.bIsTarget:
                T_num += 1
            else:
                F_num += 1
    return (F_num, T_num)

FDR_parameter = 1.0

# # Division error handling
divide = sipros_post_module.divide

# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(FP) + float(TP)
    FDR_accept = True

    if  FDR_denominator == 0:
        FDR_value = 1.0
        FDR_accept = False
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)
        FDR_accept = True

    return (FDR_accept, float(FDR_value))

# # list_to_bracket
frange = sipros_post_module.frange


def get_cutoff_data_np(FDR_threshold_1, FDR_threshold_2, lPsm, iScoreIndex, F_list, T_list):

    # save cutoff score and other data to cutoff_data
    cutoff_data = []

    # get info of lists
    (max_score, min_score, avg_score) = get_max_min_avg(lPsm, iScoreIndex)
    step_size = float(avg_score / 10000)

    # initialize
    final_cutoff_score_1 = 0.0
    _final_accept_1 = False
    prev_TF_num_1 = 0
    final_cutoff_score_2 = 0.0
    prev_TF_num_2 = 0
    _final_accept_2 = False

    # denominator is zero or not
    FDR_accept = False

    F_num_best = 0
    T_num_best = 0
    TF_num_best = 0

    # it max_score=min_score=0.0
    if (max_score == 0.0) and (min_score == 0.0):
        # cutoff_data list
        cutoff_data[min_score, max_score]
    else:
        # get list from max to min with decrement stepsize
        cutoff_score_range = frange(max_score, min_score, -step_size)
        # for loop of cutoff_score_range
        for cutoff_score in cutoff_score_range:
            F_num = (F_list >= cutoff_score).sum()
            T_num = (T_list >= cutoff_score).sum()
            # (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
            TF_num = F_num + T_num
            (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)

            # update final values if conditions satisfies
            # 1) FDR_accept is True
            # 2) FDR_value should be less than or equal to FDR_threshold_1
            # 3) TF_num is greater than to previous TF_num
            if (FDR_accept is True) and (FDR_value <= FDR_threshold_1) and (TF_num > prev_TF_num_1) :
                final_cutoff_score_1 = cutoff_score
                _final_accept_1 = FDR_accept
                # previous TF_num
                prev_TF_num_1 = TF_num

            if (FDR_accept is True) and (FDR_value <= FDR_threshold_2) and (TF_num > prev_TF_num_2) :
                final_cutoff_score_2 = cutoff_score
                _final_accept_2 = FDR_accept
                # previous TF_num
                prev_TF_num_2 = TF_num
                F_num_best = F_num
                T_num_best = T_num
                TF_num_best = TF_num

        # cutoff_data list
        # score 1 should less than score 2
        cutoff_data = [final_cutoff_score_1, final_cutoff_score_2]
        sScore = 'Score'
        if iScoreIndex == 0:
            sScore = 'MVH'
        elif iScoreIndex == 1:
            sScore = 'Xcorr'
        elif iScoreIndex == 2:
            sScore = 'WDP'
        print "%s:\t%f\t%d\t%d\t%d" % (sScore, final_cutoff_score_2, T_num_best, F_num_best, TF_num_best)

    return cutoff_data
# # Find cut-off score (range) and other data using given FDR threshold and lists of TP and FT
def get_cutoff_data(FDR_threshold_1, FDR_threshold_2, lPsm, iScoreIndex):

    # save cutoff score and other data to cutoff_data
    cutoff_data = []

    # get info of lists
    (max_score, min_score, avg_score) = get_max_min_avg(lPsm, iScoreIndex)
    step_size = float(avg_score / 10000)

    # initialize
    final_cutoff_score_1 = 0.0
    _final_accept_1 = False
    prev_TF_num_1 = 0
    final_cutoff_score_2 = 0.0
    prev_TF_num_2 = 0
    _final_accept_2 = False

    # denominator is zero or not
    FDR_accept = False

    # it max_score=min_score=0.0
    if (max_score == 0.0) and (min_score == 0.0):
        # cutoff_data list
        cutoff_data[min_score, max_score]
    else:
        # get list from max to min with decrement stepsize
        cutoff_score_range = frange(max_score, min_score, -step_size)

        # for loop of cutoff_score_range
        for cutoff_score in cutoff_score_range:
            (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
            TF_num = F_num + T_num
            (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)

            # update final values if conditions satisfies
            # 1) FDR_accept is True
            # 2) FDR_value should be less than or equal to FDR_threshold_1
            # 3) TF_num is greater than to previous TF_num
            if (FDR_accept is True) and (FDR_value <= FDR_threshold_1) and (TF_num > prev_TF_num_1) :
                final_cutoff_score_1 = cutoff_score
                _final_accept_1 = FDR_accept
                # previous TF_num
                prev_TF_num_1 = TF_num

            if (FDR_accept is True) and (FDR_value <= FDR_threshold_2) and (TF_num > prev_TF_num_2) :
                final_cutoff_score_2 = cutoff_score
                _final_accept_2 = FDR_accept
                # previous TF_num
                prev_TF_num_2 = TF_num

        # cutoff_data list
        # score 1 should less than score 2
        cutoff_data = [final_cutoff_score_1, final_cutoff_score_2]

    return cutoff_data

# # get the false and true hit given multiple scores, as long as two scores pass the threshold
def get_target_decoy_hit_multi_scores_np(F_list_np, T_list_np, F_scores_np, T_scores_np):
    F_mat = (F_list_np >= F_scores_np)
    T_mat = (T_list_np >= T_scores_np)
    iNumScores = T_mat.shape[1]
    aFactor = [1] * iNumScores
    Factor_np = np.array(aFactor)
    F_list = np.dot(F_mat, Factor_np)
    T_list = np.dot(T_mat, Factor_np)
    F_num = (F_list > ISCOREAGREE).sum()
    T_num = (T_list > ISCOREAGREE).sum()
    return (F_num, T_num)

# # get the false and true hit given multiple scores, as long as two scores pass the threshold
def get_target_decoy_hit_multi_scores(scores, lPsm):
    F_num = 0.0
    T_num = 0.0
    iNumScores = lPsm[0].iNumScores
    for oPsm in lPsm:
        iPassCount = 0
        for i in range(iNumScores):
            if oPsm.lfScores[i] >= scores[i]:
                iPassCount += 1
        if iPassCount >= 2:
            if oPsm.bIsTarget:
                T_num += 1
            else:
                F_num += 1
    return (F_num, T_num)

def pass_or_not(final_cutoff_score, scores):
    iCount = 0
    for i in range(len(scores)):
        if scores[i] >= final_cutoff_score[i]:
            iCount += 1
    if iCount > ISCOREAGREE:
        return True
    else:
        return False

def decode(iNumber, iBase, iSize, lCodes):
    del lCodes[:]
    for _i in range(iSize):
        lCodes.append(iNumber % iBase)
        iNumber = int(iNumber / iBase)

def next_code(lCodes, iBase):
    iCarry = 0
    lCodes[0] += 1
    for i in range(len(lCodes)):
        lCodes[i] += iCarry
        iCarry = 0
        if lCodes[i] >= iBase:
            iCarry = 1
            lCodes[i] = 0
    return iCarry

# using the q value rank
def get_cutoff_q_rank_product(FDR_threshold, lPsm):
    iNumScores = lPsm[0].iNumScores
    # # the cutoff score range for each score
    lUpperLowerbound = []
    for i in range(iNumScores):
        F_list = []
        T_list = []
        for oPsm in lPsm:
            if oPsm.bIsTarget:
                T_list.append(oPsm.lfScores[i])
            else:
                F_list.append(oPsm.lfScores[i])
        F_list_np = np.array(F_list)
        T_list_np = np.array(T_list)
        if i == 0:
            print "Before Filtering:\t\t%d\t%d\t%d" % (T_list_np.size, F_list_np.size, T_list_np.size + F_list_np.size)
        lUpperLowerbound.append(get_cutoff_data_np(10 * FDR_threshold, FDR_threshold, lPsm, i, F_list_np, T_list_np))
    
    for i in range(iNumScores):
        newlist = sorted(lPsm, key=lambda x: x.lfScores[i], reverse=True)
        T_num = 0
        F_num = 0
        for j in range(len(newlist)):
            if newlist[j].bIsTarget:
                T_num += 1
            else:
                F_num += 1
            (_FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
            newlist[j].lRanks.append(FDR_value)
        fSmallestQ = 1
        for j in range(len(newlist) - 1, -1, -1):
            if fSmallestQ > newlist[j].lRanks[i]:
                fSmallestQ = newlist[j].lRanks[i]
            if newlist[j].lRanks[i] > fSmallestQ:
                newlist[j].lRanks[i] = fSmallestQ
    
    for oPsm in lPsm:
        #fTemp = oPsm.lRanks[0] * oPsm.lRanks[1] * oPsm.lRanks[2]
        #oPsm.fRankProduct = np.power(float(fTemp), 1.0 / 3.0)
        fTemp = (1 - ((1 - oPsm.lRanks[0]) * (1 - oPsm.lRanks[1]) * (1 - oPsm.lRanks[2])))
        oPsm.fRankProduct = fTemp
    
    final_cutoff_score = get_cutoff_rank_product2(FDR_threshold, lPsm)
    return final_cutoff_score
    
    

# using the global rank
def get_cutoff_rank_product2(FDR_threshold, lPsm):
    F_list = []
    T_list = []
    for oPsm in lPsm:
        if oPsm.bIsTarget:
            T_list.append(oPsm.fRankProduct)
        else:
            F_list.append(oPsm.fRankProduct)
    F_list = np.array(F_list)
    T_list = np.array(T_list)
    prev_TF_num = 0
    final_cutoff_score = np.amin(T_list, axis=0)
    _final_accept = False
    T_num_best = 0
    F_num_best = 0
    TF_num_best = 0
    for cutoff_score in T_list:
        F_num = (F_list <= cutoff_score).sum()
        T_num = (T_list <= cutoff_score).sum()
        # (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
        TF_num = F_num + T_num
        (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
        # update final values if conditions satisfies
        # 1) FDR_accept is True
        # 2) FDR_value should be less than or equal to FDR_threshold_1
        # 3) TF_num is greater than to previous TF_num
        if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
            final_cutoff_score = cutoff_score
            _final_accept = FDR_accept
            # previous TF_num
            prev_TF_num = TF_num
            T_num_best = T_num
            F_num_best = F_num
            TF_num_best = TF_num
            
    print "Q-value_cutoff:\t%f\t%d\t%d\t%d" % (final_cutoff_score, T_num_best, F_num_best, TF_num_best)
    return final_cutoff_score

# rank inside each category
def get_cutoff_rank_product(FDR_threshold, lPsm):
    # set the ranks for each scores
    iNumScores = lPsm[0].iNumScores
    lRanks = []
    for i in range(iNumScores):
        TF_list = []
        for oPsm in lPsm:
            TF_list.append((0 - oPsm.lfScores[i]))
        TF_list_np = np.array(TF_list)
        TF_list_rank = rankdata(TF_list_np, method='max')
        lRanks.append(TF_list_rank)
    Psm_rank_np = np.array(([1] * len(lPsm)))
    for i in range(iNumScores):
        Psm_rank_np = np.multiply(Psm_rank_np, lRanks[i])
    # find the cutoff rank
    F_list = []
    T_list = []
    for i in range(len(Psm_rank_np)):
        lPsm[i].fRankProduct = Psm_rank_np[i]
        if lPsm[i].bIsTarget:
            T_list.append(Psm_rank_np[i])
        else:
            F_list.append(Psm_rank_np[i])
    F_list = np.array(F_list)
    T_list = np.array(T_list)
    prev_TF_num = 0
    final_cutoff_score = np.amin(T_list, axis=0)
    _final_accept = False
    T_num_best = 0
    F_num_best = 0
    TF_num_best = 0
    for cutoff_score in T_list:
        F_num = (F_list <= cutoff_score).sum()
        T_num = (T_list <= cutoff_score).sum()
        # (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
        TF_num = F_num + T_num
        (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
        # update final values if conditions satisfies
        # 1) FDR_accept is True
        # 2) FDR_value should be less than or equal to FDR_threshold_1
        # 3) TF_num is greater than to previous TF_num
        if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
            final_cutoff_score = cutoff_score
            _final_accept = FDR_accept
            # previous TF_num
            prev_TF_num = TF_num
            T_num_best = T_num
            F_num_best = F_num
            TF_num_best = TF_num
            
    print str(final_cutoff_score)
    print "\t\t\t# Target:\t%d\t# Reverse:\t%d\t# Total:\t%d" % (T_num_best, F_num_best, TF_num_best)
    return final_cutoff_score

def get_cutoff_scores(FDR_threshold, lPsm, iStepSize):
    iNumScores = lPsm[0].iNumScores

    # # the cutoff score range for each score
    lUpperLowerbound = []
    for i in range(iNumScores):
        F_list = []
        T_list = []
        for oPsm in lPsm:
            if oPsm.bIsTarget:
                T_list.append(oPsm.lfScores[i])
            else:
                F_list.append(oPsm.lfScores[i])
        F_list_np = np.array(F_list)
        T_list_np = np.array(T_list)
        if i == 0:
            print "\t%d\t%d\t%d" % (T_list_np.size, F_list_np.size, T_list_np.size + F_list_np.size)
        lUpperLowerbound.append(get_cutoff_data_np(10 * FDR_threshold, FDR_threshold, lPsm, i, F_list_np, T_list_np))

    '''
    for i in range(iNumScores):
        lUpperLowerbound.append(get_cutoff_data(2*FDR_threshold, FDR_threshold, lPsm, i))
    '''

    # brute force search
    # initialize
    final_cutoff_score = [ 0.0 for i in range(iNumScores) ]
    final_accept = False
    prev_TF_num = 0
    F_num_best = 0
    T_num_best = 0
    TF_num_best = 0

    # try np array
    F_list = []
    T_list = []
    for oPsm in lPsm:
        if oPsm.bIsTarget:
            T_list.append(oPsm.lfScores)
        else:
            F_list.append(oPsm.lfScores)
    F_list_np = np.array(F_list)
    T_list_np = np.array(T_list)

    lMaxScoresTarget = np.amax(T_list_np, axis=0)
    lMaxScoresDecoy = np.amax(F_list_np, axis=0)
    lMaxScores = np.array([lMaxScoresTarget, lMaxScoresDecoy])
    lMaxScores = np.amax(lMaxScores, axis=0)

    lStepSize = []
    for i in range(iNumScores):
        lStepSize.append((lUpperLowerbound[i][1] - lUpperLowerbound[i][0]) / float(iStepSize - 1))
    # iTotalSearch = int(math.pow(iStepSize, iNumScores))
    lScores = []
    bReachUpperBound = False
    lPreCodes = [-1] * iNumScores
    lCodes = [0] * iNumScores
    lPreCodes = np.array(lPreCodes)
    lCodes = np.array(lCodes)
    iBase = iStepSize
    iTry = 0
    while (not final_accept) and (iTry < 2):
        iCarry = 0
        while iCarry == 0:
            if np.all((lCodes <= lPreCodes)):
                iCarry = next_code(lCodes, iBase)
                continue
            del lScores[:]
            for j in range(iNumScores):
                lScores.append(lUpperLowerbound[j][0] + lCodes[j] * lStepSize[j])
            F_scores = [lScores] * len(F_list)
            F_scores_np = np.array(F_scores)
            T_scores = [lScores] * len(T_list)
            T_scores_np = np.array(T_scores)
            (F_num, T_num) = get_target_decoy_hit_multi_scores_np(F_list_np, T_list_np, F_scores_np, T_scores_np)
            # (F_num, T_num) = get_target_decoy_hit_multi_scores(lScores, lPsm)
            TF_num = F_num + T_num
            (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
            if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
                final_cutoff_score = list(lScores)
                final_accept = True
                # previous TF_num
                prev_TF_num = TF_num
                F_num_best = F_num
                T_num_best = T_num
                TF_num_best = TF_num

            iCarry = next_code(lCodes, iBase)
        print '\t'.join(str(x) for x in final_cutoff_score)
        print "\t\t\t# Target:\t%d\t# Reverse:\t%d\t# Total:\t%d" % (T_num_best, F_num_best, TF_num_best)

        bReachUpperBound = True
        for i in range(iNumScores):
            if lScores[i] < lMaxScores[i]:
                bReachUpperBound = False
        lPreCodes = np.array(([iBase - 1] * iNumScores))
        iBase += iStepSize
        iTry += 1
    '''
    while (not final_accept) and (not bReachUpperBound):
        for i in range(iTotalSearch):
            
            decode(i, iStepSize, iNumScores, lCodes)
            del lScores[:]
            for j in range(iNumScores):
                lScores.append(lUpperLowerbound[j][0]+lCodes[j]*lStepSize[j] + lPreCodes[j] * lStepSize[j])
            F_scores = [lScores] * len(F_list)
            F_scores_np = np.array(F_scores)
            T_scores = [lScores] * len(T_list)
            T_scores_np = np.array(T_scores)
            (F_num, T_num) = get_target_decoy_hit_multi_scores_np(F_list_np, T_list_np, F_scores_np, T_scores_np)
            #(F_num, T_num) = get_target_decoy_hit_multi_scores(lScores, lPsm)
            TF_num = F_num + T_num
            (FDR_accept, FDR_value) = FDR_calculator(F_num, T_num)
        
            if (FDR_accept is True) and (FDR_value <= FDR_threshold) and (TF_num > prev_TF_num) :
                final_cutoff_score = list(lScores)
                final_accept = True
                # previous TF_num
                prev_TF_num = TF_num
                F_num_best = F_num
                T_num_best = T_num
                TF_num_best = TF_num
        
        print '\t'.join(str(x) for x in final_cutoff_score)
        print "\t\t\t# Target:\t%d\t# Reverse:\t%d\t# Total:\t%d" % (T_num_best, F_num_best, TF_num_best)
    
        bReachUpperBound = True
        for i in range(iNumScores):
            if lScores[i] < lMaxScores[i]:
                bReachUpperBound = False
        lPreCodes = lPreCodes + np.array(lCodes)
    '''
    if iTry >= 2 and all(v == 0 for v in final_cutoff_score):
        for i in range(iNumScores):
            final_cutoff_score[i] = lMaxScores[i]
    return final_cutoff_score

def report_output_rank_product(psm_tab_file,
                  lPsm,
                  pep_sub_dict,
                  final_cutoff_score,
                  config_dict):
    # filtering parameters
    decoy_prefix = str(config_dict[decoy_prefix_str])

    with open((os.path.splitext(psm_tab_file)[0] + '.psm.txt'), 'w') as f:
        for oPsm in lPsm:
            if final_cutoff_score >= oPsm.fRankProduct:
                f.write(repr(oPsm))
                f.write('\n')
                # pep ID is unique with IdentifiedPeptide and ParentCharge
                pep_ID = oPsm.IdentifiedPeptide + '_+_' + oPsm.ParentCharge
                if pep_ID in pep_sub_dict:
                    pep_sub_dict[pep_ID].PassThreshold = True

    with open((os.path.splitext(psm_tab_file)[0] + '.pep.txt'), 'w') as f:
        for _pep_id, oPeptide in pep_sub_dict.iteritems():
            oPeptide.check_decoy_match(decoy_prefix)
            if oPeptide.PassThreshold:
                f.write(repr(oPeptide))
                f.write('\n')

# # Report output files
def report_output(psm_tab_file,
                  lPsm,
                  pep_sub_dict,
                  final_cutoff_score,
                  config_dict):
    # filtering parameters
    decoy_prefix = str(config_dict[decoy_prefix_str])

    with open((os.path.splitext(psm_tab_file)[0] + '.psm.txt'), 'w') as f:
        for oPsm in lPsm:
            if pass_or_not(final_cutoff_score, oPsm.lfScores):
                f.write(repr(oPsm))
                f.write('\n')
                # pep ID is unique with IdentifiedPeptide and ParentCharge
                pep_ID = oPsm.IdentifiedPeptide + '_+_' + oPsm.ParentCharge
                if pep_ID in pep_sub_dict:
                    pep_sub_dict[pep_ID].PassThreshold = True

    with open((os.path.splitext(psm_tab_file)[0] + '.pep.txt'), 'w') as f:
        for _pep_id, oPeptide in pep_sub_dict.iteritems():
            oPeptide.check_decoy_match(decoy_prefix)
            if oPeptide.PassThreshold:
                f.write(repr(oPeptide))
                f.write('\n')

get_file_list_with_ext = sipros_post_module.get_file_list_with_ext

def parse_filename(input_file):
    filename = input_file.split('/')[-1]
    asWords = filename[:-4].split('_')
    l = []
    temp = 0
    for s in asWords:
        try:
            temp = int(s)
            if temp < 5:
                l.append(temp)
        except ValueError:
            pass
    if len(l)==2:
        print 'Charge\t%d\tAgreement\t%d' % (l[0], l[1])
    if len(l)==4:
        print 'Charge\t%d\tAgreement\t%d\tPep\t%d\tPro\t%d' % (l[0], l[1], l[2], l[3])
    

def main2(argv=None):

    if argv is None:
        argv = sys.argv

    # parse options
    (psm_tab_folder, sConfig, _output_folder) = parse_options(argv)

    # read the configuration file
    config_dict = parse_config(sConfig)

    psm_tab_folder = os.path.join(psm_tab_folder, '')
    lFileList = get_file_list_with_ext(psm_tab_folder, 'tab')

    for input_file in lFileList:
        # read the PSM
        #sys.stderr.write('Processing %s:\tRunning -> \n' % input_file)
        parse_filename(input_file)
        (lPsm, pep_sub_dict) = read_tab_file(input_file, config_dict[decoy_prefix_str])
        # find the enumerate range for each score
        #final_cutoff_score = get_cutoff_scores(float(config_dict[FDR_threshold_str]), lPsm, 10)
        #final_cutoff_score = get_cutoff_rank_product2(float(config_dict[FDR_threshold_str]), lPsm)
        final_cutoff_score = get_cutoff_q_rank_product(float(config_dict[FDR_threshold_str]), lPsm)
        # Report output
        #report_output(input_file, lPsm, pep_sub_dict, final_cutoff_score, config_dict)
        report_output_rank_product(input_file, lPsm, pep_sub_dict, final_cutoff_score, config_dict)
        sys.stdout.write('\n')
    
    sys.stdout.write('Done!\n')

def main(argv=None):

    if argv is None:
        argv = sys.argv

    # parse options
    (psm_tab_file, sConfig, _output_folder) = parse_options(argv)

    # read the configuration file
    config_dict = parse_config(sConfig)

    # read the PSM
    (lPsm, pep_sub_dict) = read_tab_file(psm_tab_file, config_dict[decoy_prefix_str])

    # find the enumerate range for each score
    final_cutoff_score = get_cutoff_scores(float(config_dict[FDR_threshold_str]), lPsm, 10)
    #final_cutoff_score = get_cutoff_rank_product(float(config_dict[FDR_threshold_str]), lPsm)
    # Report output
    report_output(psm_tab_file, lPsm, pep_sub_dict, final_cutoff_score, config_dict)
    #report_output_rank_product(psm_tab_file, lPsm, pep_sub_dict, final_cutoff_score, config_dict)

if __name__ == '__main__':
    sys.exit(main2())
