'''
Created on Sep 7, 2016

@author: xgo
'''

import getopt, sys, os
import numpy as np
import csv
import math
import re
from sklearn.tree import DecisionTreeClassifier
from collections import namedtuple
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
from sklearn import preprocessing
# from sklearn.neural_network import MLPClassifier
from subprocess import call
from multiprocessing import Process
from multiprocessing import Queue, cpu_count
import scipy.stats as spst

# # Import Sipros package modules
import sipros_post_module
import sipros_peptides_assembling
import parseconfig

# # Class for ignoring comments '#' in sipros file
CommentedFile = sipros_post_module.CommentedFile

#feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'NSI', 'NSM', 'NSP', 'NSPS', 'pep_psm', 'pro_pep']
feature_name_list = ['ParentCharge', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'MassDifferent', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3', 'NMC', 'NSI', 'NSM', 'NSP', 'NSPS', 'MassWindow']
feature_selection_list = [0, 1, 2, 3, 4, 5]

ptm_str = ['~', '!', '@', '>', '<', '%', '^', '&', '*', '(', ')', '/', '$']
# ptm_selection_list = [0]
ptm_selection_list = [0, 2, 3, 4]
for x in ptm_selection_list:
    feature_name_list.append(ptm_str[x])

rev_str = 'Sh1_' # 'Rev_'
shu_str = 'Sh2_' # 'Shu_'

# rev_str = 'Rev_'
# shu_str = 'Dec_'

mass_tolerance = 0.09

# # Class for PepOutFields object
class PsmFields(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'ScoreAgreement'])):  # 13

    def __init__(self):
        self.data = self
        
# # Class for PepOutFields object
class PsmFields2(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'ScoreAgreement', # 13
         'DeltaRP1',
         'DeltaRP2',
         'DeltaRP3',
         'DeltaRS1',
         'DeltaRS2',
         'DeltaRS3',
         'DiffRP1',
         'DiffRP2',
         'DiffRP3',
         'DiffRS1',
         'DiffRS2',
         'DiffRS3',
         'DiffNorRP1',
         'DiffNorRP2',
         'DiffNorRP3',
         'DiffNorRS1',
         'DiffNorRS2',
         'DiffNorRS3'])):  

    def __init__(self):
        self.data = self

        
# # Class for PepOutFields object
class PsmFields3(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'ScoreAgreement', # 13
         'DeltaRP1',
         'DeltaRP2',
         'DeltaRP3',
         'DeltaRS1',
         'DeltaRS2',
         'DeltaRS3',
         'DiffRP1',
         'DiffRP2',
         'DiffRP3',
         'DiffRS1',
         'DiffRS2',
         'DiffRS3',
         'DiffNorRP1',
         'DiffNorRP2',
         'DiffNorRP3',
         'DiffNorRS1',
         'DiffNorRS2',
         'DiffNorRS3',
         'pep_psm',
         'pro_pep'])):  

    def __init__(self):
        self.data = self

# # Class for PepOutFields object
class PsmFields4(namedtuple('PsmFields',
        ['FileName',  # 0
         'ScanNumber',  # 1
         'ParentCharge',  # 2
         'MeasuredParentMass',  # 3
         'ScanType',  # 4
         'SearchName',  # 5
         'IdentifiedPeptide',  # 6
         'OriginalPeptide',  # 7
         'CalculatedParentMass',  # 8
         'MVH',  # 9
         'Xcorr',  # 10
         'WDP',  # 11
         'ProteinNames',  # 12
         'ScoreAgreement', # 13
         'DeltaRP1',
         'DeltaRP2',
         'DeltaRP3',
         'DeltaRS1',
         'DeltaRS2',
         'DeltaRS3',
         'DiffRP1',
         'DiffRP2',
         'DiffRP3',
         'DiffRS1',
         'DiffRS2',
         'DiffRS3',
         'DiffNorRP1',
         'DiffNorRP2',
         'DiffNorRP3',
         'DiffNorRS1',
         'DiffNorRS2',
         'DiffNorRS3',
         'RetentionTime',
         'Rank',
         'DeltaP'])): # 33 ,
          # 'Rank'])): # 33 

    def __init__(self):
        self.data = self

LabelUnknown = 2
LabelPositive = 1
LabelNegative = 0
LabelFiltered = 3

LabelFwd = 1
LabelRev = 2
LabelShu = 3

shuffle_num_control = {'3_1': 3, 
                       '3_MX': 31,
                       '3_XW': 11,
                       '3_MW': 10,
                       '3_A': 414,
                       '2_1': 5,
                       '2_MX': 16,
                       '2_XW': 1,
                       '2_MW': 13,
                       '2_A': 671}

'''
shuffle_num_control = {'3_1': 3, 
                       '3_2': 37,
                       '3_3': 414,
                       '2_1': 5,
                       '2_2': 14,
                       '2_3': 671}
'''
feature_list_str = ['ParentCharge', 'MassDiff', 'MVH', 'Xcorr', 'WDP', 'ScoreAgreement', 'DeltaRP1', 'DeltaRP2', 'DeltaRP3', 'DeltaRS1', 'DeltaRS2', 'DeltaRS3', 'DiffRP1', 'DiffRP2', 'DiffRP3', 'DiffRS1', 'DiffRS2', 'DiffRS3', 'DiffNorRP1', 'DiffNorRP2', 'DiffNorRP3', 'DiffNorRS1', 'DiffNorRS2', 'DiffNorRS3']

class PSM:

    iNumScores = 3
    fProtonMass = 1.007276466812
    pattern = re.compile('[^\w\[\]]')

    def __init__(self, psm_field):
        if type(psm_field).__name__ == 'PsmFields3':
            self.pep_psm = int(psm_field.pep_psm)
            self.pro_pep = int(psm_field.pro_pep)
        else:
            self.pep_psm = -1
            self.pro_pep = -1
        self.FileName = psm_field.FileName
        self.bFileNameChanged = False
        self.ScanNumber = int(psm_field.ScanNumber)
        self.ParentCharge = int(psm_field.ParentCharge)
        self.ScanType = psm_field.ScanType
        self.SearchName = psm_field.SearchName
        self.lfScores = [float(psm_field.MVH), float(psm_field.Xcorr), float(psm_field.WDP)]
        self.ProteinNames = psm_field.ProteinNames
        self.ScoreAgreement = int(psm_field.ScoreAgreement)
        self.IdentifiedPeptide = psm_field.IdentifiedPeptide
        self.OriginalPeptide = psm_field.OriginalPeptide
        self.OriginalPeptide = PSM.pattern.sub('', self.IdentifiedPeptide)
        self.protein_list = []
        self.RealLabel = protein_type(self.ProteinNames, self.protein_list)
        self.lRanks = []
        self.TrainingLabel = LabelNegative  # 0: negative 1: positive 2: unknown
        self.fRankProduct = 0.0
        self.iInnerId = 0
        self.fPredictProbability = 0.0
        self.fMassDiff = 0.0
        self.MeasuredParentMass = float(psm_field.MeasuredParentMass)
        self.CalculatedParentMass = float(psm_field.CalculatedParentMass)
        self.iMassWindow = 0
        self.set_mass_diff(self.MeasuredParentMass, self.CalculatedParentMass)
        
        self.score_differential_list = []
        self.sRTime = '-1.000'
        self.fRtMeasured = 0.0
        self.fRtPredict = 0.0
        self.fRtPvalue = 0.0
        self.iLocalRank = 0
        self.DeltaP = 'NA'
        if type(psm_field).__name__ == 'PsmFields3':
            self.score_differential_list.extend(float(i) for i in psm_field[14:-2])
        elif type(psm_field).__name__ == 'PsmFields4':
            self.score_differential_list.extend(float(i) for i in psm_field[14:-1])
            self.sRTime = psm_field.RetentionTime
            self.fRtMeasured = float(self.sRTime)
            self.DeltaP = psm_field.DeltaP 
            # self.iLocalRank = int(psm_field.Rank)
        else:
            self.score_differential_list.extend(float(i) for i in psm_field[14:])
        
        
        self.NMC = 0
        self.NSI = 0
        self.NSM = 0
        self.NSP = 0 # unique peptide
        self.NSPS = 0 # shared peptide
        self.feature_list = []
        
        self.ML_feature = []
        self.FDR_feature = []
        self.fdr_product = 0
        
    def get_feature_list(self):
        del self.feature_list[:]
        if self.ParentCharge <= 2:
            self.feature_list.append(2)
        else:
            self.feature_list.append(3)
        self.feature_list.extend(self.lfScores)
        self.feature_list.append(self.ScoreAgreement)
        self.feature_list.append(self.fMassDiff)
        self.feature_list.extend(self.score_differential_list)
        self.feature_list.append(self.NMC)
        self.feature_list.append(self.NSI)
        self.feature_list.append(self.NSM)
        self.feature_list.append(self.NSP)
        self.feature_list.append(self.NSPS)
        
        if self.pep_psm != -1:
            self.feature_list.append(self.pep_psm)
            self.feature_list.append(self.pro_pep)
            
        self.feature_list.append(self.iMassWindow)
        
        for c in ptm_selection_list:
            self.feature_list.append(self.IdentifiedPeptide.count(ptm_str[c]))
        
    
    def set_protein_names(self):
        self.ProteinNames = '{' + ','.join(self.protein_list) + '}'

    def set_feature(self, feature_list):
        del feature_list[:]
        '''
        if self.ParentCharge <= 2:
            feature_list.append(2)
        else:
            feature_list.append(3)
        '''
        #feature_list.append(self.ParentCharge)
        feature_list.extend(self.lfScores)
        '''
        for x in self.lRanks:
            if x == 0:
                feature_list.append(math.log(0.000001))
            else:
                feature_list.append(math.log(x))
        '''
        #feature_list.append(self.lfScores[0])
        #feature_list.append(self.lfScores[1])
        #feature_list.append(self.lfScores[2])
        #feature_list.append(self.ScoreAgreement)
        #feature_list.append(self.fMassDiff)
        
    def set_mass_diff(self, measured_mass, calculated_mass):
        MassDiffOriginal = measured_mass - calculated_mass
        MassDiff = MassDiffOriginal
        for i in range(-1, 4):
            if abs(MassDiffOriginal - i*PSM.fProtonMass) < abs(MassDiff):
                MassDiff = MassDiffOriginal - i*PSM.fProtonMass
                self.iMassWindow = i
        self.fMassDiff = MassDiff
        
    def set_fdr_product(self):
        val = 1.0
        for x in self.FDR_feature:
            val *= (1.0 - x)
        val = 1.0 - pow(val, 1.0/float(len(self.FDR_feature)))
        self.fdr_product = val
        
    def clean_protein_name(self):
        self.ProteinNames = ""
        l = []
        for sProtein in self.protein_list:
            if not (sProtein.startswith(rev_str)):
                l.append(sProtein)
        self.ProteinNames = '{'+','.join(l) + '}'
        self.protein_list = l

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python xxx.py [options]

Inputs:
    input yyy
    output zzz

Options:
    -h/--help
    -v/--version

Outputs:
    output zzz
'''

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:n:",
                                    ["help",
                                     "version",
                                     "input",
                                     "config",
                                     "output",
                                     "negative"])

    # Default working dir and config file
    input_file = ""
    output_folder = ""
    config_file = ""
    negative_file = None

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "xxx.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input"):
            input_file = value
        if option in ("-o", "--output"):
            output_folder = value
        if option in ("-c", "--config"):
            config_file = value
        if option in ("-n", "--negative"):
            negative_file = value

    if input_file == "" or output_folder == "":
        print help_message
        sys.exit(0)

    output_folder = os.path.join(output_folder, '')

    return (input_file, config_file, output_folder, negative_file)

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    sProteins = protein_sequence.replace('{', '')
    sProteins = sProteins.replace('}', '')
    asProteins = sProteins.split(',')
    if lProtein != None:
        del lProtein[:]
        lProtein.extend(asProteins[:])
    for sProtein in asProteins:
        if not (sProtein.startswith(rev_str) or sProtein.startswith(shu_str)):
            return LabelFwd
    for sProtein in asProteins:
        if sProtein.startswith(shu_str):
            return LabelShu
    return LabelRev

# # split based on score agreement
def categorize_score_agreement(val):
    if val == 3:
        return '3'
    elif val == 2:
        return '2'
    else:
        return '1'

one_top_list = [0, 1, 2, 4]
# # split based on score agreement
def categorize_score_agreement_more_info(val):
    if val in one_top_list:
        return '1'
    elif val == 7:
        return 'A'
    elif val == 3:
        return 'XW'
    elif val == 5:
        return 'MW'
    elif val == 6:
        return 'MX'
    else:
        die("error")

two_top_list = [3, 5, 6]
# # split based on score agreement
def categorize_score_agreement_one_score_agreed(val):
    if val in two_top_list:
        return '2'
    elif val == 7:
        return '3'
    elif val == 0:
        return 'N'
    elif val == 1:
        return 'W'
    elif val == 2:
        return 'X'
    elif val == 4:
        return 'M'
    else:
        die("error")

# # split based on charge
def categorize_parent_charge(val):
    if val == 1:
        return '1'
    elif val == 2:
        return '2'
    elif val >= 3:
        return '3'

# # get the categories
def categorize(input_file, negative_file):

    psm_dict = {}
    # psm_set = set()
    
    # read line with csv
    f = csv.reader(CommentedFile(open(input_file, 'rb')),
                            delimiter='\t')
    # skip header
    _sHeader = f.next()
    # get data
    for sLine in f:
        #PsmFields_obj = PsmFields._make(sLine)
        if len(sLine) == 32:
            PsmFields_obj = PsmFields2._make(sLine)
        elif len(sLine) == 35:
            PsmFields_obj = PsmFields4._make(sLine)
        else:
            PsmFields_obj = PsmFields3._make(sLine)
        sParentCharge = categorize_parent_charge(int(PsmFields_obj.ParentCharge))
        sScoreAgreement = categorize_score_agreement(int(PsmFields_obj.ScoreAgreement))
        # sScoreAgreement = categorize_score_agreement_more_info(int(PsmFields_obj.ScoreAgreement))
        # sScoreAgreement = categorize_score_agreement_one_score_agreed(int(PsmFields_obj.ScoreAgreement))
        '''
        if sParentCharge != '2' or  sScoreAgreement != '3':
            continue
        '''
        sKey = sParentCharge + '_' + sScoreAgreement
        psm_obj = PSM(PsmFields_obj)
        if psm_obj.fMassDiff > mass_tolerance:
            continue
        '''
        if psm_obj.ScoreAgreement == 1 and psm_obj.iLocalRank > 0:
            continue
        '''
        '''
        if 'I' in PsmFields_obj.OriginalPeptide:
            psm_set.add(PsmFields_obj.OriginalPeptide.replace('I', 'L'))
        else:
            psm_set.add(PsmFields_obj.OriginalPeptide)
        '''
        if not psm_dict.has_key(sKey):
            psm_list = []
            psm_list.append(psm_obj)
            psm_dict[sKey] = psm_list
        else:
            psm_list = psm_dict[sKey]
            psm_list.append(psm_obj)
    
    psm_neg_list = []
    '''
    if negative_file != None:
        f = csv.reader(CommentedFile(open(negative_file, 'rb')),
                            delimiter='\t')
        # skip header
        _sHeader = f.next()
        # get data
        for sLine in f:
            #PsmFields_obj = PsmFields._make(sLine)
            if len(sLine) == 32:
                PsmFields_obj = PsmFields2._make(sLine)
            else:
                PsmFields_obj = PsmFields3._make(sLine)
            
            psm_obj = PSM(PsmFields_obj)
            if psm_obj.fMassDiff > 0.09:
                continue
            
            if psm_obj.RealLabel != LabelShu:
                continue
            
            if PsmFields_obj.OriginalPeptide in psm_set:
                continue
            
            if 'I' in PsmFields_obj.OriginalPeptide:
                if (PsmFields_obj.OriginalPeptide.replace('I', 'L')) in psm_set:
                    continue
            psm_obj.RealLabel = LabelRev
            psm_obj.TrainingLabel = LabelNegative
            psm_neg_list.append(psm_obj)
    '''

    return (psm_dict, psm_neg_list)

# # Division error handling
divide = sipros_post_module.divide
FDR_parameter = 1.0

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


# using the global rank
def get_cutoff_rank_product2(FDR_threshold, lPsm):
    F_list = []
    T_list = []
    TT_list = []
    S_list = []
    for oPsm in lPsm:
        if oPsm.RealLabel != LabelRev:
            T_list.append(oPsm.fRankProduct)
            if oPsm.RealLabel == LabelShu:
                S_list.append(oPsm.fRankProduct)
            else:
                TT_list.append(oPsm.fRankProduct)
        else:
            F_list.append(oPsm.fRankProduct)
    F_list = np.array(F_list)
    T_list = np.array(T_list)
    TT_list = np.array(TT_list)
    S_list = np.array(S_list)
    prev_TF_num = 0
    final_cutoff_score = np.amin(T_list, axis=0)
    _final_accept = False
    T_num_best = 0
    TT_num_best = 0
    F_num_best = 0
    S_num_best = 0
    TF_num_best = 0
    for cutoff_score in T_list:
        F_num = (F_list <= cutoff_score).sum()
        T_num = (T_list <= cutoff_score).sum()
        TT_num = (TT_list <= cutoff_score).sum()
        S_num = (S_list <= cutoff_score).sum()
        # (F_num, T_num) = get_target_decoy_hit(lPsm, cutoff_score, iScoreIndex)
        TF_num = F_num + T_num
        TF_num = S_num + TT_num
        (FDR_accept, FDR_value) = FDR_calculator(S_num, TT_num)
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
            S_num_best = S_num
            TF_num_best = TF_num

    print "Q-value_cutoff:\t%f\t%d\t%d\t%d" % (final_cutoff_score, T_num_best - S_num_best, F_num_best, S_num_best)
    return final_cutoff_score


# using the q value rank
def get_cutoff_q_rank_product(FDR_threshold, lPsm):
    iNumScores = lPsm[0].iNumScores

    T_num = 0
    F_num = 0
    S_num = 0
    for oPsm in lPsm:
        if oPsm.RealLabel != LabelRev:
            T_num += 1
            if oPsm.RealLabel == LabelShu:
                S_num += 1
        else:
            F_num += 1
    print "Before Filtering:\t\t%d\t%d\t%d" % (T_num - S_num, F_num, S_num)

    for i in range(iNumScores):
        newlist = sorted(lPsm, key=lambda x: x.lfScores[i], reverse=True)
        T_num = 0
        F_num = 0
        for j in range(len(newlist)):
            if newlist[j].RealLabel != LabelRev:
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
        # fTemp = oPsm.lRanks[0] * oPsm.lRanks[1] * oPsm.lRanks[2]
        # oPsm.fRankProduct = np.power(float(fTemp), 1.0 / 3.0)
        fTemp = (1 - ((1 - oPsm.lRanks[0]) * (1 - oPsm.lRanks[1]) * (1 - oPsm.lRanks[2])))
        oPsm.fRankProduct = fTemp

    final_cutoff_score = get_cutoff_rank_product2(FDR_threshold, lPsm)
    return final_cutoff_score

def mark_training_label(lPsm, final_cutoff_score, count_list=None, sKey=None):
    num_positive = 0
    for oPsm in lPsm:
        if oPsm.RealLabel == LabelRev:
            oPsm.TrainingLabel = LabelNegative
        else:
            '''
            if oPsm.ScoreAgreement < 2:
                if oPsm.iLocalRank == 0:
                    oPsm.TrainingLabel = LabelPositive
                else:
                    oPsm.TrainingLabel = LabelUnknown
                continue
            '''
            if oPsm.fRankProduct <= final_cutoff_score:
                oPsm.TrainingLabel = LabelPositive
                num_positive += 1
            else:
                oPsm.TrainingLabel = LabelUnknown
                oPsm.TrainingLabel = LabelPositive
                if sKey != None:
                    if sKey == '3_2' or sKey == '2_3' or sKey == '3_3':
                        oPsm.TrainingLabel = LabelUnknown
    '''            
    if num_positive >= 10000:
        list_sorted = sorted(lPsm, key=lambda x: x.fRankProduct)
        num_positive *= 0.95
        for oPsm in list_sorted:
            if num_positive > 0:
                if oPsm.TrainingLabel == LabelPositive:
                    num_positive -= 1
            else:
                if oPsm.TrainingLabel == LabelPositive:
                    oPsm.TrainingLabel = LabelUnknown
    '''
    if count_list != None:
        for oPsm in lPsm:
            if oPsm.fRankProduct <= final_cutoff_score:
                if oPsm.RealLabel == LabelFwd:
                    count_list[0] += 1
                elif oPsm.RealLabel == LabelRev:
                    count_list[1] += 1
                else:
                    count_list[2] += 1

def show_TP_TN_FP_FN(label_np, predict_np):
    true_np = (label_np == predict_np)
    TP = label_np[true_np].sum()
    TN = (label_np[true_np] == 0).sum()
    false_np = (label_np != predict_np)
    FP = (label_np[false_np] == 0).sum()
    FN = label_np[false_np].sum()
    print "TP\tTN\tFP\tFN"
    print "%d\t%d\t%d\t%d" % (TP, TN, FP, FN)
    
def show_Fdr_category(psm_dict):
    Best_last_list = [0, 0, 0]
    for _key, lPsm in psm_dict.iteritems():
        #print _key
        list_sorted = sorted(lPsm, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
        T_num = 0
        F_num = 0
        Fwd_num = 0
        Rev_num = 0
        Shu_num = 0
        Best_list = [0, 0, 0]
        for oPsm in list_sorted:
            if oPsm.RealLabel == LabelFwd:
                T_num += 1
                Fwd_num += 1
            elif oPsm.RealLabel == LabelRev:
                F_num += 1
                Rev_num += 1
            else:
                T_num += 1
                Shu_num += 1
            (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
            if (FDR_accept is True) and (FDR_value <= 0.01): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
                if Best_list[0] < Fwd_num:
                    Best_list = [Fwd_num, Rev_num, Shu_num]
        #sys.stdout.write('%d\t%d\t%d\n' % (Best_list[0], Best_list[1], Best_list[2]))
        Best_last_list = [Best_last_list[0] + Best_list[0], Best_last_list[1] + Best_list[1], Best_last_list[2] + Best_list[2]]
    sys.stdout.write('%d\t%d\t%d\t' % (Best_last_list[0], Best_last_list[1], Best_last_list[2]))
    pass

def filter_Fdr(psm_list, fdr):
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    psm_filtered_list = []
    psm_left_list = []
    cutoff_probability = 0.0
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
    
    sys.stdout.write('\n%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
            oPsm.TrainingLabel = LabelFiltered
        else:
            psm_left_list.append(oPsm)
    
    return (psm_filtered_list, psm_left_list)

def filter_Fdr2(psm_list, fdr):
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    psm_filtered_list = []
    psm_left_list = []
    cutoff_probability = 0.0
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fRtPvalue < 0.0001:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
    
    sys.stdout.write('\n%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
            oPsm.TrainingLabel = LabelFiltered
        else:
            psm_left_list.append(oPsm)
    
    return (psm_filtered_list, psm_left_list)

def show_Fdr_charge(psm_list):
    psm_new_list_1 = []
    psm_new_list_2 = []
    for oPsm in psm_list:
        if oPsm.ParentCharge <= 2:
            psm_new_list_1.append(oPsm)
        else:
            psm_new_list_2.append(oPsm)
    show_Fdr(psm_new_list_1, None, None)
    show_Fdr(psm_new_list_2, None, None)
                                  
def show_Fdr_varied(psm_list, fdr):
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]

    psm_filtered_list = []
    
    cutoff_probability = 0.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
            
    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
    
    # sys.stdout.write('\n'+str(len(psm_filtered_list))+'\t')
    
    sys.stdout.write('%d\t%d\t%d\t\n' % (Best_list[0], Best_list[1], Best_list[2]))
    
    return psm_filtered_list

def show_Fdr(psm_list, sKey, fdr=None):
    
    fdr_float = 0.01
    if fdr !=None:
        fdr_float = fdr
    
    # list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability, 1 - x.fRankProduct) , reverse=True)
    list_sorted = sorted(psm_list, key=lambda x: (x.fPredictProbability) , reverse=True)
    T_num = 0
    F_num = 0
    Fwd_num = 0
    Rev_num = 0
    Shu_num = 0
    Best_list = [0, 0, 0]
    '''
    fwr = open("fwr.txt", "a")
    rev = open("rev.txt", "a")
    '''
    psm_filtered_list = []
    cutoff_probability = 0.0
    # without considering training label
    for oPsm in list_sorted:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        if oPsm.RealLabel == LabelFwd:
            T_num += 1
            Fwd_num += 1
        elif oPsm.RealLabel == LabelRev:
            F_num += 1
            Rev_num += 1
        else:
            T_num += 1
            Shu_num += 1
        (FDR_accept, FDR_value) = FDR_calculator(Shu_num, Fwd_num)
        if (FDR_accept is True) and (FDR_value <= fdr_float): #(Shu_num <= shuffle_num_control[sKey]):  # and (Shu_num <= 618) and (FDR_value <= 0.01)
            if (Best_list[0] + Best_list[2]) < (Fwd_num + Shu_num):
                Best_list = [Fwd_num, Rev_num, Shu_num]
                cutoff_probability = oPsm.fPredictProbability
            
            '''
            if oPsm.RealLabel == LabelFwd:
                fwr.write(str(oPsm.NMC))
                fwr.write('\n')
            else:
                rev.write(str(oPsm.NMC))
                rev.write('\n')
        else:
            rev.write(str(oPsm.NMC))
            rev.write('\n')
           
    fwr.close()
    rev.close()
    '''
    for oPsm in list_sorted:
        if oPsm.fPredictProbability >= cutoff_probability:
            psm_filtered_list.append(oPsm)
    '''
    # # debug
    file_temp_str = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/Ensamble/DP_fwr.txt"
    file_temp_Dec_str = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/Ensamble/DP_dec.txt"
    fwD = open(file_temp_Dec_str, 'w')
    with open(file_temp_str,'w') as fw:
        for oPsm in psm_filtered_list:
            if oPsm.RealLabel == LabelFwd:
                fw.write(oPsm.FileName)
                fw.write('_')
                fw.write(str(oPsm.ScanNumber))
                fw.write('\n')
            else:
                fwD.write(oPsm.FileName)
                fwD.write('_')
                fwD.write(str(oPsm.ScanNumber))
                fwD.write('\n')
    fwD.close()
    # # debug
    '''        
    sys.stdout.write('\t'+str(len(psm_filtered_list))+'\n')
    
    sys.stdout.write('%d\t%d\t%d\t' % (Best_list[0], Best_list[1], Best_list[2]))
    
    return psm_filtered_list

def get_RT_data(psm_filtered_list, output_folder, psm_list):
    
    shift_threshold_float = 7.5
    
    # shift the RT
    ft2_filename_list = []
    ft2_filename_max_rt_dict = {}
    for oPsm in psm_list:
        if oPsm.FileName not in ft2_filename_max_rt_dict:
            ft2_filename_list.append(oPsm.FileName)
            ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
        else:
            if oPsm.fRtMeasured > ft2_filename_max_rt_dict[oPsm.FileName]:
                ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
    
    ft2_filename_list = sorted(ft2_filename_list)
    
    for oPsm in psm_list:
        if oPsm.fRtMeasured <= shift_threshold_float and oPsm.FileName != ft2_filename_list[0]:
            oPsm.bFileNameChanged = True
            ind_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[ind_filename_int - 1]
            oPsm.fRtMeasured += ft2_filename_max_rt_dict[oPsm.FileName]
            oPsm.sRTime = "%.5f" % oPsm.fRtMeasured
        
    ft2_train_test_model_result_list = []
    train_file_dict = {}
    test_file_dict = {}
    train_num_dict = {}
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRev:
            continue
        if oPsm.FileName in train_file_dict:
            if train_num_dict[oPsm.FileName] > 1000:
                continue
            fw = train_file_dict[oPsm.FileName]
            fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '(Oxidation)'))
            fw.write('\t')
            fw.write(oPsm.sRTime)
            fw.write('\n')
            train_num_dict[oPsm.FileName] += 1
        else:
            train_file_str = output_folder + oPsm.FileName.split('.')[0] + '_train.txt'
            test_file_str = output_folder + oPsm.FileName.split('.')[0] + '_test.txt'
            model_file_str = output_folder + oPsm.FileName.split('.')[0] + '_model.txt'
            result_file_str = output_folder + oPsm.FileName.split('.')[0] + '_result.csv'
            ft2_train_test_model_result_list.append((oPsm.FileName, train_file_str, test_file_str, model_file_str, result_file_str))
            fw_test = open(test_file_str, 'w')
            test_file_dict[oPsm.FileName] = fw_test
            fw = open(train_file_str, 'w')
            train_file_dict[oPsm.FileName] = (fw)
            fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '(Oxidation)'))
            fw.write('\t')
            fw.write(oPsm.sRTime)
            fw.write('\n')
            train_num_dict[oPsm.FileName] = 1
    
    for _key, value in train_file_dict.iteritems():
        value.close()
    
    for oPsm in psm_list:
        if oPsm.FileName in test_file_dict:
            fw_test = test_file_dict[oPsm.FileName]
            fw_test.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '(Oxidation)'))
            fw_test.write('\n')
            
    for _key, value in test_file_dict.iteritems():
        value.close()
        
    PGM = "/home/xgo/ORNL/Sipros/MachineLearning/RetentionTime/OpenMS/OpenMS-build/bin/"
    iNumThreads = cpu_count() - 1
    qUnprocessed = Queue(len(ft2_train_test_model_result_list)*2)
    qProcessed = Queue(len(ft2_train_test_model_result_list)*2)
    
    # run RTModel
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"RTModel", '-in', train_file_str, '-out', model_file_str, '-cv:skip_cv']
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTModel(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    '''
    '''
    # run RTPredict
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"RTPredict", '-in_text', test_file_str, '-svm_model', model_file_str, '-out_text:file', result_file_str]
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTPredict(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    
    # collect RTtime
    ft_pep_rtime_dict = {}
    for (ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        pep_rtime_dict = {}
        ft_pep_rtime_dict[ft_file_str] = pep_rtime_dict
        with open(result_file_str, 'r') as fr:
            for line_str in fr:
                words = line_str.strip().split()
                pep_str = '[' + words[0].replace('(Oxidation)', '~') + ']'
                if pep_str in pep_rtime_dict:
                    print 'error in pep_str in pep_rtime_dict'
                else:
                    rt_float = float(words[1])
                    pep_rtime_dict[pep_str] = rt_float
    
    # assign rtime to psm, and roll back the file name
    psm_na_rt_list = []
    for oPsm in psm_list:
        if oPsm.FileName not in ft_pep_rtime_dict:
            psm_na_rt_list.append(oPsm)
        else:
            pep_rtime_dict = ft_pep_rtime_dict[oPsm.FileName]
            oPsm.fRtPredict = float(pep_rtime_dict[oPsm.IdentifiedPeptide])
        
        if oPsm.bFileNameChanged:
            idx_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[idx_filename_int + 1]
            oPsm.bFileNameChanged = False
        
        
    range_list = []
    threshold_list = []
    distribuation_list = []
    num_range = 6
    scale_range = 20
    for i in range(1, num_range+1):
        range_list.append([])
        threshold_list.append(i* scale_range)
    threshold_list[num_range - 1] = 125
    for oPsm in psm_filtered_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        range_list[idx].append(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        
    for value_list in range_list:
        value_list_np = np.array(value_list)
        std_float = np.std(value_list_np)
        mean_float = np.mean(value_list_np)
        distribuation_list.append((mean_float, std_float))
    
    for oPsm in psm_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        oPsm.fRtPvalue = spst.norm(distribuation_list[idx][0], distribuation_list[idx][1]).cdf(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        oPsm.fRtPvalue = 1 - oPsm.fRtPvalue

def get_RT_elude_data(psm_filtered_list, output_folder, psm_list):
    
    shift_threshold_float = 7.5
    
    # shift the RT
    ft2_filename_list = []
    ft2_filename_max_rt_dict = {}
    for oPsm in psm_list:
        if oPsm.FileName not in ft2_filename_max_rt_dict:
            ft2_filename_list.append(oPsm.FileName)
            ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
        else:
            if oPsm.fRtMeasured > ft2_filename_max_rt_dict[oPsm.FileName]:
                ft2_filename_max_rt_dict[oPsm.FileName] = oPsm.fRtMeasured
    
    ft2_filename_list = sorted(ft2_filename_list)
    
    for oPsm in psm_list:
        if oPsm.fRtMeasured <= shift_threshold_float and oPsm.FileName != ft2_filename_list[0]:
            oPsm.bFileNameChanged = True
            ind_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[ind_filename_int - 1]
            oPsm.fRtMeasured += ft2_filename_max_rt_dict[oPsm.FileName]
            oPsm.sRTime = "%.5f" % oPsm.fRtMeasured
        
    ft2_train_test_model_result_list = []
    train_file_dict = {}
    test_file_dict = {}
    train_num_dict = {}
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRev:
            continue
        if oPsm.FileName in train_file_dict:
            if train_num_dict[oPsm.FileName] > 400:
                continue
            fw = train_file_dict[oPsm.FileName]
            fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '[Oxidation]'))
            fw.write('\t')
            fw.write(oPsm.sRTime)
            fw.write('\n')
            train_num_dict[oPsm.FileName] += 1
        else:
            train_file_str = output_folder + oPsm.FileName.split('.')[0] + '_train.txt'
            test_file_str = output_folder + oPsm.FileName.split('.')[0] + '_test.txt'
            model_file_str = output_folder + oPsm.FileName.split('.')[0] + '_model.txt'
            result_file_str = output_folder + oPsm.FileName.split('.')[0] + '_result.txt'
            ft2_train_test_model_result_list.append((oPsm.FileName, train_file_str, test_file_str, model_file_str, result_file_str))
            fw_test = open(test_file_str, 'w')
            test_file_dict[oPsm.FileName] = fw_test
            fw = open(train_file_str, 'w')
            train_file_dict[oPsm.FileName] = (fw)
            fw.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '[Oxidation]'))
            fw.write('\t')
            fw.write(oPsm.sRTime)
            fw.write('\n')
            train_num_dict[oPsm.FileName] = 1
    
    for _key, value in train_file_dict.iteritems():
        value.close()
    
    for oPsm in psm_list:
        if oPsm.FileName in test_file_dict:
            fw_test = test_file_dict[oPsm.FileName]
            fw_test.write((oPsm.IdentifiedPeptide[1:-1]).replace('~', '[Oxidation]'))
            fw_test.write('\n')
            
    for _key, value in test_file_dict.iteritems():
        value.close()
        
    PGM = "/usr/bin/"
    iNumThreads = cpu_count() - 1
    qUnprocessed = Queue(len(ft2_train_test_model_result_list)*2)
    qProcessed = Queue(len(ft2_train_test_model_result_list)*2)
    
    # run Elude
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"elude", '-t', train_file_str, '-s', model_file_str]
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTModel(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    '''
    '''
    # run RTPredict
    for (_ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        call_list = [PGM+"elude", '-e', test_file_str, '-l', model_file_str, '-o', result_file_str]
        qUnprocessed.put(call_list)
    for _i in range(iNumThreads):
        PsmProcessor = RTPredict(qUnprocessed, qProcessed)
        PsmProcessor.daemon = True
        PsmProcessor.start()
        qUnprocessed.put(None)
    
    iNumRankers = iNumThreads
    while True:
        call_list = qProcessed.get(True)
        if call_list is None:
            iNumRankers -= 1
            if iNumRankers == 0:
                break
            else:
                continue
    
    # collect RTtime
    ft_pep_rtime_dict = {}
    for (ft_file_str, train_file_str, test_file_str, model_file_str, result_file_str) in ft2_train_test_model_result_list:
        pep_rtime_dict = {}
        ft_pep_rtime_dict[ft_file_str] = pep_rtime_dict
        with open(result_file_str, 'r') as fr:
            fr.next()
            fr.next()
            fr.next()
            for line_str in fr:
                words = line_str.strip().split()
                pep_str = '[' + words[0].replace('[Oxidation]', '~') + ']'
                if pep_str in pep_rtime_dict:
                    pass
                    # print 'error in pep_str in pep_rtime_dict'
                else:
                    rt_float = float(words[1])
                    pep_rtime_dict[pep_str] = rt_float
    
    # assign rtime to psm, and roll back the file name
    psm_na_rt_list = []
    for oPsm in psm_list:
        if oPsm.FileName not in ft_pep_rtime_dict:
            psm_na_rt_list.append(oPsm)
        else:
            pep_rtime_dict = ft_pep_rtime_dict[oPsm.FileName]
            if oPsm.IdentifiedPeptide not in pep_rtime_dict:
                print 'error in pep_rtime_dict.'
            oPsm.fRtPredict = float(pep_rtime_dict[oPsm.IdentifiedPeptide])
        
        if oPsm.bFileNameChanged:
            idx_filename_int = ft2_filename_list.index(oPsm.FileName)
            oPsm.FileName = ft2_filename_list[idx_filename_int + 1]
            oPsm.bFileNameChanged = False
    
    if len(psm_na_rt_list) > 0:
        print 'error in psm_na_rt_list.'
        
    range_list = []
    threshold_list = []
    distribuation_list = []
    num_range = 6
    scale_range = 20
    for i in range(1, num_range+1):
        range_list.append([])
        threshold_list.append(i* scale_range)
    threshold_list[num_range - 1] = 125
    for oPsm in psm_filtered_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        range_list[idx].append(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        
    for value_list in range_list:
        value_list_np = np.array(value_list)
        std_float = np.std(value_list_np)
        mean_float = np.mean(value_list_np)
        distribuation_list.append((mean_float, std_float))
    
    for oPsm in psm_list:
        idx = 0
        for i in range(num_range):
            if oPsm.fRtMeasured <= threshold_list[i]:
                idx = i
                break
        oPsm.fRtPvalue = spst.norm(distribuation_list[idx][0], distribuation_list[idx][1]).cdf(abs(oPsm.fRtMeasured - oPsm.fRtPredict))
        oPsm.fRtPvalue = 1 - oPsm.fRtPvalue

# # thread class for ranking the PSM
class RTModel(Process):

    def __init__(self, qUnprocessed, qProcessed):
        super(RTModel, self).__init__()
        self.qUnprocessed = qUnprocessed
        self.qProcessed = qProcessed
        return

    def run(self):
        while True:
            call_list = self.qUnprocessed.get(True)
            if call_list is None:
                break
            call(call_list)
        self.qProcessed.put(None)
        return

# # thread class for ranking the PSM
class RTPredict(Process):

    def __init__(self, qUnprocessed, qProcessed):
        super(RTPredict, self).__init__()
        self.qUnprocessed = qUnprocessed
        self.qProcessed = qProcessed
        return

    def run(self):
        while True:
            call_list = self.qUnprocessed.get(True)
            if call_list is None:
                break
            call(call_list)
        self.qProcessed.put(None)
        return
'''
from keras.models import Sequential    
from keras.layers import Dense, Dropout

def test_DeepLearning(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    num_positive = float((label_np==LabelPositive).sum())
    num_negative = float((label_np==LabelNegative).sum())
    class_weight_dict = {0: (num_positive/(num_negative+num_positive)), 1:(num_negative/(num_negative+num_positive))}
   
    clf = Sequential()
    
    
    clf.add(Dense(32, input_dim=len(feature_selection_list), init='uniform', activation='sigmoid'))
    clf.add(Dropout(0))
    clf.add(Dense(32, init='uniform', activation='sigmoid'))
    clf.add(Dropout(0))
    clf.add(Dense(1, init='uniform', activation='sigmoid'))
    
    # sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)
    clf.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['accuracy'])
    # clf.compile(loss='categorical_crossentropy', optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True))
    
    clf.fit(train_data_np, label_np, batch_size=32, nb_epoch=2, verbose=1)#,class_weight=class_weight_dict)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i]
        psm_list[id_list[i]].ML_feature.append(predict_np[i])
    fdr_rank(psm_list, 3) 
    return None
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list
'''

from sklearn.ensemble import AdaBoostClassifier

def test_AdaBoost(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = AdaBoostClassifier(# base_estimator=DecisionTreeClassifier(min_samples_split=800,min_samples_leaf=100),
                             n_estimators=200,
                             random_state=50)
    
    clf.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        psm_list[id_list[i]].ML_feature.append(predict_np[i, 1])
    fdr_rank(psm_list, 2)
    return None
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    
    sys.stdout.write('\n')
        
    return psm_filtered_list
    
    
from sklearn.svm import SVC

def test_SVM(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = SVC(random_state=50,
              kernel='linear',
              class_weight='balanced',
              cache_size=4000,
              C=1,
              gamma='auto')
    clf.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.decision_function(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    '''
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    '''    
    sys.stdout.write('\n')
        
    return psm_filtered_list

from sklearn.naive_bayes import GaussianNB

def test_naive_bayes(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    train_data_np = preprocessing.normalize(train_data_np, norm='l2',axis=0)
    
    clf = GaussianNB()
    clf.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    test_unknown_np = preprocessing.normalize(test_unknown_np, norm='l2',axis=0)
    
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    '''
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    '''    
    sys.stdout.write('\n')
        
    return psm_filtered_list

def test_random_forest(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = RandomForestClassifier(n_estimators=200, 
                                 n_jobs=-1, 
                                 class_weight='balanced', 
                                 oob_score=True,
                                 bootstrap=True,
                                 random_state=50,
                                 criterion="gini",
                                 max_features="auto",
                                 min_samples_leaf=50,
                                 min_samples_split=800)
    
    clf.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        psm_list[id_list[i]].ML_feature.append(predict_np[i, 1])
    fdr_rank(psm_list, 1)
    return None
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % clf.feature_importances_[idx])
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
        
    return psm_filtered_list

def test_stacking_MOE(psm_list):
    # # construct training data
    psm_filtered_list = []
    data_list = []
    label_list = []
    id_list = []
    D_list = []
    U_list = []
    D_id_list = []
    U_id_list = []     
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
        if oPsm.TrainingLabel == LabelNegative:
            D_list.append(oPsm.feature_list)
            D_id_list.append(oPsm.iInnerId)
        else:
            U_list.append(oPsm.feature_list)
            U_id_list.append(oPsm.iInnerId)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\n")
    label_np = np.array(label_list)
    
    feature_list_local = [1, 2, 3, 16]
    D_np = np.array(D_list, dtype=np.float64)
    U_np = np.array(U_list, dtype=np.float64)
    
    # features and models
    feature_matrix = [[1, 2, 3, 15, 16, 17],    # score section
                      [5, 29],                  # mass section
                      [24],                     # digestion section
                      [30],                     # PTM section
                      [25, 26, 27, 28]]         # pep pro support section
    tier_1 = []
    D_phix = []
    U_phix = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
        predict_np = tier_1[i].predict_proba(D_np[:, feature_matrix[i]])
        D_phix.append(predict_np[:, 0])
        predict_np = tier_1[i].predict_proba(U_np[:, feature_matrix[i]])
        U_phix.append(predict_np[:, 1])
    D_phi_np = np.transpose(np.array(D_phix))
    U_phi_np = np.transpose(np.array(U_phix))
    D_np = D_np[:, feature_list_local]
    U_np = U_np[:, feature_list_local]
    (D_pyx, U_pyx) = test_MoE_semi(D_np, U_np, D_phi_np, U_phi_np)
    # D_pyx = D_phi_np[:,1]
    # U_pyx = U_phi_np[:,1]
    
    for i in range(len(D_id_list)):
        psm_list[D_id_list[i]].fPredictProbability = 1 - D_pyx[i]
    
    for i in range(len(U_id_list)):
        psm_list[U_id_list[i]].fPredictProbability = U_pyx[i]
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
    
    return psm_filtered_list

def test_MoE_super(D_np, F_np, D_phi_np, F_phi_np):
    # normalize 
    preprocessing.normalize(D_np, axis=0,copy=False, norm='max')
    preprocessing.normalize(F_np, axis=0,copy=False, norm='max')
    # preprocessing.normalize(D_phi_np, axis=0,copy=False, norm='max')
    # preprocessing.normalize(F_phi_np, axis=0,copy=False, norm='max')
    # n sample #, m expert #, d feature #
    # initialize parameters: alpha, mu, sigma_square
    numD = np.float128(D_np.shape[0])
    numF = np.float128(F_np.shape[0])
    num = numD + numF
    d = np.float128(D_np.shape[1])
    m = np.float128(D_phi_np.shape[1])
    alpha = np.zeros(m, dtype=np.float128) + 1.0/float(m)
    mu_D = np.average(D_np, axis=0)
    mu_U = np.average(F_np, axis=0)
    mu_temp = (mu_D + mu_U)/2.0
    mu = np.tile(mu_temp, (m, 1))
    sigma_square = np.zeros(m, dtype=np.float128) + 1.0
    pi = math.pi
    
    # alpha = p(z,y|theta) 
    alpha_old = np.copy(alpha)
    mu_old = np.copy(mu)
    sigma_square_old = np.copy(sigma_square)
    
    # the iteration of E-M algorithm
    for _ite in range(100):
        # E-step
        D_minus_mu = []
        F_minus_mu = []
        for mu_j in mu:
            D_minus_mu.append(np.sum(np.multiply((D_np - mu_j), (D_np - mu_j)), axis=1)) # (1, n)
            F_minus_mu.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
        D_minus_mu = np.array(D_minus_mu, dtype=np.float128) # (m, n)
        D_minus_mu = np.transpose(D_minus_mu) # (n, m)
        F_minus_mu = np.array(F_minus_mu, dtype=np.float128) # (m, n)
        F_minus_mu = np.transpose(F_minus_mu) # (n, m)
        
        D_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(D_minus_mu, (-2*sigma_square))))
        D_alpha_px = np.multiply(alpha, D_px) # (n, m)
        D_phi_alpha_px = np.multiply(D_alpha_px , D_phi_np) # (n, m)
        D_phi_alpha_px_sum = np.sum(D_phi_alpha_px, axis=1)
        # # Q(z=j)
        D_hzj = np.divide(D_phi_alpha_px , D_phi_alpha_px_sum[:, None]) # (n, m)
        
        F_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(F_minus_mu, (-2*sigma_square))))
        F_alpha_px = np.multiply(alpha, F_px) # (n, m)
        F_phi_alpha_px = np.multiply(F_alpha_px , F_phi_np) # (n, m)
        F_phi_alpha_px_sum = np.sum(F_phi_alpha_px, axis=1)
        # # Q(z=j)
        F_hzj = F_phi_alpha_px / F_phi_alpha_px_sum[:,None] # (n, m)
        
        # log-likelihood
        part1 = np.sum(np.log(D_phi_alpha_px_sum), axis=0)
        part2 = np.sum(np.log(F_phi_alpha_px_sum), axis=0)
        ftheta = part1 + part2
        print str(ftheta)
        
        # M -step
        # # alpha
        alpha = np.divide((np.sum(D_hzj, axis=0) + np.sum(F_hzj, axis=0)), (num))
        # # mu_j, sigma_j
        for j in range(int(m)):
            # sigma square
            sigma_square[j] = (np.sum(D_hzj[:,j]*np.sum(np.multiply((D_np - mu[j]), (D_np - mu[j])), axis = 1)) + 
                               np.sum(F_hzj[:,j]*np.sum(np.multiply((F_np - mu[j]), (F_np - mu[j])), axis = 1)))
            sigma_square[j] /= (d*alpha[j] * num)
            # mu
            mu[j] = (np.sum(D_np*((D_hzj[:,j])[:,None]), axis=0)+np.sum(F_np*((F_hzj[:,j])[:,None]), axis=0))/(alpha[j] * num)
        
        # check if converge
        error_sum = np.sum(np.absolute(mu-mu_old)) + np.sum(np.absolute(alpha-alpha_old)) + np.sum(np.absolute(sigma_square-sigma_square_old))
        error_sum /= float(mu.size + alpha.size + sigma_square.size)
        alpha_old = np.copy(alpha)
        mu_old = np.copy(mu)
        sigma_square_old = np.copy(sigma_square)
        
        if error_sum < 0.0001:
            print 'Converged:\t%d,\tError:\t%f' % (_ite, error_sum)
            '''
            print alpha
            print mu
            print sigma_square
            '''
            break
        else:
            print 'Iteration:\t%d,\tError:\t%f' % (_ite, error_sum)
            '''
            print alpha
            print mu
            print sigma_square
            '''
    
    # assign probability
    D_minus_mu = []
    F_minus_mu = []
    for mu_j in mu:
        D_minus_mu.append(np.sum((D_np - mu_j)*(D_np - mu_j), axis=1)) # (1, n)
        F_minus_mu.append(np.sum((F_np - mu_j)*(F_np - mu_j), axis=1))
    D_minus_mu = np.array(D_minus_mu, dtype=np.float128) # (m, n)
    D_minus_mu = np.transpose(D_minus_mu) # (n, m)
    F_minus_mu = np.array(F_minus_mu, dtype=np.float128) # (m, n)
    F_minus_mu = np.transpose(F_minus_mu) # (n, m)
    D_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(D_minus_mu, (-2*sigma_square))))
    D_alpha_px = np.multiply(alpha, D_px) # (n, m)
    D_alpha_px_sum = np.sum(D_alpha_px, axis=1)
    D_gx = D_alpha_px / D_alpha_px_sum[:,None]
    D_phi_gx = D_phi_np * D_gx
    D_pyx = np.sum(D_phi_gx, axis=1)
    
    F_px = np.multiply((np.power(2*pi*sigma_square, -d/2)), np.exp(np.divide(F_minus_mu, (-2*sigma_square))))
    F_alpha_px = np.multiply(alpha, F_px) # (n, m)
    F_alpha_px_sum = np.sum(F_alpha_px, axis=1)
    F_gx = F_alpha_px / F_alpha_px_sum[:, None]
    F_phi_gx = F_phi_np * F_gx
    U_pyx = np.sum(F_phi_gx, axis=1)
    
    return (D_pyx, U_pyx)

def test_MoE_semi(D_np, F_np, D_phi_np, F_phi_np):
    # normalize 
    preprocessing.normalize(D_np, axis=0,copy=False, norm='max')
    preprocessing.normalize(F_np, axis=0,copy=False, norm='max')
    F_phi_2_np = F_phi_np
    F_phi_1_np = 1.0 - F_phi_np 
    # preprocessing.normalize(D_phi_np, axis=0,copy=False, norm='max')
    # preprocessing.normalize(F_phi_np, axis=0,copy=False, norm='max')
    # n sample #, m expert #, d feature #
    # initialize parameters: alpha_j1, mu_j1, sigma_square_j1
    numD = np.float128(D_np.shape[0])
    numF = np.float128(F_np.shape[0])
    num = numD + numF
    d = np.float128(D_np.shape[1])
    m = np.float128(D_phi_np.shape[1])
    alpha_j1 = np.zeros(m, dtype=np.float128) + 1.0/float(m)
    alpha_j2 = np.zeros(m, dtype=np.float128) + 1.0/float(m)
    mu_D = np.average(D_np, axis=0)
    mu_U = np.average(F_np, axis=0)
    mu_temp = (mu_D + mu_U)/2.0
    mu_j1 = np.tile(mu_temp, (m, 1))
    mu_j2 = np.tile(mu_temp, (m, 1))
    sigma_square_j1 = np.zeros(m, dtype=np.float128) + 1.0
    sigma_square_j2 = np.zeros(m, dtype=np.float128) + 1.0
    pi = math.pi
    
    # alpha = p(z,y|theta) 
    alpha_j1_old = np.copy(alpha_j1)
    alpha_j2_old = np.copy(alpha_j2)
    mu_j1_old = np.copy(mu_j1)
    mu_j2_old = np.copy(mu_j2)
    sigma_square_j1_old = np.copy(sigma_square_j1)
    sigma_square_j2_old = np.copy(sigma_square_j2)
    
    # the iteration of E-M algorithm
    for _ite in range(100):
        # E-step
        D_minus_mu_j1 = []
        F_minus_mu_j1 = []
        F_minus_mu_j2 = []
        for mu_j in mu_j1:
            D_minus_mu_j1.append(np.sum(np.multiply((D_np - mu_j), (D_np - mu_j)), axis=1)) # (1, n)
            F_minus_mu_j1.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
        for mu_j in mu_j2:
            F_minus_mu_j2.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
        D_minus_mu_j1 = np.array(D_minus_mu_j1, dtype=np.float128) # (m, n)
        D_minus_mu_j1 = np.transpose(D_minus_mu_j1) # (n, m)
        F_minus_mu_j1 = np.array(F_minus_mu_j1, dtype=np.float128) # (m, n)
        F_minus_mu_j1 = np.transpose(F_minus_mu_j1) # (n, m)
        F_minus_mu_j2 = np.array(F_minus_mu_j2, dtype=np.float128) # (m, n)
        F_minus_mu_j2 = np.transpose(F_minus_mu_j2) # (n, m)
        D_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(D_minus_mu_j1, (-2*sigma_square_j1))))
        D_alpha_px_zy1 = np.multiply(alpha_j1, D_px_zy1) # (n, m)
        D_phi_alpha_px_zy1 = np.multiply(D_alpha_px_zy1 , D_phi_np) # (n, m)
        D_phi_alpha_px_zy1_sum = np.sum(D_phi_alpha_px_zy1, axis=1)
        # # Q(z=j)
        Qzj = np.divide(D_phi_alpha_px_zy1 , D_phi_alpha_px_zy1_sum[:, None]) # (n, m)
        # D_alphaPx_log = np.exp(np.log(alpha_j1) + np.log(pow(2*pi*sigma_square_j1, d/2)) + (D_minus_mu_j1*(-1)/(2*sigma_square_j1))) # (n, m)
        F_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(F_minus_mu_j1, (-2*sigma_square_j1))))
        F_px_zy2 = np.multiply((np.power(2*pi*sigma_square_j2, -d/2)), np.exp(np.divide(F_minus_mu_j2, (-2*sigma_square_j2))))
        F_alpha_px_zy1 = np.multiply(alpha_j1, F_px_zy1) # (n, m)
        F_alpha_px_zy2 = np.multiply(alpha_j2, F_px_zy2) # (n, m)
        F_phi_alpha_px_zy1 = np.multiply(F_alpha_px_zy1 , F_phi_1_np) # (n, m)
        F_phi_alpha_px_zy2 = np.multiply(F_alpha_px_zy2 , F_phi_2_np) # (n, m)
        F_phi_alpha_px_zy1_sum = np.sum(F_phi_alpha_px_zy1, axis=1)
        F_phi_alpha_px_zy2_sum = np.sum(F_phi_alpha_px_zy2, axis=1)
        F_phi_alpha_px_zy_sum = F_phi_alpha_px_zy1_sum + F_phi_alpha_px_zy2_sum
        Qzjy1 = F_phi_alpha_px_zy1 / F_phi_alpha_px_zy_sum[:,None] # (n, m)
        Qzjy2 = F_phi_alpha_px_zy2 / F_phi_alpha_px_zy_sum[:,None] # (n, m)
        
        # log-likelihood
        part1 = np.sum(np.log(D_phi_alpha_px_zy1_sum), axis=0)
        part2 = np.sum(np.log(F_phi_alpha_px_zy_sum), axis=0)
        ftheta = part1 + part2
        print str(ftheta)
        
        # M -step
        # # alpha
        alpha_j1 = np.divide((np.sum(Qzj, axis=0) + np.sum(Qzjy1, axis=0)), (num))
        alpha_j2 = np.divide(np.sum(Qzjy2, axis=0), (num))
        # # mu_j, sigma_j
        for j in range(int(m)):
            # sigma square
            sigma_square_j1[j] = (np.sum(Qzj[:,j]*np.sum(np.multiply((D_np - mu_j1[j]), (D_np - mu_j1[j])), axis=1)) + 
                               np.sum(Qzjy1[:,j]*np.sum(np.multiply((F_np - mu_j1[j]), (F_np - mu_j1[j])), axis = 1)))
            sigma_square_j1[j] /= (d*alpha_j1[j] * num)
            sigma_square_j2[j] = np.sum(Qzjy2[:,j]*np.sum(np.multiply((F_np - mu_j2[j]), (F_np - mu_j2[j])), axis = 1))
            sigma_square_j2[j] /= (d*alpha_j2[j] * num)
            # mu
            mu_j1[j] = (np.sum(D_np*((Qzj[:,j])[:,None]), axis=0)+np.sum(F_np*((Qzjy1[:,j])[:,None]), axis=0))/(alpha_j1[j] * num)
            mu_j2[j] = (np.sum(F_np*((Qzjy2[:,j])[:,None]), axis=0))/(alpha_j2[j] * num)
        
        # check if converge
        error_sum = np.sum(np.absolute(mu_j1-mu_j1_old)) + np.sum(np.absolute(alpha_j1-alpha_j1_old)) + np.sum(np.absolute(sigma_square_j1-sigma_square_j1_old))
        error_sum += np.sum(np.absolute(mu_j2-mu_j2_old)) + np.sum(np.absolute(alpha_j2-alpha_j2_old)) + np.sum(np.absolute(sigma_square_j2-sigma_square_j2_old))
        error_sum /= 2.0*float(mu_j1.size + alpha_j1.size + sigma_square_j1.size)
        alpha_j1_old = np.copy(alpha_j1)
        mu_j1_old = np.copy(mu_j1)
        sigma_square_j1_old = np.copy(sigma_square_j1)
        alpha_j2_old = np.copy(alpha_j2)
        mu_j2_old = np.copy(mu_j2)
        sigma_square_j2_old = np.copy(sigma_square_j2)
        
        if error_sum < 0.0001:
            print 'Converged:\t%d,\tError:\t%f' % (_ite, error_sum)
            '''
            print alpha_j1
            print mu_j1
            print sigma_square_j1
            '''
            break
        else:
            print 'Iteration:\t%d,\tError:\t%f' % (_ite, error_sum)
            '''
            print alpha_j1
            print mu_j1
            print sigma_square_j1
            '''
    
    # assign probability
    D_minus_mu_j1 = []
    F_minus_mu_j1 = []
    F_minus_mu_j2 = []
    for mu_j in mu_j1:
        D_minus_mu_j1.append(np.sum((D_np - mu_j)*(D_np - mu_j), axis=1)) # (1, n)
        F_minus_mu_j1.append(np.sum((F_np - mu_j)*(F_np - mu_j), axis=1))
    for mu_j in mu_j2:
        F_minus_mu_j2.append(np.sum(np.multiply((F_np - mu_j), (F_np - mu_j)), axis=1))
    D_minus_mu_j1 = np.array(D_minus_mu_j1, dtype=np.float128) # (m, n)
    D_minus_mu_j1 = np.transpose(D_minus_mu_j1) # (n, m)
    F_minus_mu_j1 = np.array(F_minus_mu_j1, dtype=np.float128) # (m, n)
    F_minus_mu_j1 = np.transpose(F_minus_mu_j1) # (n, m)
    F_minus_mu_j2 = np.array(F_minus_mu_j2, dtype=np.float128) # (m, n)
    F_minus_mu_j2 = np.transpose(F_minus_mu_j2) # (n, m)
    D_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(D_minus_mu_j1, (-2*sigma_square_j1))))
    D_alpha_px_zy1 = np.multiply(alpha_j1, D_px_zy1) # (n, m)
    D_alpha_px_zy1_sum = np.sum(D_alpha_px_zy1, axis=1)
    D_pzj_x = D_alpha_px_zy1 / D_alpha_px_zy1_sum[:,None]
    D_phi_pzj_x = D_phi_np * D_pzj_x
    D_pyx = np.sum(D_phi_pzj_x, axis=1)
    
    F_px_zy1 = np.multiply((np.power(2*pi*sigma_square_j1, -d/2)), np.exp(np.divide(F_minus_mu_j1, (-2*sigma_square_j1))))
    F_px_zy2 = np.multiply((np.power(2*pi*sigma_square_j2, -d/2)), np.exp(np.divide(F_minus_mu_j2, (-2*sigma_square_j2))))
    F_alpha_px_zy1 = np.multiply(alpha_j1, F_px_zy1) # (n, m)
    F_alpha_px_zy2 = np.multiply(alpha_j2, F_px_zy2) # (n, m)
    F_alpha_px_zy = F_alpha_px_zy1 + F_alpha_px_zy2
    F_alpha_px_zy_sum = np.sum(F_alpha_px_zy, axis=1)
    F_pzj_x = F_alpha_px_zy / F_alpha_px_zy_sum[:, None]
    F_phi_pzj_x = F_phi_2_np * F_pzj_x
    U_pyx = np.sum(F_phi_pzj_x, axis=1)
    
    return (D_pyx, U_pyx)

    
    
def logistic_regression(psm_dict, psm_list):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    psm_list_selected = []
    for key, lPsm in psm_dict.iteritems():
        sys.stdout.write(key + "\t")
        psm_list_selected = lPsm
        data_list = []
        label_list = []
        unknown_list = []
        id_list = []
        for oPsm in lPsm:
            #feature_list = []
            #oPsm.set_feature(feature_list)
            if len(oPsm.feature_list) == 0:
                print 'check'
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
            if oPsm.TrainingLabel != LabelUnknown:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)

        data_np = np.array(data_list)
        sys.stdout.write(str(len(data_list)) + "\t")
        label_np = np.array(label_list)
        
        '''
        for i in range(len(feature_name_list)):
            del feature_selection_list[:]
            feature_selection_list.append(i)
            for j in range(len(feature_name_list)):
                if j != i:
                    feature_selection_list.append(j)
        '''    
        '''
        del feature_selection_list[:]
        for i in range(len(feature_name_list)):
            if i < 6 or i > 23 or ( i>=15 and i<=17):
                feature_selection_list.append(i)
        #feature_selection_list.extend([1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28])
        '''
        train_data_np = data_np[:, feature_selection_list]
    
    # # training
        class_dict = {0: 1, 1:1000}
        logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
        logreg.fit(train_data_np, label_np)
        predict_np = logreg.predict(train_data_np)
        #np.savetxt("xxx.txt", train_data_np)
        #return

    #show_TP_TN_FP_FN(label_np, predict_np)
    #print 'Actual number of iterations for all classes.'
        #print logreg.n_iter_
    # # test
        unknown_np = np.array(unknown_list)
        test_unknown_np = unknown_np[:, feature_selection_list]
        predict_np = logreg.predict_proba(test_unknown_np)

        for i in range(len(predict_np)):
            psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]

    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
        psm_filtered_list_local = show_Fdr(psm_list_selected, key)
        psm_filtered_list.extend(psm_filtered_list_local)
    
    #print 'Coefficient of the features in the decision function:'
    #print logreg.coef_
        idx = 0
        for i in range(len(feature_name_list)):
            if i in feature_selection_list:
                sys.stdout.write('%.3f' % logreg.coef_[0][idx])
                idx += 1
            sys.stdout.write('\t')
        sys.stdout.write('\n')
        #print str(logreg.intercept_)
    '''
    for x in logreg.coef_[0]:
        sys.stdout.write(str(x))
        sys.stdout.write('\t')
    '''
    return psm_filtered_list

def fdr_rank(psm_list, feature_id):
    list_sorted = sorted(psm_list, key=lambda x: (x.ML_feature[feature_id]) , reverse=True)
    num_Fwd = 0
    num_Rev = 0
    num_Shu = 0
    for oPsm in list_sorted:
        if oPsm.RealLabel == LabelFwd:
            num_Fwd += 1
        elif oPsm.RealLabel == LabelNegative:
            num_Rev += 1
        else:
            num_Shu += 1
        (_FDR_accept, FDR_value) = FDR_calculator(num_Shu, num_Fwd)
        oPsm.FDR_feature.append(FDR_value)

def re_rank(psm_list):
    psm_new_list = []
    psm_dict = {}
    for oPsm in psm_list:
        sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if sId in psm_dict:
            if oPsm.fPredictProbability > psm_dict[sId].fPredictProbability:
                psm_dict[sId] = oPsm
        else:
            psm_dict[sId] = oPsm
    
    for _key, value in psm_dict.iteritems():
        psm_new_list.append(value)
    
    return psm_new_list

def test_MoE_ensamble(psm_list):
    psm_filtered_list = []
    id_list = []
    D_list = []
    U_list = []
    D_id_list = []
    U_id_list = []
    D_phi_list = []
    U_phi_list = []     
    for oPsm in psm_list:
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel == LabelNegative:
            D_list.append(oPsm.feature_list)
            D_id_list.append(oPsm.iInnerId)
            D_phi_list.append(oPsm.ML_feature)
        else:
            U_list.append(oPsm.feature_list)
            U_id_list.append(oPsm.iInnerId)
            U_phi_list.append(oPsm.ML_feature)

    D_np = np.array(D_list, dtype=np.float64)
    U_np = np.array(U_list, dtype=np.float64)
    
    
    D_phi_np = np.array(D_phi_list, dtype=np.float128)
    U_phi_np = np.array(U_phi_list, dtype=np.float128)
    feature_list_local = [1, 2, 3, 5, 15, 16, 17, 24, 25, 26, 27, 28, 29]
    D_np = D_np[:, feature_list_local]
    U_np = U_np[:, feature_list_local]
    # (D_pyx, U_pyx) = test_MoE_semi(D_np, U_np, D_phi_np, U_phi_np)
    (D_pyx, U_pyx) = test_MoE_super(D_np, U_np, D_phi_np, U_phi_np)
    # D_pyx = D_phi_np[:,1]
    # U_pyx = U_phi_np[:,1]
    
    for i in range(len(D_id_list)):
        psm_list[D_id_list[i]].fPredictProbability = 1 - D_pyx[i]
    
    for i in range(len(U_id_list)):
        psm_list[U_id_list[i]].fPredictProbability = U_pyx[i]
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
    
    return psm_filtered_list

def test_LR_ensamble(psm_list):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []        
    for oPsm in psm_list:
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.ML_feature)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.ML_feature)
            label_list.append(oPsm.TrainingLabel)
                
    train_data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)

    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    test_unknown_np = np.array(unknown_list)
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)

    psm_filtered_list_local = show_Fdr(psm_list, None, 0.01)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    for i in range(4):
        sys.stdout.write('%.3f' % logreg.coef_[0][i])
        sys.stdout.write('\t')
    sys.stdout.write('\n')
        
    return psm_filtered_list

def test_voting(psm_list, fdr_given=None):
    # re_rank
    psm_new_list = []
    psm_dict = {}
    for oPsm in psm_list:
        oPsm.set_fdr_product()
        sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
        if sId in psm_dict:
            if oPsm.fdr_product < psm_dict[sId].fdr_product:
                psm_dict[sId] = oPsm
        else:
            psm_dict[sId] = oPsm
    
    for _key, value in psm_dict.iteritems():
        psm_new_list.append(value)
    
    num_threshold = 20
    vote_threshold = 4
    fdr_threshold_list = [0.01] * num_threshold
    for i in range(num_threshold):
        fdr_threshold_list[i] += i * 0.001
    
    num_fwd_list = [0] * num_threshold
    num_rev_list = [0] * num_threshold
    num_shu_list = [0] * num_threshold
    for oPsm in psm_new_list:
        for i in range(num_threshold):
            if sum(val <= fdr_threshold_list[i] for val in oPsm.FDR_feature) >= vote_threshold:
                if oPsm.RealLabel == LabelFwd:
                    num_fwd_list[i] += 1
                elif oPsm.RealLabel == LabelRev:
                    num_rev_list[i] += 1
                else:
                    num_shu_list[i] += 1
    print '\nFDR\t# FWD\t# REV\t#SHU'
    for i in range(num_threshold):
        print str(fdr_threshold_list[i]) + '\t' + str(num_fwd_list[i]) + '\t' + str(num_rev_list[i]) + '\t' + str(num_shu_list[i])
    
    return None

def logistic_regression_no_category(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    '''
    for key, lPsm in psm_dict.iteritems():
        #sys.stdout.write(key + "\t")
        for oPsm in lPsm:
            if len(oPsm.feature_list) == 0:
                print 'check'
            if oPsm.TrainingLabel == LabelFiltered:
                continue
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
            if oPsm.TrainingLabel != LabelUnknown:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)
    '''           
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
        
    train_data_np = data_np[:, feature_selection_list]
    '''
    if not psm_left_list is None:
        np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/train_data_np.txt", train_data_np)
    '''
    # # training
    # num_positive = float((label_np==LabelPositive).sum())
    # num_negative = float((label_np==LabelNegative).sum())
    # num_positive = 100.0
    # num_negative = 1.0
    # class_weight_dict = {0: (num_positive/(num_negative+num_positive)), 1:(num_negative/(num_negative+num_positive))}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        # psm_list[id_list[i]].ML_feature.append(predict_np[i, 1])
    # fdr_rank(psm_list, 0)
    # return None
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    '''    
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_list, fdr_f)
    '''    
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % logreg.coef_[0][idx])
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
    '''
    iC = 0
    for oPsm in psm_filtered_list:
        if oPsm.ScoreAgreement == 1 and oPsm.RealLabel == LabelFwd:
            iC += 1
    print str(iC)
    '''        
        
    return psm_filtered_list

import random

def test_train_bias(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    train_data_list = []
    train_label_list = []
    train_id_list = []
    test_data_list = []
    test_id_list = []
    # # psm with the same name
    psm_sa_1_dict = {}
    
    for oPsm in psm_list:
        if oPsm.ScoreAgreement <= 1:
            sId = oPsm.FileName + '_' + str(oPsm.ScanNumber)
            if sId in psm_sa_1_dict:
                if psm_sa_1_dict[sId] == True:
                    train_data_list.append(oPsm.feature_list)
                    if oPsm.RealLabel == LabelRev:
                        train_label_list.append(0)
                    else:
                        train_label_list.append(1)
                    train_id_list.append(oPsm.iInnerId)
                else:
                    test_data_list.append(oPsm.feature_list)
                    test_id_list.append(oPsm.iInnerId)
            else:
                f_rand = random.random()
                if f_rand >= 0.5:
                    psm_sa_1_dict[sId] = True
                    train_data_list.append(oPsm.feature_list)
                    if oPsm.RealLabel == LabelRev:
                        train_label_list.append(0)
                    else:
                        train_label_list.append(1)
                    train_id_list.append(oPsm.iInnerId)
                else:
                    psm_sa_1_dict[sId] = False
                    test_data_list.append(oPsm.feature_list)
                    test_id_list.append(oPsm.iInnerId)
        else:
            f_rand = random.random()
            if f_rand >= 0.5:
                train_data_list.append(oPsm.feature_list)
                if oPsm.RealLabel == LabelRev:
                    train_label_list.append(0)
                else:
                    train_label_list.append(1)
                train_id_list.append(oPsm.iInnerId)
            else:
                test_data_list.append(oPsm.feature_list)
                test_id_list.append(oPsm.iInnerId)
                
    train_data_np = np.array(train_data_list)
    sys.stdout.write(str(len(train_data_list)) + "\t")
    train_label_np = np.array(train_label_list)
        
    train_data_np = train_data_np[:, feature_selection_list]
   
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             # class_weight='balanced',
                                             # class_weight=class_weight_dict, 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, train_label_np)

    # # test
    test_data_np = np.array(test_data_list)
    test_unknown_np = test_data_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[test_id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    print 'testing results:'
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_new_list, fdr_f)
    
    for i in range(len(predict_np)):
        psm_list[test_id_list[i]].fPredictProbability = 0 
    predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[train_id_list[i]].fPredictProbability = predict_np[i, 1]
    print 'train results:'
    psm_new_list = re_rank(psm_list)
    for i in range(1, 11):
        fdr_f = 0.001 * i
        show_Fdr_varied(psm_new_list, fdr_f)       
        
    return None

from sklearn import svm

def svm_one_class_test(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = [] 
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        if oPsm.TrainingLabel == LabelFiltered:
            continue
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel == LabelNegative:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
    
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
    train_data_np = data_np[:, feature_selection_list]
    
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
    
    clf.fit(train_data_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = clf.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
        
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, None, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, None, fdr=fdr_given)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list

def logistic_regression_no_category_rt(psm_dict, psm_list, psm_neg_list, fdr_given=None, psm_left_list=None):
    # machine learning
    # # construct training data
    psm_filtered_list = []
    #psm_list_selected = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []
    for key, lPsm in psm_dict.iteritems():
        #sys.stdout.write(key + "\t")
        for oPsm in lPsm:
            if len(oPsm.feature_list) == 0:
                print 'check'
                
            if oPsm.RealLabel == LabelRev:
                data_list.append(oPsm.feature_list)
                label_list.append(oPsm.TrainingLabel)
            if oPsm.TrainingLabel == LabelFiltered:
                if oPsm.RealLabel != LabelRev:
                    data_list.append(oPsm.feature_list)
                    label_list.append(oPsm.TrainingLabel)
                continue
            unknown_list.append(oPsm.feature_list)
            id_list.append(oPsm.iInnerId)
           
    for oPsm in psm_neg_list:
        data_list.append(oPsm.feature_list)
        label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\t")
    label_np = np.array(label_list)
        
    train_data_np = data_np[:, feature_selection_list]
    '''
    if not psm_left_list is None:
        np.savetxt("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/train_data_np.txt", train_data_np)
    '''
    # # training
    #class_dict = {0: 1, 1:1000}
    logreg = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    logreg.fit(train_data_np, label_np)
    #predict_np = logreg.predict(train_data_np)

    # # test
    unknown_np = np.array(unknown_list)
    test_unknown_np = unknown_np[:, feature_selection_list]
    predict_np = logreg.predict_proba(test_unknown_np)
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]

    # np.savetxt("predict_probability.txt", predict_np)
        #print str(len(psm_list_selected))
    if not psm_left_list is None:
        psm_filtered_list_local = show_Fdr(psm_left_list, key, fdr=fdr_given)
    else:
        psm_filtered_list_local = show_Fdr(psm_list, key, fdr=fdr_given)
    psm_filtered_list.extend(psm_filtered_list_local)
    # sys.stdout.write('\n')
    # show_Fdr_category(psm_dict)
    #print 'Coefficient of the features in the decision function:'
    #print logreg.coef_
    idx = 0
    for i in range(len(feature_name_list)):
        if i in feature_selection_list:
            sys.stdout.write('%.3f' % logreg.coef_[0][idx])
            idx += 1
        sys.stdout.write('\t')
    sys.stdout.write('\n')
        
    return psm_filtered_list
    
def get_num_missed_cleavage_sites(sIdentifiedSeq, sResiduesBeforeCleavage):
    count = 0
    for a in sResiduesBeforeCleavage:
        count += sIdentifiedSeq[:-2].count(a)
    return count

# generate NSM (# sibling modification), NSI (# sibling ions, charge), NSP (# sibling peptides), NMC (# missed cleavage sites)
def generate_Prophet_features(lPsm, config_dict):
    # peptide with PTM dictionary is for NSI
    peptide_with_modification_dict = {}
    # peptide without PTM dictionary is for NSM
    peptide_dict = {}
    peptide_protein_dict = {}
    for oPsm in lPsm:
        oPsm.NMC = get_num_missed_cleavage_sites(oPsm.OriginalPeptide, config_dict[cleave_after_residues_str])
        if oPsm.IdentifiedPeptide in peptide_with_modification_dict:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] += 1
        else:
            peptide_with_modification_dict[oPsm.IdentifiedPeptide] = 1
        if oPsm.OriginalPeptide in peptide_protein_dict:
            pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
            for protein in oPsm.protein_list:
                if not protein in pro_list:
                    pro_list.append(protein)
        else:
            pro_list = []
            pro_list.extend(oPsm.protein_list)
            peptide_protein_dict[oPsm.OriginalPeptide] = pro_list
            

    pattern = re.compile('[^\w\[\]]')
    for key, _value in peptide_with_modification_dict.iteritems():
        peptide_str = pattern.sub('', key)
        if peptide_str in peptide_dict:
            peptide_dict[peptide_str] += 1
        else:
            peptide_dict[peptide_str] = 1
    
    # # sibling peptides
    pro_unique_dict = {}
    pro_shared_dict = {}
    changed_flag = False
    for oPsm in lPsm:
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        changed_flag = False
        for protein in pro_list:
            if not protein in oPsm.protein_list:
                oPsm.protein_list.append(protein)
                changed_flag = True
        if changed_flag:
            oPsm.set_protein_names()
            oPsm.RealLabel = protein_type(oPsm.ProteinNames)
        if len(pro_list) > 1:
            for protein in pro_list:
                if protein in pro_shared_dict: 
                    pro_shared_dict[protein] += 1
                else:
                    pro_shared_dict[protein] = 1
        else:
            if pro_list[0] in pro_unique_dict:
                pro_unique_dict[pro_list[0]] += 1
            else:
                pro_unique_dict[pro_list[0]] = 1
    
    # collect features
    for oPsm in lPsm:
        oPsm.NSI = peptide_with_modification_dict[oPsm.IdentifiedPeptide] - 1
        oPsm.NSM = peptide_dict[oPsm.OriginalPeptide] - 1
        pro_list = peptide_protein_dict[oPsm.OriginalPeptide]
        for protein in pro_list:
            if protein in pro_unique_dict:
                oPsm.NSP += pro_unique_dict[protein]
                if len(pro_list) == 1:
                    oPsm.NSP -= 1
            if protein in pro_shared_dict:
                oPsm.NSPS += pro_shared_dict[protein]
                if len(pro_list) > 1:
                    oPsm.NSPS -= 1
                
# # Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)

# # Check file exist
def check_file_exist(filename):

    try:
        with open(filename) as _f: pass
    except IOError as _e:
        print >> sys.stderr, '\nCannot open', filename
        die("Program exit!")

# defaul value
decoy_prefix = 'Rev_'
min_peptide_per_protein = 2
min_unique_peptide_per_protein = 1
remove_decoy_identification = 'No'

pep_iden_str = '[Peptide_Identification]'
fasta_database_str = 'FASTA_Database'
pro_iden_str = '[Protein_Identification]'
decoy_prefix_str = 'Decoy_Prefix'
min_peptide_per_protein_str = 'Min_Peptide_Per_Protein'
min_unique_peptide_per_protein_str = 'Min_Unique_Peptide_Per_Protein'
remove_decoy_identification_str = 'Remove_Decoy_Identification'
cleave_after_residues_str = 'Cleave_After_Residues'

## Parse config file
def parse_config(config_filename):

    # Save config values to dictionary
    config_dict = {}    # initialize dictionay

    # Call Yinfeng's parseconfig.py module
    check_file_exist(config_filename)
    # Save all config values to dictionary
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
        elif key == (pep_iden_str + cleave_after_residues_str):
            config_dict[cleave_after_residues_str] = value
        else:
            continue

    # return config dictionary
    return config_dict


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
        
    def add(self, oPsm):
        self.SpectralCount += 1
        if self.BestScore < oPsm.lfScores[2]:
            self.BestScore = oPsm.lfScores[2]
        self.PSMs.append(oPsm.FileName+'['+str(oPsm.ScanNumber) +']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        
        
    def set(self, oPsm):
        self.IdentifiedPeptide = oPsm.IdentifiedPeptide
        self.ParentCharge = oPsm.ParentCharge
        self.OriginalPeptide = oPsm.OriginalPeptide
        self.ProteinNames = oPsm.ProteinNames
        self.ProteinCount = len(oPsm.protein_list)
        self.SpectralCount = 1
        self.BestScore = oPsm.lfScores[2]
        self.PSMs.append(oPsm.FileName+'['+str(oPsm.ScanNumber) +']')
        self.ScanType.append(oPsm.ScanType)
        self.SearchName.append(oPsm.SearchName)
        if oPsm.RealLabel == LabelFwd:
            self.TargetMatch = 'T'
        else:
            self.TargetMatch = 'F'
        
    def __repr__(self):
        l = [self.IdentifiedPeptide,
             str(self.ParentCharge),
             self.OriginalPeptide,
             self.ProteinNames,
             str(self.ProteinCount),
             self.TargetMatch,
             str(self.SpectralCount),
             str(self.BestScore),
             ('{'+','.join(self.PSMs)+'}'),
             ('{'+','.join(self.ScanType)+'}'),
             ('{'+','.join(self.SearchName)+'}')]
        
        return '\t'.join(l) 


def generate_psm_pep_txt(input_file, out_folder, psm_filtered_list):
    base_out_filename = input_file.split('/')[-1]
    base_out = out_folder + base_out_filename
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
        for oPsm in psm_filtered_list:
            if oPsm.RealLabel == LabelRev:
                continue 
            oPsm.clean_protein_name()
            fw.write(oPsm.FileName)
            fw.write('\t')
            fw.write(str(oPsm.ScanNumber))
            fw.write('\t')
            fw.write(str(oPsm.ParentCharge))
            fw.write('\t')
            fw.write('%.3f' % oPsm.MeasuredParentMass)
            fw.write('\t')
            fw.write('%.3f' % oPsm.CalculatedParentMass)
            fw.write('\t')
            fw.write('%.3f' % (oPsm.MeasuredParentMass - oPsm.CalculatedParentMass))
            fw.write('\t')
            fw.write('%.3f' % ((oPsm.MeasuredParentMass - oPsm.CalculatedParentMass)/oPsm.CalculatedParentMass))
            fw.write('\t')
            fw.write(oPsm.ScanType)
            fw.write('\t')
            fw.write(oPsm.SearchName)
            fw.write('\t')
            fw.write('Sipros10')
            fw.write('\t')
            fw.write('NA')
            fw.write('\t')
            fw.write('NA')
            fw.write('\t')
            fw.write(oPsm.DeltaP)
            fw.write('\t')
            fw.write(oPsm.IdentifiedPeptide)
            fw.write('\t')
            fw.write(oPsm.OriginalPeptide)
            fw.write('\t')
            fw.write(oPsm.ProteinNames)
            fw.write('\t')
            fw.write(str(len(oPsm.protein_list)))
            fw.write('\t')
            if oPsm.RealLabel == LabelFwd:
                fw.write('T')
            else:
                fw.write('F')
            fw.write('\n')
            
    # pep_sub_dict for preparing pep_out
    pep_sub_dict = {}    # initialize dict of list
    for oPsm in psm_filtered_list:
        if oPsm.RealLabel == LabelRev:
            continue
        pep_ID = oPsm.IdentifiedPeptide + '_+_' + str(oPsm.ParentCharge)
        if pep_ID in pep_sub_dict:
            pep_sub_dict[pep_ID].add(oPsm)
        else:
            oPeptide = Peptide()
            oPeptide.set(oPsm)
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
            fw.write(repr(oPeptide))
            fw.write('\n')


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

PepOutFields = sipros_post_module.PepOutFields
# # need pep.txt file and pro.txt file, and big PSM table file
# # Spectral Count for original/identified peptide
# # Unique peptide counts and total peptide counts for protein
def feature_update(run_num_dict, pro_file, psm_tab, psm_tab_new):
    
    # save the pep file and pro file data to the defaultdict
    id_pep_data_dict = {}
    or_pep_data_dict = {}
    pro_data_dict = {}
    
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
        f.write('\t'.join(headline))
        f.write('\tpep_psm\t')
        f.write('pro_pep')
        f.write('\n')
        
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
                f.write('2')
            elif iNumTotalPep > 1:
                f.write('1')
            else:
                f.write('0')
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

## Glboal variables
pep_file_ext = '.pep.txt'
psm_file_ext = '.psm.txt'
read_fasta_file = sipros_peptides_assembling.read_fasta_file
get_file_list_with_ext = sipros_post_module.get_file_list_with_ext
get_base_out = sipros_post_module.get_base_out
read_run_files = sipros_peptides_assembling.read_run_files
greedy_alg = sipros_peptides_assembling.greedy_alg
report_output = sipros_peptides_assembling.report_output

def assembly(output_folder, config_dict, psm_tab_file):
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
    return None
    # Feature Update
    sys.stderr.write('[Step 6] Feature update: Running -> ')
    base_out_filename = psm_tab_file.split('/')[-1]
    psm_tab_new = output_folder + base_out_filename
    feature_update(run_num_dict, base_out+'.pro.txt', psm_tab_file, psm_tab_new)
    return psm_tab_new

def show_measured_predicted_rt(psm_list, filename_prefix):
    fw_fwr = open(filename_prefix+"_fwr.txt", 'w')
    fw_fwr.write("measuread\tpredicted\n")
    fw_shu = open(filename_prefix+"_shu.txt", 'w')
    fw_shu.write("measuread\tpredicted\n")
    
    for oPsm in psm_list:
        if oPsm.RealLabel == LabelFwd:
            fw_fwr.write(str(oPsm.fRtMeasured))
            fw_fwr.write('\t')
            fw_fwr.write(str(oPsm.fRtPredict))
            fw_fwr.write('\n')
        else:
            fw_shu.write(str(oPsm.fRtMeasured))
            fw_shu.write('\t')
            fw_shu.write(str(oPsm.fRtPredict))
            fw_shu.write('\n')
    
    fw_fwr.close()
    fw_shu.close()
    
def stacking(block_num_int, psm_list):
    # # construct training data
    psm_filtered_list = []
    data_list = []
    label_list = []
    unknown_list = []
    id_list = []     
    for oPsm in psm_list:
        if len(oPsm.feature_list) == 0:
            print 'check'
        unknown_list.append(oPsm.feature_list)
        id_list.append(oPsm.iInnerId)
        if oPsm.TrainingLabel != LabelUnknown:
            data_list.append(oPsm.feature_list)
            label_list.append(oPsm.TrainingLabel)
                
    data_np = np.array(data_list)
    sys.stdout.write(str(len(data_list)) + "\n")
    label_np = np.array(label_list)
    
    # split the data into B blocks
    train_data_block = [[], []]
    train_label_block = [[], []]
    block_idx_list = np.random.randint(0, block_num_int, size=len(data_np))
    data_idx_int = 0
    for c in block_idx_list:
        if c == 0:
            train_data_block[0].append(data_np[data_idx_int])
            train_label_block[0].append(label_np[data_idx_int])
        else:
            train_data_block[1].append(data_np[data_idx_int])
            train_label_block[1].append(label_np[data_idx_int])
        data_idx_int += 1
        
    # features and models
    feature_matrix = [[1, 2, 3, 15, 16, 17],    # score section
                      [5, 29],                  # mass section
                      [24],                     # digestion section
                      [25, 26, 27, 28],         # pep pro support section
                      [30]]                     # PTM section
    all_feature_list = [1, 2, 3, 5, 15, 16, 17, 24, 25, 26, 27, 28, 29, 30]
    tier_1 = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        # tier_1[i].fit(np.array(train_data_block[1])[:, feature_matrix[i]], np.array(train_label_block[1]))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
        predict_np = tier_1[i].predict(data_np[:, feature_matrix[i]])
        # show_TP_TN_FP_FN(label_np, predict_np)
    
    
    
    
    tier_2 = linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1)
    # # train tier_2
    train_tier_2_np = []
    for i in range(len(feature_matrix)):
        predict_np = tier_1[i].predict_proba(data_np[:, feature_matrix[i]])
        train_tier_2_np.append(predict_np[:, 1])
    train_tier_2_np = np.transpose(np.array(train_tier_2_np))
    tier_2.fit(train_tier_2_np, label_np)
    '''
    
    # # re-train tier_1
    tier_1 = []
    for i in range(len(feature_matrix)):
        tier_1.append(linear_model.LogisticRegression(penalty='l2', 
                                             dual=False, 
                                             C=1, 
                                             fit_intercept=True, 
                                             intercept_scaling=1, 
                                             class_weight='balanced', 
                                             random_state=None, 
                                             solver='liblinear', 
                                             max_iter=100, 
                                             multi_class='ovr', 
                                             verbose=0, 
                                             warm_start=False, 
                                             n_jobs=-1))
        tier_1[i].fit(data_np[:, feature_matrix[i]], label_np)
    '''
    # # testing
    unknown_np = np.array(unknown_list)
    test_tier_2_np = []
    for i in range(len(feature_matrix)):
        predict_np = tier_1[i].predict_proba(unknown_np[:, feature_matrix[i]])
        test_tier_2_np.append(predict_np[:, 1])
    test_tier_2_np = np.transpose(np.array(test_tier_2_np))
    predict_np = tier_2.predict_proba(test_tier_2_np)
    
    #predict_np = logreg.predict_proba(train_data_np)
    for i in range(len(predict_np)):
        psm_list[id_list[i]].fPredictProbability = predict_np[i, 1]
    
    
    psm_new_list = re_rank(psm_list)
    del psm_list[:]
    psm_list.extend(psm_new_list)
    psm_filtered_list_local = show_Fdr(psm_list, None, fdr=0.01)
    # show_Fdr_charge(psm_list)
    psm_filtered_list.extend(psm_filtered_list_local)
    sys.stdout.write('\n')
        
    return psm_filtered_list

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse options
    (input_file, config_file, output_folder, negative_file) = parse_options(argv)
    
    # get the configuration parameters
    config_dict = parse_config(config_file)
    min_pep_shared = config_dict[min_peptide_per_protein_str]
    min_pep_unique = config_dict[min_unique_peptide_per_protein_str]
    psm_dict = {}
    psm_neg_list = []
    # get the categorized data
    (psm_dict, psm_neg_list) = categorize(input_file, negative_file)
    psm_list = []
    for key, lPsm in psm_dict.iteritems():
        psm_list.extend(lPsm)
    
    i = 0
    for oPsm in psm_list:
        oPsm.iInnerId = i
        i += 1
    # generate features
    generate_Prophet_features(psm_list, config_dict)
    generate_Prophet_features(psm_neg_list, config_dict)
    
    # q product filtering
    count_list = [0, 0, 0]
    for key, lPsm in psm_dict.iteritems():
        # print key
        # final_cutoff_score = get_cutoff_q_rank_product(0.01, lPsm)
        final_cutoff_score = 0.0
        mark_training_label(lPsm, final_cutoff_score, count_list, sKey = None)

    print 'Before Machine Learning:\n\t# FWD\t# REV\t# SHU'
    print '\t%d\t%d\t%d' % (count_list[0], count_list[1], count_list[2])
    
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_list()
        
    for oPsm in psm_neg_list:
        oPsm.get_feature_list()
    '''
    # ensamble learning
    psm_filtered_list = test_stacking_MOE(psm_list)
    # stacking(2, psm_list)
    print 'Done.'
    return
    '''
    # machine learning
    #test_random_forest(psm_dict, psm_list)
    del feature_selection_list[:]
    feature_selection_list.extend([1, 2, 3, 5, 15, 16, 17, 24, 25, 26, 27, 28, 29]) #
    # psm_filtered_list = logistic_regression(psm_dict, psm_list)
    psm_filtered_list = logistic_regression_no_category(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # train bias test
    # psm_filtered_list = test_train_bias(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # neural network
    # psm_filtered_list = svm_one_class_test(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # random forest
    # psm_filtered_list = test_random_forest(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # naive bayes
    # psm_filtered_list = test_naive_bayes(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # SVM
    # psm_filtered_list = test_SVM(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # ada boost
    # psm_filtered_list = test_AdaBoost(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # deep learning
    # psm_filtered_list = test_DeepLearning(psm_dict, psm_list, psm_neg_list, 0.01, None)
    # ensemble voting
    # psm_filtered_list = test_voting(psm_list, 0.01)
    # ensemble LR
    # psm_filtered_list = test_LR_ensamble(psm_list)
    # generate_psm_pep_txt(input_file, output_folder, psm_filtered_list)
    # ensemble MoE
    # psm_filtered_list = test_MoE_ensamble(psm_list)
    print 'Done.'
    return
    # RT data
    # get_RT_data(psm_filtered_list, output_folder, psm_list)
    #get_RT_elude_data(psm_filtered_list, output_folder, psm_list)
    # (psm_selected_list, psm_left_list) = filter_Fdr(psm_list, 0.01)
    # (psm_selected_list, psm_left_list) = filter_Fdr2(psm_list, 0.01)
    # show_measured_predicted_rt(psm_selected_list, '/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/fdrless0.01')
    # show_measured_predicted_rt(psm_left_list, '/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/check/fdrlarger0.01')
    # psm_filtered_list = logistic_regression_no_category_rt(psm_dict, psm_list, psm_neg_list, 0.01, psm_left_list)
    
    # return
    
    # clean_folder(output_folder)
    # generate_psm_pep_txt(input_file, output_folder, psm_filtered_list)
    
    config_dict[min_peptide_per_protein_str] = 1
    config_dict[min_unique_peptide_per_protein_str] = 1
    psm_table_new_file = assembly(output_folder, config_dict, input_file)
    print 'Done.'
    return
    # # second round #########################
    # get the categorized data
    psm_dict = categorize(psm_table_new_file)
    psm_list = []
    for key, lPsm in psm_dict.iteritems():
        psm_list.extend(lPsm)
    
    i = 0
    for oPsm in psm_list:
        oPsm.iInnerId = i
        i += 1
    # generate features
    generate_Prophet_features(psm_list, config_dict)
    
    # q product filtering
    count_list = [0, 0, 0]
    for key, lPsm in psm_dict.iteritems():
        print key
        final_cutoff_score = get_cutoff_q_rank_product(0.01, lPsm)
        mark_training_label(lPsm, final_cutoff_score, count_list, key)

    print 'Before Machine Learning:\n\t# FWD\t# REV\t# SHU'
    print '\t%d\t%d\t%d' % (count_list[0], count_list[1], count_list[2])
    
    # set feature all PSMs
    for oPsm in psm_list:
        oPsm.get_feature_list()

    # machine learning
    #test_random_forest(psm_dict, psm_list)
    del feature_selection_list[:]
    feature_selection_list.extend([0, 1, 2, 3, 4, 5, 15, 16, 17, 24, 25, 26, 27, 28])
    psm_filtered_list = logistic_regression(psm_dict, psm_list)
    
    print "Done."
    
def main2(argv=None):
    return
    if argv is None:
        argv = sys.argv

    # parse options
    (input_file, config_file, output_folder, negative_file) = parse_options(argv)
    
    # get the configuration parameters
    config_dict = parse_config(config_file)
    config_dict[min_peptide_per_protein_str] = 1
    config_dict[min_unique_peptide_per_protein_str] = 1
    assembly(output_folder, config_dict, input_file)
    print 'Done.'
    return
    

if __name__ == '__main__':
    main()
    sys.exit(main2())
