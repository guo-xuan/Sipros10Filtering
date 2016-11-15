'''
Created on Aug 4, 2016

@author: xgo
'''
import getopt, sys
import numpy as np
from scipy.stats import rankdata
import sipros_post_module
from audioop import reverse

# # Version control
def get_version():
    return "1.0.1 (Alpha)"

# # Help message
help_message = '''
Usage:
    python sipros_psm_rank.py [options]

Inputs:
    tab file
    output file

Options:
    -h/--help
    -v/--version
    -i/--input-file
    -o/--output-file

Outputs:
    output PSM table + rank info
'''


# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVi:c:o:",
                                    ["help",
                                     "version",
                                     "input-file",
                                     "output-file"])

    # Default working dir and config file
    input_file = ""
    output_file = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print help_message
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print "sipros_peptides_filtering.py V%s" % (get_version())
            sys.exit(0)
        if option in ("-i", "--input-file"):
            input_file = value
        if option in ("-o", "--output-file"):
            output_file = value

    if input_file == "" or output_file == "":
        print help_message
        sys.exit(0)

    return (input_file, output_file)

class PSM:

    def __init__(self):
        self.sLine = ''
        self.bIsTarget = False
        self.fRank = 0.0
        self.lScores = []
        self.lRanks = []

    def readInfo(self, sLine):
        self.sLine = sLine.strip()
        asWords = sLine.strip().split('\t')
        self.check_decoy_match(asWords[12], 'Rev_')
        self.lScores.append(float(asWords[9]))
        self.lScores.append(float(asWords[10]))
        self.lScores.append(float(asWords[11]))
        return (float(asWords[9]), float(asWords[10]), float(asWords[11]))

    def check_decoy_match(self, protein_sequence, decoy_prefix):
        sProteins = protein_sequence.replace('{', '')
        sProteins = sProteins.replace('}', '')
        asProteins = sProteins.split(',')
        for sProtein in asProteins:
            if not sProtein.startswith(decoy_prefix):
                self.bIsTarget = True
                return
        self.bIsTarget = False

def read_tab_file(input_file):
    lPsm = []
    lMvh = []
    lXcorr = []
    lWdp = []
    with open(input_file, 'r') as f:
        # skip the header
        _sLine = f.readline()
        # read the PSM
        for sLine in f:
            oPsm = PSM()
            (Mvh, Xcorr, Wdp) = oPsm.readInfo(sLine)
            lPsm.append(oPsm)
            lMvh.append(-Mvh)
            lXcorr.append(-Xcorr)
            lWdp.append(-Wdp)

    return (lPsm, lMvh, lXcorr, lWdp)

def read_tab_file_q_value(input_file):
    lPsm = []

    with open(input_file, 'r') as f:
        # skip the header
        _sLine = f.readline()
        # read the PSM
        for sLine in f:
            oPsm = PSM()
            oPsm.readInfo(sLine)
            lPsm.append(oPsm)

    return (lPsm)

# # Division error handling
divide = sipros_post_module.divide
FDR_parameter = 1.0

# # FDR calculator
def FDR_calculator(FP, TP):
    FDR_numerator = float(FP) * float(FDR_parameter)
    FDR_denominator = float(FP) + float(TP)

    if  FDR_denominator == 0:
        FDR_value = 1.0
    else:
        FDR_value = divide(FDR_numerator, FDR_denominator)

    return float(FDR_value)


def q_value_product(lPsm):

    for i in range(3):
        newlist = sorted(lPsm, key=lambda x: x.lScores[i], reverse=True)
        T_num = 0
        F_num = 0
        for j in range(len(newlist)):
            if newlist[j].bIsTarget:
                T_num += 1
            else:
                F_num += 1
            newlist[j].lRanks.append(FDR_calculator(F_num, T_num))
        fSmallestQ = 1
        for j in range(len(newlist) - 1, -1, -1):
            if fSmallestQ > newlist[j].lRanks[i]:
                fSmallestQ = newlist[j].lRanks[i]
            if newlist[j].lRanks[i] > fSmallestQ:
                newlist[j].lRanks[i] = fSmallestQ

    for oPsm in lPsm:
        #fTemp = oPsm.lRanks[0] * oPsm.lRanks[1] * oPsm.lRanks[2]
        fTemp = (1 - ((1 - oPsm.lRanks[0]) * (1 - oPsm.lRanks[1]) * (1 - oPsm.lRanks[2])))
        oPsm.fRank = np.power(float(fTemp), 1.0 / 3.0)

def rank_product(lPsm, lMvh, lXcorr, lWdp):
    lMvhRank = np.array(lMvh)
    lMvhRank = rankdata(lMvhRank, method='max')

    lXcorrRank = np.array(lXcorr)
    lXcorrRank = rankdata(lXcorrRank, method='max')

    lWdpRank = np.array(lWdp)
    lWdpRank = rankdata(lWdpRank, method='max')

    Psm_rank_np = np.multiply(lMvhRank, lXcorrRank)
    Psm_rank_np = np.multiply(Psm_rank_np, lWdpRank)

    for i in range(len(Psm_rank_np)):
        lPsm[i].fRank = np.power(Psm_rank_np[i], 1.0 / 3.0)

def output(lPsm, output_file):
    with open(output_file, 'w') as f:
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
        f.write('RankProduct')
        f.write('\n')
        for oPsm in lPsm:
            f.write(oPsm.sLine)
            f.write('\t')
            f.write(str(oPsm.fRank))
            f.write('\n')

def main(argv=None):

    if argv is None:
        argv = sys.argv

    # parse options
    (psm_tab_file, output_file) = parse_options(argv)

    # (lPsm, lMvh, lXcorr, lWdp) = read_tab_file(psm_tab_file)

    # rank_product(lPsm, lMvh, lXcorr, lWdp)

    # output(lPsm, output_file)


    (lPsm) = read_tab_file_q_value(psm_tab_file)

    q_value_product(lPsm)

    output(lPsm, output_file)

    print "Done."

if __name__ == '__main__':
    sys.exit(main())
