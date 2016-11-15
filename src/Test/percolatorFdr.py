'''
Created on Sep 21, 2016

@author: xgo
'''

import sys, re

import sipros_post_module

# # Division error handling
divide = sipros_post_module.divide
FDR_parameter = 1.0

rev_str = 'Sh1_'
shu_str = 'Sh2_'

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

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    asProteins = protein_sequence.split(',')
    if lProtein != None:
        del lProtein[:]
        lProtein.extend(asProteins[:])
    for sProtein in asProteins:
        if not (sProtein.startswith(rev_str) or sProtein.startswith(shu_str)):
            return 1
    for sProtein in asProteins:
        if sProtein.startswith(rev_str):
            return 2
        elif sProtein.startswith(shu_str):
            return 3
        
def get_peptide(pep_str):
    pos_start = pep_str.find('.') + 1
    pos_end = pep_str.rfind('.')
    return pep_str[pos_start:pos_end]

def try_fdr(fdr_f):
    
    input_file = '/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Comet/09272013/pin/all.txt'
    
    psm_list = []
    with open(input_file, 'r') as f:
        f.next()
        for line in f:
            word_list = line.split()
            score = float(word_list[1])
            qvalue = float(word_list[2])
            protein = protein_type(word_list[5])
            psm_list.append((score, qvalue, protein))
    
    psm_list_sorted = sorted(psm_list, key = lambda psm: psm[0], reverse = True)
    
    num_shu = 0
    num_rev = 0
    num_fwr = 0
    best_fwr = 0
    best_rev = 0
    best_shu = 0
    for (_a, _c, b) in psm_list_sorted:
        if b == 1:
            num_fwr += 1
        elif b == 2:
            num_rev += 1
        elif b == 3:
            num_shu += 1
        else:
            print "Error."
        (FDR_accept, FDR_value) = FDR_calculator(num_shu, num_fwr)
        if (FDR_accept is True) and (FDR_value <= fdr_f) and ((num_fwr + num_shu) > (best_fwr + best_shu)) :
            best_fwr = num_fwr
            best_shu = num_shu
            best_rev = num_rev
        
    # print "# FWD\t# REV\t# SHU"
    print "%d\t%d\t%d" %(best_fwr, best_rev, best_shu)

    # print("Done.")

def main(argv=None):
    
    input_file = '/media/xgo/Seagate/Proteomics/Experiments/Benchmark/Comet/08202013/pin/all.txt'
    
    psm_list = []
    with open(input_file, 'r') as f:
        f.next()
        for line in f:
            word_list = line.split()
            score = float(word_list[1])
            qvalue = float(word_list[2])
            protein = protein_type(word_list[5])
            psm_list.append((score, qvalue, protein))
    
    psm_list_sorted = sorted(psm_list, key = lambda psm: psm[0], reverse = True)
    
    num_shu = 0
    num_rev = 0
    num_fwr = 0
    best_fwr = 0
    best_rev = 0
    best_shu = 0
    for (_a, _c, b) in psm_list_sorted:
        if b == 1:
            num_fwr += 1
        elif b == 2:
            num_rev += 1
        elif b == 3:
            num_shu += 1
        else:
            print "Error."
        (FDR_accept, FDR_value) = FDR_calculator(num_shu, num_fwr)
        if (FDR_accept is True) and (FDR_value <= 0.01) and ((num_fwr + num_shu) > (best_fwr + best_shu)) :
            best_fwr = num_fwr
            best_shu = num_shu
            best_rev = num_rev
        
    print "# FWD\t# REV\t# SHU"
    print "%d\t%d\t%d" %(best_fwr, best_rev, best_shu)

    print("Done.")
    
def main2(argv=None):
    print "#Total\t# FWD\t# REV\t# SHU"
    for i in range(1, 11):
        fdr = 0.001*i
        try_fdr(fdr)
    print 'Done.'

if __name__ == '__main__':
    sys.exit(main2())