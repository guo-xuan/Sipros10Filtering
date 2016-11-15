'''
Created on Aug 10, 2016

@author: xgo
'''

import sys, re

import sipros_post_module

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

# # Decoy Reverse Forward protein
def protein_type(protein_sequence, lProtein=None):
    asProteins = protein_sequence.split(',')
    if lProtein != None:
        del lProtein[:]
        lProtein.extend(asProteins[:])
    for sProtein in asProteins:
        if not (sProtein.startswith('Rev_') or sProtein.startswith('Dec_')):
            return 1
    for sProtein in asProteins:
        if sProtein.startswith('Rev_'):
            return 2
        elif sProtein.startswith('Dec_'):
            return 3
        
def get_peptide(pep_str):
    pos_start = pep_str.find('.') + 1
    pos_end = pep_str.rfind('.')
    return pep_str[pos_start:pos_end]

def main(argv=None):
    
    input_file = '/home/xgo/Temp/Angelo/0.09Da/Comet/Angelo_Comet_Best_NoMass_iProphet.xls'
    
    scan_id_str = ''
    probability_float = 0.0
    mass_diff = 0.0
    ppm = 0.0
    psm_dict = {}
    pep_str = ''
    regex = re.compile('[^a-zA-Z]')
    with open(input_file, 'r') as f:
        f.next()
        for line in f:
            word_list = line.split()
            scan_id_str = (word_list[1][word_list[1].rfind('_')+1:]).split('.')
            scan_id_str = scan_id_str[0] + '_' + scan_id_str[1]
            probability_float = float(word_list[3])
            protein_type_int = protein_type(word_list[6])
            if len(word_list) == 10:
                mass_diff = float(word_list[8])
                temp_mass_diff = mass_diff
                for i in range(-5, 6):
                    if abs(mass_diff - (i * 1.00727647)) < abs(temp_mass_diff):
                        temp_mass_diff = (mass_diff - (i * 1.00727647))
                mass_diff = temp_mass_diff
                ppm = (abs(mass_diff) / float(word_list[7]))*1000000
                pep_str = get_peptide(word_list[5])
                pep_str = regex.sub('', pep_str)
            if scan_id_str in psm_dict:
                if probability_float > psm_dict[scan_id_str][1]:
                    psm_dict[scan_id_str] = (protein_type_int, probability_float, ppm, pep_str)
            else:
                psm_dict[scan_id_str] = (protein_type_int, probability_float, ppm, pep_str)
    
    psm_list = []   
    
    for _key, value in psm_dict.iteritems():
        psm_list.append(value)
    
    psm_list_sorted = sorted(psm_list, key = lambda psm: psm[1], reverse = True)
    
    num_shu = 0
    num_rev = 0
    num_fwr = 0
    best_fwr = 0
    best_rev = 0
    best_shu = 0
    psm_filitered_list = []
    for (b, _c, d, e) in psm_list_sorted:
        if b == 1:
            num_fwr += 1
        elif b == 2:
            num_rev += 1
        elif b == 3:
            num_shu += 1
        else:
            print "Error."
        psm_filitered_list.append([_c, d, e])
        (FDR_accept, FDR_value) = FDR_calculator(num_shu, num_fwr)
        if (FDR_accept is True) and (FDR_value <= 0.01) and ((num_fwr + num_shu) > (best_fwr + best_shu)) :
            best_fwr = num_fwr
            best_shu = num_shu
            best_rev = num_rev
        
    print "# FWD\t# REV\t# SHU"
    print "%d\t%d\t%d" % (best_fwr, best_rev, best_shu)
    
    '''
    fw1 = open('10ppm.txt', 'w')
    fw2 = open('rest.txt', 'w')
    for psm in psm_filitered_list:
        if psm[1] <= 10:
            fw1.write(psm[2])
            fw1.write('\t')
            fw1.write(str(psm[1]))
            fw1.write('\n')
        else:
            fw2.write(psm[2])
            fw2.write('\t')
            fw2.write(str(psm[1]))
            fw2.write('\n')
    fw1.close()
    fw2.close()
    '''
    print("Done.")

if __name__ == '__main__':
    sys.exit(main())