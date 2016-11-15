'''
Created on Aug 30, 2016

@author: xgo
'''

import os, sys
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

def main():
    folder_str = '/home/xgo/Temp/Angelo/0.09Da/Myrimatch/AfterPeptideProphet/'
    file_list = []
    for file_name in os.listdir(folder_str):
        file_path_name = folder_str + file_name
        if file_path_name.endswith('xls'):
            file_list.append(file_path_name)
    
    probability_float = 0.0
    scan_id_str = ''
    protein_type_int = 0
    psm_list = []
    for file_str in file_list:
        psm_dict = {}
        with open(file_str, 'r') as f:
            f.next()
            for line_str in f:
                word_list = line_str.split()
                if any(x == '[unavailable]' for x in word_list) :
                    continue
                probability_float = float(word_list[0])
                scan_id_str = word_list[1][word_list[1].find('.')+1:word_list[1].rfind('.')]
                protein_type_int = protein_type(word_list[6])
                if scan_id_str in psm_dict:
                    if probability_float > psm_dict[scan_id_str][1]:
                        psm_dict[scan_id_str] = (protein_type_int, probability_float)
                else:
                    psm_dict[scan_id_str] = (protein_type_int, probability_float)
        
        print '# PSM:\t%i' % len(psm_dict)
        
        output_file_str = ''
        if output_file_str != '':
            file_writer = open(output_file_str, 'w')
        for _key, value in psm_dict.iteritems():
            psm_list.append(value)
            if output_file_str != '':
                file_writer.write(_key+'\n')
        if output_file_str != '':
            file_writer.close()
            
    psm_list_sorted = sorted(psm_list, key=lambda tup: tup[1], reverse = True)
    num_shu = 0
    num_rev = 0
    num_fwr = 0
    best_fwr = 0
    best_rev = 0
    best_shu = 0
    for (b, _c) in psm_list_sorted:
        if b == 1:
            num_fwr += 1
        elif b == 2:
            num_rev += 1
            print _c
            break
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
    print "%d\t%d\t%d" % (best_fwr, best_rev, best_shu)
    

if __name__ == '__main__':
    sys.exit(main())