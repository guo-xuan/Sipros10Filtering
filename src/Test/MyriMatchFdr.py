'''
Created on Aug 30, 2016

@author: xgo
'''
import sys

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

def main(argv=None):
    
    peptide_file_str = '/media/xgo/Seagate/Proteomics/Experiments/Angelo/Myrimatch/Peptide.tsv'
    peptide_protain_dict = {}
    
    with open(peptide_file_str, 'r') as f:
        f.next()
        for line_str in f:
            word_list = line_str.split('\t')
            peptide_str = (word_list[0].split(' '))[0]
            protein_str = word_list[5].strip()
            protein_type_int = protein_type(protein_str[1:-1])
            if peptide_str in peptide_protain_dict:
                if protein_type_int == 1:
                    peptide_protain_dict[peptide_str] = protein_type_int
            else:
                peptide_protain_dict[peptide_str] = protein_type_int
            
    psm_file_str = '/media/xgo/Seagate/Proteomics/Experiments/Angelo/Myrimatch/psm.tsv'
    
    psm_dict = {}
    
    file_id_str = ''
    scan_id_str = ''
    Q_value_float = 0.0
    with open(psm_file_str, 'r') as f:
        f.next()
        f.next()
        f.next()
        for line_str in f:
            line_str = line_str.strip()
            word_list = line_str.split()
            if word_list[0].startswith('Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_'):
                file_id_str = word_list[0][-2:]
            else:
                scan_id_str = file_id_str + '_' + word_list[0]
                peptide_str = word_list[-1]
                Q_value_float = float(word_list[-2])
                protein_type_int = peptide_protain_dict[peptide_str]
                if scan_id_str in psm_dict:
                    if protein_type_int == 1:
                        if Q_value_float < psm_dict[scan_id_str][1]:
                            psm_dict[scan_id_str] = (protein_type_int, Q_value_float)
                        elif psm_dict[scan_id_str][0] != 1:
                            psm_dict[scan_id_str] = (protein_type_int, Q_value_float)
                else:
                    psm_dict[scan_id_str] = (protein_type_int, Q_value_float)
    psm_list = []    
    for key, value in psm_dict.iteritems():
        psm_list.append((key, value[0], value[1]))
    
    psm_list_sorted = sorted(psm_list, key=lambda tup: tup[2])
    num_shu = 0
    num_rev = 0
    num_fwr = 0
    for (_a, b, _c) in psm_list_sorted:
        if b == 1:
            num_fwr += 1
        elif b == 2:
            num_rev += 1
        elif b == 3:
            num_shu += 1
        else:
            print "Error."
        
        if num_shu >= 618:
            break
        
    print "# FWD\t# REV\t# SHU"
    print "%d\t%d\t%d" %(num_fwr, num_rev, num_shu)
    

if __name__ == '__main__':
    sys.exit(main())