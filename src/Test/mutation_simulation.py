'''
Created on Oct 24, 2016

@author: xgo
'''

import sys, os, getopt, random

all_aa_list = ['J', 'I', 'L', 'A', 'S', 'G', 'V', 'E', 'K', 'T', 'D', 'R', 'P', 'N', 'F', 'Q', 'Y', 'M', 'H', 'C', 'W']
mark_beg_str = "(((("
mark_end_str = "))))"

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hi:d:o:",
                                    ["help",
                                     "input",
                                     "database"
                                     "output"])

    # Default working dir and config file
    pep_txt_file = ""
    new_file = ""
    protein_db_file = ""

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            sys.exit(0)
        if option in ("-i", "--input"):
            pep_txt_file = value
        if option in ("-o", "--output"):
            new_file = value
        if option in ("-d", "--database"):
            protein_db_file = value

    return (pep_txt_file, protein_db_file, new_file)

def parse_peptide_protein(pep_txt_file):

    pep_pro_dict = {}
    pro_pep_dict = {}
    with open(pep_txt_file, 'r') as fr:
        fr.next()
        for line_str in fr:
            line_array = line_str.split()
            if line_array[5] == 'T':
                pep_pro_dict[line_array[2]] = line_array[3][1:-1]

    # only choose the 10%, one peptide one protein
    for pep, pro in pep_pro_dict.iteritems():
        if random.randint(0, 9) == 0:
            line_array = pro.split(',')
            for protein in line_array:
                    pro_pep_dict[protein] = pep[1:-1]

    return pro_pep_dict

def mutate_peptide_db(pro_pep_dict, protein_db_file, new_file):
    seq_str = ""
    protein_id_str = ""
    pos_pep_beg_int = 0
    pep_str = ""
    mutation_pos_relative_int = 0
    mutation_pos_absolute_int = 0
    subsitute_char = ''
    with open(protein_db_file, 'r') as fr:
        with open(new_file, 'w') as fw:
            for line_str in fr:
                if line_str[0] == '>':
                    if seq_str != "" and protein_id_str != "":
                        if protein_id_str in pro_pep_dict:
                            pep_str = pro_pep_dict[protein_id_str]
                            pos_pep_beg_int = seq_str.find(pep_str)
                            mutation_pos_relative_int = random.randint(0, len(pep_str) - 1)
                            subsitute_char = all_aa_list[random.randint(0, len(all_aa_list) - 1)]
                            while subsitute_char == pep_str[mutation_pos_relative_int]:
                                subsitute_char = all_aa_list[random.randint(0, len(all_aa_list) - 1)]
                            mutation_pos_absolute_int = pos_pep_beg_int + mutation_pos_relative_int
                            protein_id_str += mark_beg_str + str(mutation_pos_absolute_int) + pep_str[mutation_pos_relative_int] + mark_end_str
                            new_list = list(seq_str)
                            new_list[mutation_pos_absolute_int] = subsitute_char
                            seq_str = ''.join(new_list)
                            fw.write('>')
                            fw.write(protein_id_str)
                            fw.write('\n')
                            fw.write(seq_str)
                            fw.write('\n')
                        else:
                            fw.write('>')
                            fw.write(protein_id_str)
                            fw.write('\n')
                            fw.write(seq_str)
                            fw.write('\n')
                    if line_str.startswith('>Rev_') or line_str.startswith('>Dec_'):
                        protein_id_str = ""
                        protein_id_str = line_str[1:].split()[0]
                    else:
                        # get the protein id
                        protein_id_str = line_str[1:].split()[0]

                    seq_str = ""
                else:
                    seq_str += line_str.strip()

            if seq_str != "" and protein_id_str != "":
                if protein_id_str in pro_pep_dict:
                    pep_str = pro_pep_dict[protein_id_str]
                    pos_pep_beg_int = seq_str.find(pep_str)
                    mutation_pos_relative_int = random.randint(0, len(pep_str) - 1)
                    subsitute_char = all_aa_list[random.randint(0, len(all_aa_list) - 1)]
                    while subsitute_char == pep_str[mutation_pos_relative_int]:
                        subsitute_char = all_aa_list[random.randint(0, len(all_aa_list) - 1)]
                    mutation_pos_absolute_int = pos_pep_beg_int + mutation_pos_relative_int
                    protein_id_str += mark_beg_str + str(mutation_pos_absolute_int) + pep_str[mutation_pos_relative_int] + mark_end_str
                    new_list = list(seq_str)
                    new_list[mutation_pos_absolute_int] = subsitute_char
                    seq_str = ''.join(new_list)
                    fw.write('>')
                    fw.write(protein_id_str)
                    fw.write('\n')
                    fw.write(seq_str)
                    fw.write('\n')
                else:
                    fw.write('>')
                    fw.write(protein_id_str)
                    fw.write('\n')
                    fw.write(seq_str)
                    fw.write('\n')

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse options
    (pep_txt_file, protein_db_file, new_file) = parse_options(argv)
    # read the peptides and proteins from the psm.txt file
    pro_pep_dict = parse_peptide_protein(pep_txt_file)
    # read in the database and mark the mutations
    mutate_peptide_db(pro_pep_dict, protein_db_file, new_file)
    # done
    print 'Done.'

if __name__ == '__main__':
    sys.exit(main())
