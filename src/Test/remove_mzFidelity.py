'''
Created on Aug 29, 2016

@author: xgo
'''

import sys

def main():
    input_file_str = '/home/xgo/Temp/pepXML/MyriMatch/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_11.pepXML'
    output_file_str = '/home/xgo/Temp/pepXML/MyriMatch(MVH_only)/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_11.pep.xml'
    
    output = open(output_file_str, 'w')
    with open(input_file_str, 'r') as f:
        for line_str in f:
            new_line_str = line_str.lstrip()
            if not new_line_str.startswith('<search_score name="mzFidelity"'):
                output.write(line_str)
            else:
                pass
    output.close()
    
    print 'Done.'

if __name__ == '__main__':
    sys.exit(main())