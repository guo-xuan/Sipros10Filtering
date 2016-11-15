'''
Created on Sep 21, 2016

@author: xgo
'''

import sys

def main(argv=None):
    base_num = 275391
    ms2_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/qranker/comet/0.09Da/ms2/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_11.ms2"
    sqt_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/qranker/comet/0.09Da/sqt/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_11.sqt"
    
    ms2_new_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/qranker/comet/0.09Da/ms2/reindex/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_11.ms2"
    sqt_new_file = '/media/xgo/Seagate/Proteomics/Experiments/Angelo/qranker/comet/0.09Da/sqt/reindex/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_11.sqt'
    
    scan_num = 1
    scan_dict_1 = {}
    scan_dict_2 = {}
    
    with open(ms2_file, 'r') as fr:
        with open(ms2_new_file, 'w') as fw:
            for line_str in fr:
                if line_str[0] == 'H':
                    continue
                if line_str[0] == 'S':
                    words = line_str.split()
                    scan_dict_1[int(words[1])] = scan_num + base_num
                    scan_dict_2[int(words[2])] = scan_num + base_num
                    fw.write(words[0])
                    fw.write('\t')
                    fw.write(str(scan_num+base_num))
                    fw.write('\t')
                    fw.write(str(scan_num+base_num))
                    fw.write('\t')
                    fw.write('\t'.join(words[3:]))
                    fw.write('\n')
                    scan_num += 1
                else:
                    fw.write(line_str)
    
    print scan_num + base_num
    
    with open(sqt_file, 'r') as fr:
        with open(sqt_new_file, 'w') as fw:
            for line_str in fr:
                if line_str[0] == 'H':
                    continue
                if line_str[0] == 'S':
                    words = line_str.split()
                    if int(words[1]) not in scan_dict_1:
                        print "error 1"
                    if int(words[2]) not in scan_dict_2:
                        print "error 2"
                    fw.write(words[0])
                    fw.write('\t')
                    fw.write(str(scan_dict_1[int(words[1])]))
                    fw.write('\t')
                    fw.write(str(scan_dict_2[int(words[2])]))
                    fw.write('\t')
                    fw.write('\t'.join(words[3:]))
                    fw.write('\n')
                else:
                    fw.write(line_str)
    
    print 'Done.'

if __name__ == '__main__':
    sys.exit(main())