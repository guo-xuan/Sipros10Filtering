'''
Created on Sep 21, 2016

@author: xgo
'''

import sys, random

def get_cross_validate_data(sId):
    input_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Data/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_"+sId+".txt"
    output_folder = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Data/train_test/"
    
    shift_float = 0
    
    pep_list = []
    
    with open(input_file, 'r') as fr:
        for line_str in fr:
            words = line_str.split()
            time_float = float(words[1])
            # time_float += shift_float + 120
            # time_float %= 120 
            pep_list.append((words[0], str(time_float), words[2]))
    
    fw_train = open(output_folder+"train_"+sId+".txt", 'w')
    fw_test = open(output_folder+"test_"+sId+".txt", 'w')
    
    test_id = -1
    for pep in pep_list:
        if pep[2] == '2':
            fw_test.write(pep[0])
            fw_test.write('\n')
        else:
            test_id =  random.randint(0, 1)
            if test_id == 0:
                fw_train.write(pep[0])
                fw_train.write('\t')
                fw_train.write(pep[1])
                fw_train.write('\n')
            else:
                fw_test.write(pep[0])
                fw_test.write('\n')

    fw_train.close()
    fw_test.close()

def get_elude_RTpredict_data():
    elude_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Elude/result_1.txt"
    rt_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/RTmodel/result_1.csv"
    elude_rt_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/elude_RTpredict.txt"
    
    psm_dict = {}
    
    with open(elude_file, 'r') as fr:
        fr.next()
        fr.next()
        fr.next
        for line_str in fr:
            words = line_str.strip().split()
            psm_dict[words[0]] = ((words[1], ""))
        
    with open(rt_file, 'r') as fr:
        for line_str in fr:
            words = line_str.strip().split()
            if words[0] not in psm_dict:
                print 'error.'
            else:
                psm_dict[words[0]] = ((psm_dict[words[0]][0], words[1]))
    
    with open(elude_rt_file, 'w') as fw:
        for _key, value in psm_dict.iteritems():
            fw.write(value[0])
            fw.write('\t')
            fw.write(value[1])
            fw.write('\n')
            
def get_RTPredict_scatter_data(sId):
    input_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/RTmodel/result_"+sId+".csv"
    varify_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Data/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_"+sId+".txt"
    fwr_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/RTmodel/r/fwr_"+sId+".txt"
    shu_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/RTmodel/r/shu_"+sId+".txt"
    
    pep_dict = {}
    
    shift_float = -5
    
    with open(varify_file, 'r') as fr:
        for line_str in fr:
            words = line_str.strip().split()
            time_float = float(words[1])
            # time_float += shift_float + 120
            # time_float %= 120
            pep_dict[words[0]] = (str(time_float), words[2])
    
    fw_fwr = open(fwr_file, 'w')
    fw_shu = open(shu_file, 'w')
    fw_fwr.write("real\tpredict\n")
    fw_shu.write("real\tpredict\n")
    with open(input_file, 'r') as fr:
        for line_str in fr:
            words = line_str.strip().split()
            if words[0] not in pep_dict:
                print 'error.'
            else:
                info = pep_dict[words[0]]
                if info[1] == '1':
                    fw_fwr.write(info[0])
                    fw_fwr.write('\t')
                    fw_fwr.write(words[1])
                    fw_fwr.write('\n')
                else:
                    fw_shu.write(info[0])
                    fw_shu.write('\t')
                    fw_shu.write(words[1])
                    fw_shu.write('\n')
    
    fw_fwr.close()
    fw_shu.close()

def get_elude_scatter_data(sId):
    input_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Elude/result_"+sId+".txt"
    varify_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Data/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_"+sId+".txt"
    fwr_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Elude/r/fwr_"+sId+".txt"
    shu_file = "/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/RTime/Elude/r/shu_"+sId+".txt"
    
    pep_dict = {}
    
    shift_float = -5
    
    with open(varify_file, 'r') as fr:
        for line_str in fr:
            words = line_str.strip().split()
            time_float = float(words[1])
            # time_float += shift_float + 120
            # time_float %= 120
            pep_dict[words[0].replace('(', '[').replace(')', ']')] = (str(time_float), words[2])
    
    fw_fwr = open(fwr_file, 'w')
    fw_shu = open(shu_file, 'w')
    fw_fwr.write("real\tpredict\n")
    fw_shu.write("real\tpredict\n")
    with open(input_file, 'r') as fr:
        fr.next()
        fr.next()
        fr.next()
        for line_str in fr:
            words = line_str.strip().split()
            if words[0] not in pep_dict:
                print 'error.'
            else:
                info = pep_dict[words[0]]
                if info[1] == '1':
                    fw_fwr.write(info[0])
                    fw_fwr.write('\t')
                    fw_fwr.write(words[1])
                    fw_fwr.write('\n')
                else:
                    fw_shu.write(info[0])
                    fw_shu.write('\t')
                    fw_shu.write(words[1])
                    fw_shu.write('\n')
    
    fw_fwr.close()
    fw_shu.close()

def main(argv=None):
    for i in range(1, 11):
        sId = "%02d" % i
        # get_cross_validate_data(sId)
        # get_RTPredict_scatter_data(sId)
        get_elude_scatter_data(sId)
    # get_elude_RTpredict_data()
    print 'Done.'

if __name__ == '__main__':
    sys.exit(main())