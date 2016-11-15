'''
Created on Aug 10, 2016

@author: xgo
'''

if __name__ == '__main__':
    
    sipros10_scan_list = []
    
    with open('/home/xgo/Temp/Sipros10Debug/Temp/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_01.AngeloSipros10_Spe2Pep.txt.tab', 'r') as f:
        f.next()
        for line in f:
            l = line.split('\t')
            sipros10_scan_list.append(int(l[1]))
    
    missed_scan_list = []
    
    all_scan_list = []
    
    with open('/home/xgo/Temp/Sipros10Debug/Ft2/Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_01.FT2','r') as f:
        for line in f:
            if line.startswith('S'):
                l = line.split('\t')
                all_scan_list.append(int(l[1]))
                if not int(l[1]) in sipros10_scan_list:
                    missed_scan_list.append(int(l[1]))

    fout = open('/home/xgo/Temp/Sipros10Debug/Temp/missed_scans_by_sipros10_in_iProphet.txt', 'w')
    fout2 = open('/home/xgo/Temp/Sipros10Debug/Temp/missed_scans_by_iProphet.txt', 'w')
    f_iProphet = open('/media/xgo/Seagate/Proteomics/Experiments/Angelo/iProphet/Angelo_iProphet_0.8.xls', 'r')
    f_iProphet.next()
    for line in f_iProphet:
        l = line.split('\t')
        name = l[1].split('.')
        if name[0] == 'Angelo_10022013_P2_1020cm_MB_FASP_Elite_Run2_01':
            scan_id = int(l[2])
            if scan_id in missed_scan_list:
                fout.write(l[2])
                fout.write('\n')
            if not scan_id in all_scan_list:
                fout2.write(l[2])
                fout2.write('\n')
    f_iProphet.close()
    fout.close()
    fout2.close()
    print 'Done'