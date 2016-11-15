'''
Created on Aug 11, 2016

@author: xgo
'''

if __name__ == '__main__':
    FT2_list = []
    for _i in range(12):
        FT2_list.append([])
    
    scan_number = 0
    FT2_id = 0
    count = 0
    fout = open('/home/xgo/Temp/iProphetCheck/duplicate.txt', 'w')
    with open('/media/xgo/SANDISK/iProphet/Angelo/myrimatch.interact.ipro.full.xls', 'r') as fin:
        fin.next()
        for line in fin:
            words = line.split('\t')
            scan_number = int(words[2])
            FT2_id = int((words[1].split('.')[0])[-2:])
            if scan_number in FT2_list[FT2_id]:
                fout.write(words[1])
                fout.write('\t')
                fout.write(str(scan_number))
                fout.write('\n')
            else:
                FT2_list[FT2_id].append(scan_number)
            count += 1
            
            if count % 1000 == 0:
                print count
            
    fout.close()
    
    print 'Done.'