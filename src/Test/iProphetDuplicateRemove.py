'''
Created on Aug 11, 2016

@author: xgo
'''

class PSM:
    
    def __init__(self):
        self.line = ''
        self.scan_number = 0
        self.iProb = 0.0
        
    
    def set(self, line = '', words=[]):
        self.line = line
        self.scan_number = int(words[2])
        self.iProb = float(words[3])
        
    def update(self, line='', words=[]):
        iProb = float(words[3])
        if self.scan_number != int(words[2]):
            print "Error"
        if iProb > self.iProb:
            self.line = line
            self.iProb = iProb

if __name__ == '__main__':
    FT2_list = []
    PSM_list = []
    for _i in range(12):
        FT2_list.append([])
        PSM_list.append([])
    
    
    scan_number = 0
    FT2_id = 0
    count = 0
    duplicate_count = 0
    
    with open('/home/xgo/Temp/iProphetCheck/myrimatch.interact.ipro.full.xls', 'r') as fin:
        fin.next()
        for line in fin:
            words = line.split('\t')
            scan_number = int(words[2])
            FT2_id = int((words[1].split('.')[0])[-2:])
            if scan_number in FT2_list[FT2_id]:
                index = FT2_list[FT2_id].index(scan_number)
                PSM_list[FT2_id][index].update(line, words)
                duplicate_count += 1
            else:
                FT2_list[FT2_id].append(scan_number)
                psm = PSM()
                psm.set(line, words)
                PSM_list[FT2_id].append(psm)
            count += 1
            
            if count % 1000 == 0:
                print count
    
    fout = open('/home/xgo/Temp/iProphetCheck/myrimatch.interact.ipro.full_duplicate_removed.txt', 'w')
    for l in PSM_list:
        for psm in l:
            fout.write(psm.line)        
    
    fout.close()
    
    print 'Duplicate: %s' % str(duplicate_count)
    
    print 'Done.'