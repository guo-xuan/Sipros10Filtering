'''
Created on Nov 2, 2016

@author: xgo
'''

import sys
from sets import Set

def main(argv=None):
    LR = open("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/Ensamble/LR_dec.txt", 'r')
    RF = open("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/Ensamble/RF_dec.txt", 'r')
    AB = open("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/Ensamble/AB_dec.txt", 'r')
    DP = open("/media/xgo/Seagate/Proteomics/Experiments/Angelo/Sipros10/Ensamble/DP_dec.txt", 'r')
    
    LR_set = Set()
    RF_set = Set()
    AB_set = Set()
    DP_set = Set()
    All_set = Set()
    
    for line_str in LR:
        LR_set.add(line_str.strip())
        All_set.add(line_str.strip())
        
    for line_str in RF:
        RF_set.add(line_str.strip())
        All_set.add(line_str.strip())
        
    for line_str in AB:
        AB_set.add(line_str.strip())
        All_set.add(line_str.strip())
        
    for line_str in DP:
        DP_set.add(line_str.strip())
        All_set.add(line_str.strip())
    
    print 'All:\t%d' % len(All_set)
    print 'LR:\t%d' % len(LR_set)
    print 'RF:\t%d' % len(RF_set)
    print 'AB:\t%d' % len(AB_set)
    print 'DP:\t%d' % len(DP_set)
    
    LR_RF = 0
    LR_AB = 0
    LR_DP = 0
    RF_AB = 0
    RF_DP = 0
    AB_DP = 0
    
    LR_RF_AB = 0
    LR_RF_DP = 0
    LR_AB_DP = 0
    RF_AB_DP = 0
    
    LR_RF_AB_DP = 0
    
    for oPsm in LR_set:
        if oPsm in RF_set:
            LR_RF += 1
            if oPsm in AB_set:
                LR_RF_AB += 1
                if oPsm in DP_set:
                    LR_RF_AB_DP += 1
            if oPsm in DP_set:
                LR_RF_DP += 1
        if oPsm in AB_set:
            LR_AB += 1
            if oPsm in DP_set:
                LR_AB_DP += 1
        if oPsm in DP_set:
            LR_DP += 1
    
    for oPsm in RF_set:
        if oPsm in AB_set:
            RF_AB += 1
            if oPsm in DP_set:
                RF_AB_DP += 1
        if oPsm in DP_set:
            RF_DP += 1
    
    for oPsm in AB_set:
        if oPsm in DP_set:
            AB_DP += 1
    
    print "LR&RF:\t%d" % LR_RF
    print "LR&AB:\t%d" % LR_AB
    print "LR&DP:\t%d" % LR_DP
    print "RF&AB:\t%d" % RF_AB
    print "RF&DP:\t%d" % RF_DP
    print "AB&DP:\t%d" % AB_DP
    
    print "LR&RF&AB:\t%d" % LR_RF_AB
    print "LR&RF&DP:\t%d" % LR_RF_DP
    print "LR&AB&DP:\t%d" % LR_AB_DP
    print "RF&AB&DP:\t%d" % RF_AB_DP
    
    print "LR&RF&AB&DP:\t%d" % LR_RF_AB_DP
    
    pass

if __name__ == '__main__':
    sys.exit(main())