import sys,os

for i in range(2,94):
    path = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/lagtime50/bootstrap/%dtrajs/'%i
    for j in range(1,101):
        afn = os.path.join(path,'Assignments%d_bs.h5'%j)
        ofn = os.path.join(path,'ImpliedTimescales_%d.dat'%j)
        print 'Run CalculateImpliedTimescales.py -l 50,100 -i 50 -a %s -o %s -p 8 -e 1'%(afn,ofn)
        os.system('CalculateImpliedTimescales.py -l 50,100 -i 50 -a %s -o %s -p 8 -e 1'%(afn,ofn))