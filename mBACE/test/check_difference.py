import os,sys

for i in xrange(2,900):
    os.system("diff Output_BACE_s/map{mapid}.dat Output_BACE/map{mapid}.dat".format(mapid=i))
