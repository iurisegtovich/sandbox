import os,sys

if len(sys.argv) < 3:
    print "Usage: python submitjobs.py startjobnumber endjobnumber"
    sys.exit()
for i in range(int(sys.argv[1]),int(sys.argv[2])+1):
    cmd = "qsub qsub%d.sh"%i
    print "Submitting qsub%d.sh ..."%i 
    os.system(cmd)
print "Done!"
