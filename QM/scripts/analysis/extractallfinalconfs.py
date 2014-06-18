import os,sys
import commands

logs = commands.getoutput('ls logs/*.log').split('\n')

for log in logs:
    cmd='python extractfinalconf.py %s'%log
    os.system(cmd)
