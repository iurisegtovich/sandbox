import sys,os
sys.path.append('User/tud51931/voelz/sandbox/gfzhou/test')
import bootstrap
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input',help='Input file for bootstraping',required = True)
parser.add_argument('-n','--NumOfTrajs',help="The Number of Trajectories used for bootstrapping in input file. Also, you can use 'All' to use all trajectories to bootstrap. Default: 'ALL'",default='ALL')
parser.add_argument('-b','--BootstrapNum',help="The number of bootstraped files.Default: 100",default=100,type=int)
parser.add_argument('-o','--Output_Dir',help="Output directory.Default: .InputfileDir/bs-Inputfilename/")

args = parser.parse_args()
if args.Output_Dir == None:
    args.Output_Dir = os.path.join(os.path.dirname(os.path.abspath(args.Input)),'bs-%s'%os.path.basename(args.Input).split('.')[0])
bootstrap.bootstrap(args.Input,args.NumOfTrajs, args.BootstrapNum, args.Output_Dir)
       
