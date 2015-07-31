import os, sys, glob, string
import cPickle
import numpy as np


class Experiment(object):
    """An object to store and direct adaptive sampling experiments."""

    def __init__(self, experimentFn):
        """Create an Experiment object from file, or if it doesn't exist, make a new one."""

        # The Experiment filename
        self.experimentFn = experimentFn

        # Default settings
        self.parms = {}
        self.parms['adaptive_method'] = 'JSD' #'JSD', 'JSD variance', 'random', 'continue traj', 'surprisal', 'surprisal variance', 'pi_si'
        self.parms['sampling_method'] = 'uniform' # 'uniform', 'traj'
        self.parms['nrounds'] = 100
        self.parms['nsamples'] = 100  # samples per round
        self.parms['npresamples'] = 10 # number of pre-samples for each state (before adaptive sampling begins)
        self.parms['system'] = 'hairpin-macro150' # 'toy' or 'dipep' or 'hairpin-macro150'
      
        if not os.path.exists(self.experimentFn):
            self.write_expfile(self.experimentFn)
        else:
            self.read_expfile(self.experimentFn)

        return

    def write_expfile(self, fn):
        """Write experimental setting to file."""
        fout = open(fn, 'w')
        for key, value in self.parms.iteritems(): 
            fout.write('%-24s'%key+' '+repr(value)+'\n')
        fout.close()
        
    def read_expfile(self, fn):
        """Read exp settings from file."""
        fin = open(fn, 'r')
        lines = fin.readlines() 
        for line in lines:
            fields = line.strip().split()
            self.parms[fields[0]] = eval(string.joinfields(fields[1:],' '))
        fin.close() 

