import os,sys
from sklearn.externals import joblib
import matplotlib.pyplot as plt
import numpy as np

map_id = 40
pids = range(6383,6391)
mut = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
p = 1./40.*np.ones(40)

plt.figure(figsize=(28,14))
for i,p_id in enumerate(pids):
    
    model_dir = "MSMs-%d-macro%d"%(p_id,map_id)
    if not os.path.exists(model_dir):
        print "%s doesn't exist! Exit!"%model_dir
    msm_fn = os.path.join(model_dir,"MSMs-%d-macro%d.pkl"%(p_id,map_id))
    msm = joblib.load(msm_fn)
    for j in range(5):
        ax = plt.subplot(5,8,i+1+j*8)
        #plt.plot(range(len(msm.left_eigenvectors_[:,j])),msm.left_eigenvectors_[:,j],"ro--")
        #plt.plot(range(len(msm.left_eigenvectors_[:,j])),[0]*len(msm.left_eigenvectors_[:,j]),"k--")
        plt.hold(True)
        if np.dot(p,msm.right_eigenvectors_[:,j]) > 0 :
            for stateid in range(40):
                plt.vlines(stateid,min(0,msm.left_eigenvectors_[stateid,j]),max(msm.left_eigenvectors_[stateid,j],0),linestyles='solid',linewidth="2")
        else:
            for stateid in range(40):
                plt.vlines(stateid,min(0,-msm.left_eigenvectors_[stateid,j]),max(-msm.left_eigenvectors_[stateid,j],0),linestyles='solid',linewidth="2")
        #plt.stem(msm.left_eigenvectors_[:,j],linefmt='k-', markerfmt="|", basefmt='r-')
        plt.subplots_adjust(hspace = .1)
        locs,labels = plt.yticks()
        plt.yticks(locs,map(lambda x: "%.2f" % x, locs))
        #plt.ylim(min(msm.left_eigenvectors_[:,j])-0.02,max(msm.left_eigenvectors_[stateid,j])+0.02)
        #ax.margins(x=None,y=0.02)
        if j != 4:
           ax.set_xticklabels(()) 
        if j == 0:
            plt.title("Fs-%s"%mut[i])
        if i == 0:
            plt.ylabel("$\phi_%d$"%j)
         
        #plt.show()
plt.savefig("eighenvectors_all_together_sign_fixed.pdf")
   



