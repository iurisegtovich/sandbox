from scipy import savetxt,loadtxt
for i in range(1,8):
    fn = 'surprisalvalues-weighted-ERRRMSD4.0-%d00trajs.dat'%i
    s=loadtxt(fn)
    s=s.mean(0)
    fn1='meansurprisalvalues-weighted-ERRRMSD4.0-%d00trajs.dat'%i
    savetxt(fn1,s)
