
__metaclass__= type

class Cutoff:
    def __init__(self,startcutoff,endcutoff,step):
        """
           Initialize the cutoffs list.
 
           startcutoff - the first cutoff in the list.
           endcutoff - the last cutoff in the list
           step - the difference between two adjacent values
           
           Note : 
               1.the endcutoff will not be included in the list if the startcutoff + n*step is bigger than the end cutoff.
               2.if there is no input or the startcutoff > endcutoff, empty list will be created.
        """        
        self.cutoffs = []
        i = startcutoff
        while i <= endcutoff:
            self.cutoffs.append(round(i,1))
            i += step
       
    
    def ext(self,startcutoff,endcutoff,step):
        """
        Extend the cutoffs list.
        
        startcutoff - the first cutoff in the list.
        endcutoff - the last cutoff in the list
        step - the difference between two adjacent values        
        
        """
        newcutoff = Cutoff(startcutoff,endcutoff,step)
        oldcutoffs = self.cutoffs[:]
        for x in newcutoff.cutoffs:
            if not self.cutoffs.count(x):
                self.cutoffs.append(round(x,1))
        self.cutoffs.sort()
        return self.cutoffs
    
    def dellist(self,startcutoff,endcutoff,step,delList,mode=0):
        """
        Delete the elements in the cutoffs list.
        mode = 0, use the start,end,step.(default)
        mode = 1, use the delList(delete List) user supplied.
        
        """
        if mode == 1:
            try:
                for x in delList:
                    self.cutoffs.remove(x) 
            except ValueError:
                print "%0.1f is not in the list, continue..."%x                
                pass
        elif mode == 0:
            delcutoff = Cutoff(startcutoff,endcutoff,step)
            try:            
                for x in delcutoff.cutoffs:
                    self.cutoffs.remove(x)
            except ValueError:
                print "%0.1f is not in the list, continue..."%x
                pass
        else:
            print"dellist usage", Cutoff.dellist.__doc__
            
def test():
    cutoff = Cutoff(2,3,0.5)
    cutoff.ext(3,4,0.5)
    print "cutoffs",cutoff.cutoffs
    cutoff.dellist(2,3,0.5)
    print 'cutoffs',cutoff.cutoffs    

if __name__ == "__main__":
    test()
    
