#===============================================================================
# ssaCalculator.py
#
# Please reference
# N Singhal Hinrichs, and VS Pande. J. Chem. Phys. 2007. Calculation of the 
# distribution of eigenvalues and eigenvectors in Markovian State Models for
# molecular dynamics.
#
# Written 10/1/10 by
# Dan Ensign <densign@mail.utexas.edu>
# Gregory Bowman <gregoryrbowman@gmail.com>
# Sergio Bacallado <sergiobacallado@gmail.com>
# Stanford University
# Pande group
#
# Copyright (C) 2008  Stanford University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#===============================================================================
# TODO:
#===============================================================================
# GLOBAL IMPORTS:
from math import log
from numpy import real, round, float64, common_type, ones
from os.path import exists
from os import mkdir
from pickle import Pickler
from scipy.linalg import eigvals, lu, solve, LinAlgError
from scipy import argmax, diag, identity, int32, matrix, zeros
import scipy.sparse.linalg.eigen
from scipy.io import mmread, mmwrite

import numpy as np  # VAV
import sys
#===============================================================================
# LOCAL IMPORTS:
from ssaTools import *
#===============================================================================


class ssaCalculator( object ):
        """ Class for objects to govern sensitivity calculations. """
        def __init__( self, lagTime, tcMatrix, priorCounts=.5, evalList=[1], nNewSamples=1, timeUnits="ns", recommendationScheme = 'Nina' ):
                self.lagTime = lagTime          # a number in appropriate units 

                self.tcMatrix = tcMatrix        # a regular list of lists of numbers, square

                self.priorCounts = priorCounts  # this integer is the uniform alpha parameter of the Dirichlet prior
						# which must be greater than 0. By default, it is 0.5

                self.evalList = evalList	# the list of eigenvalues (as arguments of a sorted list from largest to
						# smallest) to be used in sensitivity computations, by default
						# only the second largest eigenvalue is used, [1].

                self.nNewSamples = nNewSamples  # the number of new samples that we plan to simulate.

                self.recommendationScheme = recommendationScheme

		self.timeUnits = timeUnits	

                self.evalList.sort()

                self.qlists = []
                self.wlists = []
                self.varianceContributions = []         # a list of sundry results

		# derived properties
                try:
                    self.dim = len( self.tcMatrix )
                except:
                    self.dim = self.tcMatrix.shape[0]
		self.weights = self.computeWeights()	# the total number of transitions observed yet from each state
              
                self.compute()


	def computeWeights( self ):

                if scipy.sparse.isspmatrix(self.tcMatrix):
                    C=scipy.sparse.csr_matrix(self.tcMatrix).asfptype()
                    weights = np.asarray(C.sum(axis=1)).flatten()

                else:
                    # sum the rows of the tcMatrix and return as a matrix
		    weights = [] 
		    for i in range( self.dim ):
			try:
                            row = self.tcMatrix[ i ]
                        except:
                            row = self.tcMatrix[ i,: ].sum()

			rowsum = float( matrix( row ).sum() )
			weights.append( rowsum )

		return matrix( weights )

        def compute( self, Verbose=False):
                # 1. generate the Dirichlet matrix 
                self.dirmat = self.computeDirichletMatrix()
		del self.tcMatrix

                # 2. compute the eigenvalues
                self.avgEvals = list( eigvals( self.dirmat ) )

		# Sergio: take only the real part of the eigenvalues, such that the program
		# doesnt give weird results in computing the LU decomposition of A
                self.avgEvals = list(real(self.avgEvals))

                # VAV:  TO DO -- needs a check that the count matrix should be symmetric!
                # If not, we'll get a pair of a+bi, a-bi eigenvalues, so that the real components are the same.

		#self.avgEvals = list( arpack.eigen( self.dirmat, 1+max(self.evalList) )[0] )
		self.avgEvals.sort()
		self.avgEvals.reverse() # longest time scales first

		f = open("eigVals.dat", 'a')
		for ind in range(len(self.avgEvals)):
		  f.write(str(float(self.avgEvals[ind]))+" ")
		f.write("\n")
		f.close()

                for evalIndex in self.evalList:
                        # print '\n*** computing uncertainties for evalIndex =', evalIndex, 
                        evalue = self.avgEvals[ evalIndex ]
			
                        # 3. compute the characteristic polynomial-generating matrix for the average matrix
#                        A = self.computeCharPolGenerator( evalue )
			# to save memory will modify dirmat and then fix later
			for i in range(self.dim):
				self.dirmat[i,i] -= evalue

                        # 4. Decompose the characteristic polynomial generator
                        ( P, L, U, Q ) = self.computeDecomposition( self.dirmat )

			# fix dirmat
			for i in range(self.dim):
				self.dirmat[i,i] += evalue
			
			# 5. solve for the projection vectors x1 and x2
			( x1, x2 ) = self.computeProjectionVectors( P, L, U )
			del P
			del L
			del U

			# 6. compute the matrix of first derivatives
			normalization = x2.T * x1
			myNumerator = Q * x2 * x1.T * Q 
			del x1
			del x2
			firstDerivatives = myNumerator / normalization
			if Verbose:
                            print "myNumerator", myNumerator
			    print "normalization", normalization
			if normalization==0: 
				print "ALERT"
			del normalization
			del myNumerator
			del Q

			# 7. compute the contributions of each state to the eigenvalue uncertainty
			# and the reductions given 1 additional sample and m additional samples
			self.qlists.append( self.computeVarianceContributions( firstDerivatives ) )
			self.wlists.append( self.weights )

                        qlist, wlist = self.qlists[-1], self.wlists[-1]
			self.computeReductions( qlist, wlist )

                        # Added by Sergio June/10
                        # Store raw qlist in object
                        #self.qlist = qlist   # removed by VAV



        def computeDirichletMatrix( self ):

            if scipy.sparse.isspmatrix(self.tcMatrix):
                print "Currently sparse matrices are not supported! Converting to dense..."
                self.tcMatrix = self.tcMatrix.todense() 
 
            if (0):
                dirmat = zeros([self.dim,self.dim], float64)
                i = 0
                while i < self.dim :

                        # big note = if we didn't observe the state, then we assume that the transition
                        # probability of the state to itself is 1.0 -- therefore it won't show up as
                        # influencing the uncertainty. We also do not include the prior counts for such
                        # a state.
 
                        """
                        if self.weights[0,i] == 0 :
                                dirmat[ i, i ] = 1
                                
                        else:
                                row = arr2lst( self.tcMatrix[i] + self.priorCounts )
                                frow = DirichletDistribution( row )

                                for g in range( len( row ) ) :
                                        dirmat[i,g] = frow.mean( g )
                        """
                        # Sergio: I am removing Dan's modification to the algorithm for states that 
                        # haven't been observed.

                        # VAV: this is needed for sparse matrices (i.e. only add prior counts for observed transitons) 
                        thisrow = self.tcMatrix[i,:].todense()
                        row = arr2lst( thisrow + self.priorCounts )

                        # Sergio: I'm trying a new thing here: adding 2 prior counts to self-transition
                        # probabilities to try to avoid the possible effect of "traps"
                        # row[i] += 2

                        frow = DirichletDistribution( row )

                        for g in range( len( row ) ) :
                            dirmat[i,g] = frow.mean( g )

                        i += 1


            if (1):  
		dirmat = zeros([self.dim,self.dim], float64)
		i = 0
		while i < self.dim :
			
			# big note = if we didn't observe the state, then we assume that the transition
        		# probability of the state to itself is 1.0 -- therefore it won't show up as
        		# influencing the uncertainty. We also do not include the prior counts for such
			# a state.
			"""
			if self.weights[0,i] == 0 :
				dirmat[ i, i ] = 1
				
			else:
				row = arr2lst( self.tcMatrix[i] + self.priorCounts )
				frow = DirichletDistribution( row )

				for g in range( len( row ) ) :
					dirmat[i,g] = frow.mean( g )
			"""
			# Sergio: I am removing Dan's modification to the algorithm for states that 
			# haven't been observed.

                        # VAV: this works for dense matrices
			row = arr2lst( self.tcMatrix[i] + self.priorCounts )

			
                        # Sergio: I'm trying a new thing here: adding 2 prior counts to self-transition
			# probabilities to try to avoid the possible effect of "traps"
			# row[i] += 2

			frow = DirichletDistribution( row )

			for g in range( len( row ) ) :
			    dirmat[i,g] = frow.mean( g )

			i += 1

		return dirmat

        def computeCharPolGenerator( self, evalue ):
		A = self.dirmat.copy()
		for i in range(self.dim):
			A[i,i] -= evalue
                return A

        def computeDecomposition( self, A ):
		# first decomposition
		( P, L, U ) = lu( A )
	
		smallestIndex = findsmallestdiag( U )
		smallest = U[ smallestIndex, smallestIndex ]

		Q = identity( self.dim, dtype=float64 )
		if smallestIndex+1 != self.dim :
			del P
			del L
			del U

			# exchange smallestIndex row with dim-1 row
			# multiplying A by this matrix on both sides (Q.A.Q) will ensure that the
			# smallest element of U will be in the lower right corner.
			swaprow( Q, smallestIndex, self.dim-1 )
			
			# recompute the decomposition
			Q = matrix( Q )
			A = matrix( A )
			( P, L, U ) = lu( Q*A*Q )
	
		return ( P, L, U, Q )

	def computeProjectionVectors( self, P, L, U ) :	
		eK = matrix( identity( self.dim, float64 )[ 0: ,( self.dim - 1 ) ] ).T
		U = matrix(U, float64)
		U[ self.dim - 1, self.dim - 1 ] = 1.0
		# Sergio: I added this exception because in rare cases, the matrix
		# U is singular, which gives rise to a LinAlgError.
		try: 
			x1 = matrix( solve( U, eK ), float64 )
		except LinAlgError:
			print "Matrix U was singular, so we input a fake x1\n"
			print "U: ", U
			x1 = matrix(ones(self.dim))

		#print "x1", x1
		del U

		LT = matrix( L, float64, copy=False ).T
		PT = matrix( P, float64, copy=False ).T

		x2 = matrix( solve( LT*PT, eK ), float64 )
		del L
		del P
		del LT
		del PT
		del eK

		return ( x1, x2 )

	def computeVarianceContributions( self, firstDerivatives ) :

		qlist = []
		i = 0
		while i < self.dim :
			#print "Compute var %d" % i
			# derivatives of the eigenvalue for this state
			s = matrix( firstDerivatives[ i,0: ], float64 ).T
			##print s.shape

			# cross probability matrix for this state
			Pi = matrix( self.dirmat[ i,0: ], float64 ).T
			##print Pi.shape
			part1 = diag( arr2lst( Pi.T ) )
			part1 = matrix(part1, float64)
			##print part1.shape
			##print common_type(part1)
			Cp = matrix( part1 - Pi * Pi.T )
			##print common_type(Cp)
			##print Cp.shape
			del part1

			# degree of sensitivity for this state
			q = float( abs( s.T * Cp * s ) )
			del s
			del Pi
			del Cp

			qlist.append( q )

			i += 1

		return matrix( qlist )	

	def computeReductions( self, qlist, wlist ) :

		deltalist = qlist / ( 1. + self.weights )

		deltalist_1add = qlist / ( 2. + self.weights )
		reductionlist = deltalist - deltalist_1add
		fractionlist = reductionlist/reductionlist.sum()
		recclist = matrix( identity( self.dim )[ argmax( fractionlist ),0: ], "int" )

		deltalist_1add = arr2lst( deltalist_1add )
		reductionlist = arr2lst( reductionlist )
		fractionlist = arr2lst( fractionlist )
		recclist = arr2lst( recclist )

                # These are the recommendations for only one new sample
		currentEvalVariance = [ deltalist_1add, reductionlist, fractionlist, recclist ]

                # If you want more samples, add four more entries to this list, with this info
 		if self.nNewSamples > 1 :
			deltalist_madd = qlist / ( 1. + self.nNewSamples + self.weights ) 
			reductionlist_m = deltalist - deltalist_madd
			self.stateToSampleMore = argmax(reductionlist_m)
                        fractionlist_m = reductionlist_m / reductionlist_m.sum()        

                        if self.recommendationScheme == 'Nina':
                            # Nina's scheme: put all nNewSamples at one state
                            recclist_m = matrix( self.nNewSamples*identity( self.dim )[ argmax( fractionlist_m ),0: ], "int" )

                        elif self.recommendationScheme == 'VAV':
                            # VAV: a new possible scheme -- pick several to sample, but make sure the total counts are self.nNewSamples
			    recclist_m = matrix( self.nNewSamples* fractionlist_m , "int" )
                            while recclist_m[0,:].sum() < self.nNewSamples:
                                recclist_m[ 0, argmax( fractionlist_m ) ] += 1       

			
			deltalist_madd = arr2lst( deltalist_madd )
			reductionlist_m = arr2lst( reductionlist_m )
			fractionlist_m  = arr2lst( fractionlist_m )
			recclist_m = arr2lst( recclist_m )

			currentEvalVariance.extend( [ deltalist_madd, reductionlist_m, fractionlist_m, recclist_m] ) 

                # add two more things to the beginning of the list: deltalist,  and 0....n-1	
		currentEvalVariance.insert( 0, arr2lst( deltalist ) )
		currentEvalVariance.insert( 0, range( self.dim ) )   # 0 .... n-1 will be the firstthing

                #Now, transpose the whole thing so that this information if in column format 
		currentEvalVariance = matrix( currentEvalVariance, "float64" ).T   
	
		self.varianceContributions.append( currentEvalVariance )

	def displayContributions( self, bling = False ):
		if bling:
			sep1 = ":"
			sep2 = "|" 
		else:
			sep1 = " "
			sep2 = " "

		i = 0
		for evalIndex in self.evalList :

			evalue = self.avgEvals[ evalIndex ] 

			if bling:
				s1 = "eigenvalue %6.6f" % evalue

				try:
					timescale = -self.lagTime/log( evalue )
					s1 += "   time scale %7.4f %s" % ( timescale, self.timeUnits )

				except ZeroDivisionError :
					s1 += "   time scale Inf %s" % self.timeUnits 

				except ValueError :
					pass
				
				#print s1
				#print "-"*len( s1 ) 

			# here are the results for this eigenvalue
			s2 = ""
			contributions = self.varianceContributions[ i ] 

			frm = "%6.6e " * 3
			v = "%d" % ( int( log( self.nNewSamples ) + 1 ) ) 
			frm += "%" + v + "d" 
	
			for j in range( contributions.shape[0] ):
				row = tuple( arr2lst( contributions[ j,0: ] ) )

				rowfrm = "%3d " + sep1 + " %6.6e " + sep2 + " \t" + frm
				if self.nNewSamples  > 1 :
					rowfrm += "\t "+ sep2 + " " + frm
				#print rowfrm % row

			if bling :	
				print
			i += 1

	def resamplingDistribution( self, eigIndex ):
		# just return the list of suggested new samples, to reduce the variance of that eigenvalue
		evalVarianceData = self.varianceContributions[ eigIndex ]
     
                return [int(x) for x in list(evalVarianceData[ :, -1 ])]   # recommendations for m new samples 
                #return [self.nNewSamples*int(x) for x in list(evalVarianceData[ :, 5 ])]   # recommendations for 1 new sample

        def sensitivities( self, eigIndex ):
                evalVarianceData = self.varianceContributions[ eigIndex ].T

                # take the i=2 row, which contain the deltalist components of the variance 
                #print 'evalVarianceData', evalVarianceData
                return evalVarianceData[ -4,: ]

        def uncertainty_variance( self, eigIndex ):
                evalVarianceData = self.varianceContributions[ eigIndex ].T

                # take the i=2 row, which contain the deltalist components of the variance 
                #print 'evalVarianceData', evalVarianceData
                return np.sum(evalVarianceData[ -4,: ]) 


	def save( self, filename ):
		FILE = open( "%s/%s" % ( self.outputDir, filename ), "w" )
		p = Pickler( FILE )
		p.dump( self )
		FILE.close()	
		return True	
