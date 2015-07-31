#!/usr/bin/env python

import numpy as np
import os, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def set_rc( plt, figsize=(1, 1), linewidth=1.5, markersize=5,
            fontsize=8.0, legend_fontsize='medium', xlabelsize = 'small', ylabelsize='small',
            fontfamily={'family':'sans-serif','sans-serif':['Arial']} ):

    plt.rc('figure', figsize=figsize)
    plt.rc('lines', linewidth=linewidth, markersize=markersize)
    plt.rc('font', size=fontsize)
    plt.rc('font',**fontfamily)
    plt.rc('xtick', labelsize=xlabelsize)
    plt.rc('ytick', labelsize=ylabelsize)
    plt.rc('legend', fontsize=legend_fontsize)


class Figure(object):

    def __init__(self, plt, figsize=(1, 1), linewidth=1.5, markersize=5,
                       fontsize=8.0, legend_fontsize='medium', xlabelsize = 'small', ylabelsize='small',
                       fontfamily={'family':'sans-serif','sans-serif':['Arial']} ):

        """Notes on Parameters
        figsize = (width, height) in inches.

        For more info on plt.rc, see:  htpp://[ADD THIS]."""

        set_rc( plt, figsize=figsize,
                linewidth=linewidth,
		markersize=markersize,
                fontsize=fontsize, legend_fontsize=legend_fontsize, xlabelsize=xlabelsize, ylabelsize=ylabelsize )
 
        # Attributes
        self.fig = plt.figure(1, figsize=figsize)  # a pyplot handle to the figure
        plt.clf()

        self.panelpos = []  # panel position -- list of elements [xorig, yorig, width, height]
        self.panelmargins = []   # list of elements [left, right, bottom, top] in fractional coords, e.g.  [0.00, 0.02, 0.18, 0.1]
        self.axis_handles = []
        self.npanels = 0  # the number of panels in the figure


    def add_panel(self, xorigin, yorigin, width, height,
                        leftmargin=0.01, rightmargin=0.01, bottommargin=0.01, topmargin=0.01):
        """Adds a new panel to the figure.

        ARGUMENTS
        xorigin, yorigin, width, height = bottom left corner x, y; width, height
        	The units of these are in fractional coordinates of the figure, i.e. ranging from 0 to 1.

        RETURNS
        ax - the axis handle of the panel."""

        self.panelpos.append([xorigin, yorigin, width, height])
        self.panelmargins.append([leftmargin, rightmargin, bottommargin, topmargin])

        # make the subplot panel
        rect = xorigin+leftmargin, yorigin+bottommargin, width-leftmargin-rightmargin, height-bottommargin-topmargin
        ax = self.fig.add_axes(rect)
        self.axis_handles.append( ax )
 
        self.npanels += 1
        return ax



if __name__ == '__main__':

    # Test the functionality of the Figure object

    f = Figure()
    f.add_panel(0.05, 0.1, 0.95, 0.24)
    f.add_panel(0.05, 0.42, 0.95, 0.24)
    f.add_panel(0.05, 0.74, 0.95, 0.24)

    mutants = ['W55F', 'F26A', 'Y31N']
    colors = ['r', 'b', 'g']
    styles = [ s+'.' for s in colors ]

    for i in range(f.npanels):
         
        x = np.arange(-10,10,0.1)
        y = np.sin(x)

        ax = plt.axes(f.axis_handles[i])
        plt.plot(x, y, color=colors[i], marker='o', markerfacecolor=colors[i], markeredgecolor='k', markersize=3, linewidth=0, markeredgewidth=0.5)
        plt.hold(True)

        #ax.set_xscale('log')
        plt.axis([-15, 15, -1.5, 1.5])  # sets the x and y range of the plot
        plt.xlabel('time (${\mu}s$)')
        plt.ylabel('relative FRET')

        ax.set_xticks([-5, -2.5, 0, 2.5, 5]) 
        ax.set_yticks([])

        plt.xlabel('[GuHCl] (M)')
        plt.ylabel('relative FRET')


plt.savefig('test.eps', format='eps')
print 'Wrote: test.eps'

