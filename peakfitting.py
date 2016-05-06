# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 15:12:12 2015

@author: perumal
"""

import numpy as np
import re
import os
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
from lmfit.models import GaussianModel, LorentzianModel, PseudoVoigtModel 

###############################################################################

###############################################################################

os.chdir('D:/Synchrotron Data/2015_10/')

def readfio(filename):
    """
    READ *.fio file into dictionary of arrays
    """
    cols = []
    # Open file
    f = open(filename +'.fio')
    counter = 0
    for n, line in enumerate(f):
        #if "Col" in line:   # This command can also be replaced with the below line 
        if line.strip().startswith('Col'):
            linestoskip = n +  1
            cols.append(line.split()[2]) #To have only the column headers in a list
            counter = counter + 1 #gives the number of columns, in principle it is not necessary
    f.close() # close the file 
    # Read numerical data without header
    print (cols)
    data = np.genfromtxt(filename+'.fio',skip_header=linestoskip, skip_footer =1)
    return data, cols

###############################################################################

###############################################################################

def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(np.sqrt(2*np.pi)*wid))*np.exp(-(x-cen)**2/(2*wid**2))



###############################################################################

###############################################################################


def plotsample(filename):
    data, cols = readfio(filename)
    y = data[:,cols.index('signalcounter_atten')]
    
    #x = data[:,1] # Chi angle    
    #x = (2*math.pi/wavelength)*(math.sin(math.radians(data[:,1]))+ math.sin(math.radians(data[:,2])-math.radians(data[:,1])))
    x = (2*np.pi/wavelength)*(np.sin(np.radians(data[:,cols.index('om')]))+
    np.sin(np.radians(data[:,cols.index('tt')])-np.radians(data[:,cols.index('om')])))        
    plt.semilogy(x,y,'go-',label=filename)
    plt.xlabel(r'$Q_{z}$ ($\AA^{-1}$)')
    plt.ylabel('Intensity [counts/sec]')
    plt.legend(prop={'size':16},frameon=False,numpoints=1)
    plt.show()
    plt.clf()

###############################################################################

###############################################################################

def fitsample(filename):
    data, cols = readfio(filename)
    y = data[:,cols.index('signalcounter_atten')]
    x = (2*np.pi/wavelength)*(np.sin(np.radians(data[:,cols.index('om')]))+
    np.sin(np.radians(data[:,cols.index('tt')])-np.radians(data[:,cols.index('om')])))        
    
    # create a mask to select only interesting data range
    m = (x > 0.25) & (x < 0.45)
    x_fit = x[m]
    y_fit = y[m]
    
    mod = GaussianModel(prefix='gau1_')
    pars = mod.guess(y_fit, x=x_fit)
    pars['gau1_center'].set(0.36, min = 100, max = 10000)
    out = mod.fit(y_fit, pars, x=x_fit)
    mod.set_param_hint('cen', value = 0.36, min=0.25, max =0.5)
    
    print(out.fit_report())
    
    plt.xlim([0.2, 0.5])
    plt.ylim([10, 50000])
    plt.semilogy(x, y, 'bo')
    
    plt.semilogy (x_fit, out.init_fit, 'k--')    
    plt.semilogy(x_fit, out.best_fit, 'r-')
    plt.show()
    

###############################################################################

###############################################################################


###############################################################################

###############################################################################


filename = 'p08_126_GSTSi111_00366'
energy = 18.000 # in KeV
def wavelength(energy):
    #Energy in KeV, so that's why eV is multiplied by 1000 
    #meter has to be converted to Angstroms, so multiplied by 1e10
    #global h, c
    h = 6.626e-34  #joules
    c = 2.998e8    #m/sec
    eV = 1.602e-19 #Joules

    wavelength = h*c/(energy*1000*eV)*1e10       
    return wavelength
    
wavelength = wavelength(energy)
#plotsample(filename)
fitsample(filename)

fig = plt.figure(1)
fig.clf()
print (wavelength)
plt.show()
    
    





#filename = 'p08_126_GSTSi111_00366'
#filename = raw_input('Enter the file name: ')
readfio(filename)

