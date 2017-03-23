#load libraries, have to make sure you're in the metabolomics_gui folder
import sys
import os
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import random
import spectra
import bucket
import processing
import interpolate
import pca
import pls_old as pls

#get data
fd = '/home/louic/Desktop/vaginal_samples/vagswab/cpmg'
spectralist = spectra.BrukerDirectory(str(fd).strip())
xdata = spectralist.getAllxData()
ydata = spectralist.getAllyData()
ppm = xdata[0]

#big data file now, ydata is 637 x 65536 = 41 *10^6 datapoint. if saved becomes a 1Gb beast. 
#first thing is to get rid of shit to make it more managable
#removing parts of spectra to reduce file size
remove_list =  [ [-100, -0.4],[4.7, 5.0], [8.8, 100]]

#find indexes to remove from ppm
#idx = np.argmin(np.abs(ppm-new_ppms[i]))
#normalise data to 1000
ydata = ydata/ (np.max(ydata/1000.0))
remove_list_indexes = []
for i in remove_list:
  j = []
  idx0 = np.argmin(np.abs(ppm-[i[0]]))
  idx1 = np.argmin(np.abs(ppm-[i[1]]))
  remove_list_indexes.append([idx0,idx1])
  
  
#extremes = [remove_list_indexes[0][1], remove_list_indexes[2][0]]  

#xdata = xdata[:,extremes[1]:extremes[0]]
#ydata = ydata[:,extremes[1]:extremes[0]]
#ppm = ppm[extremes[1]:extremes[0]]

for i in remove_list_indexes:
  ydata= np.delete(ydata, range(i[1],i[0]),1)
  xdata = np.delete(xdata, range(i[1],i[0]),1)
  
#chop off last index
n,m = ydata.shape
ydata = ydata[:,:m-2]
xdata = xdata[:,:m-2]
#now data is 637 x 30922 so almost half the size, still saves as a 500Mb file

tok = np.vstack((xdata[0], ydata))
np.savetxt(fd+'/data', tok, delimiter ='\t')
#when reading this file make sure to declare the delimiter
#otherwise gives invalid literal for float error. 
data=  np.loadtxt('data', delimiter = '\t')

#get list of spectra
i=0
titles = []
for s in spectralist.spectra():
  s.setID(str(i))
  title = s.title().split('\n')[0]
  item= "%2s %20s %20s" %(s.ID(), title, s.spectrumFile())
  
  titles.append(item)
import pandas as pd
titles = pd.DataFrame(titles)
title.to_csv(fd+'/titles')
