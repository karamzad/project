#!/usr/bin env python
# This is the main program 
#
import numpy as np
from pyrocko import gf, io, util, trace, pile, model, orthodrome
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Circle, Arrow
from utilities import *

####################### BODY #####################
'''

Site geometery ( corner location): lat_min, lat_max, lon_min, lon_max 
nsta : number of stations
Main shock parameters: lat_source, lon_source, depth_source, Mag, stri70ke, dip, rake
Aftershocks parameters: delta_lat, delta_lon, delta_dep (The other parameters are selected as i.e: rake(after)= rake(main)*(1+- 0.3) 
Aftershock magnitude (should be selected based on GR power low but not considered yet!) : [3:Mag]
nsource : number of aftershocks

 
'''
lat_min, lat_max, lon_min, lon_max, nsta = 10.0, 30.0, 20.0, 40.0, 8
lat_source, lon_source, depth_source, Mag, strike, dip, rake = 10., 80., 20000.0, 6.0, 22.0, 50.0, 10.0
delta_lat, delta_lon, delta_dep, nsource = 5.0, 5.0, 5000.0, 10

# make traces  
#RSource =  RandomSource(lat_source, lon_source, depth_source, Mag, strike, dip, rake, 
 #   delta_lat, delta_lon, delta_dep, nsource)  
#sources, S_lat, S_lon, S_dpths, S_strks, S_dps, S_rks, S_mgs    = RSource.RandomS()
Slat = np.array([10.0, 12.0, 11.0, 10.5, 11.2, 12.6, 9.8, 10.2, 12.1, 10.5])
Slon = np.array([80.0, 81.5, 79.5, 80.6, 81.1, 79.4, 80.7, 79.1, 79.5, 81.4])
Sdepth = np.array([20000., 15000, 12000, 11000, 25000, 13000, 20000, 18000, 22000, 15000 ])
Sstrike = np.array([ 22., 120., 80., 25., 30., 170., 23., 56., 27., 18.])
Sdip = np.array([12.0, 20., 30., 40., 10., 50., 8., 30., 40.,50.])
Srake = np.array([56.7, 10.0, 50.0, 15.0, 20.0,10.0, 40.0,30.5,15., 48.5])
Smag = np.array([7.1, 6.7, 6.1, 5.7, 5.5, 5.3, 5.1, 4.9, 4.8, 4.7])

SE = np.zeros([nsource])
sources = [
    gf.DCSource(
        time = util.str_to_time('2008-02-17 11:06:01.10'),
        lat = Slat[i],
        lon = Slon[i],
        depth = Sdepth[i],
        strike = Sstrike[i],
        dip = Sdip[i],
        rake = Srake[i],
        magnitude = Smag[i] ) for i in range(nsource)]
        
#print 'Source parameters:', S_lat, S_lon, S_dpths, S_strks, S_dps, S_rks, S_mgs

sy = np.linspace(-0.0001, 0.0001, 50)
sx = np.linspace(-0.0001, 0.0001, 50)
a = np.array([1, 30000])

deltat = 10.0/2.


#lat = np.around(np.random.uniform(lat_min, lat_max, nsta), decimals=2)
#lon = np.around(np.random.uniform(lon_min, lon_max, nsta), decimals=2)

lines = [line.strip() for line in open('parameters')] 
data = np.genfromtxt("parameters",delimiter="\n")

ind_lat = np.arange(0, 2*nsta, 2)
ind_lon = np.arange(1, 2*nsta, 2)

lat = data[ind_lat]
lon = data[ind_lon]
for it in range(0, nsource):

   SE[it] = Cost(lat, lon, nsta, it, sources, sx, sy, deltat)  

EE = np.power(np.sum(np.power(SE,10)), 0.1)  
#EE.tofile('misfit', sep='\n')
print EE

 


