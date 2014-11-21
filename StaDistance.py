#!/usr/bin/env python 
#  this program reads stations position from an input file and calculate the minimum distance between stations.

import numpy as np
from pyrocko import  orthodrome


def all_dist(lats, lons):
    nl = len(lats)
    dis = np.zeros([nl,nl])
    for ind in  range(0, nl):  
        dis[ind] = orthodrome.distance_accurate50m_numpy( lats, lons, lats[ind], lon[ind]) 
    return dis
    
    
fname = 'vogtland.slist' 
f = open(fname)
lines = f.readlines()
f.close()

length = len(lines)
lat = []
lon = []
for ind in range(1,length):
    lat.append(lines[ind][10:19])
    lon.append(lines[ind][20:29])

lat = map(float, lat)
lon = map(float, lon)
lats = np.asarray(lat)
lons = np.asarray(lon)
distance = all_dist(lats, lons)
distance = distance[np.nonzero(np.tril(distance))]
min_distance = distance.min() # minimum distance in meter
