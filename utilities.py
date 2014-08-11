import numpy as np
from pyrocko import gf, io, util, trace, pile, model, orthodrome
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arrow
import datetime as dt
class Abeam(trace.Trace):
    def __init__(self, trc, delay, tinit, tend):
# delay_v calculated time delay for all traces
        self.tr = trc
        self.td = delay
        self.start = tinit
        self.tend = tend
        self.x = len(self.tr)
        self.delt = trc[0].deltat

    def trace_beam(self):
        c, b = self.td.shape
        P = np.zeros(c)
        for j in range(1,c):
   
             dum = trace.Trace('','','Beam','', tmin= self.start, deltat=self.delt,
                 ydata=np.zeros(int((self.tend - self.start)/self.delt))) 
             #print self.td[j]
             for i in range(1,self.x):
                   tt= self.td[j][i]
                   dd = self.tr[i].copy()   
                   dd.shift(tt)
                   dum.add(dd)
                  
             dum.chop( self.start , self.tend)
              
             P[j]= sum((dum.get_ydata())**2)/self.x 
             

      #  print  "inside = ", P.max() 
        P = P/P.max()
      #  print  P.max()
        return P
            
#####################################################################################
def array_config( n, aperture, min_distance, lat_min, lat_max, lon_min, lon_max): 
   
   rlat = np.random.uniform(lat_min, lat_max,1)
   lat = np.around(rlat, decimals=1)
   rlon = np.random.uniform(lon_min, lon_max, 1)
   lon = np.around(rlon, decimals=1)
   locs=np.array([lat[0], lon[0]])
   lat2=np.zeros([1])
   lon2= np.zeros([1])
  
   dum=0.0
   while  dum==0.0  :
           if lat2[0] <= lat_max and lat2[0] >= lat_min and lon2[0] <= lon_max  and lon2[0] >= lon_min:
              dum=1.0
              new_loc = [lat2[0],lon2[0]]
              Loc=np.vstack((locs,new_loc)) 
    
           else: 
              rphi = np.random.uniform(0, 2*np.pi ,1)
              x = aperture * np.sin(rphi)
              y = aperture * np.cos(rphi)
              lat2 = lat + x
              lon2 = lon + y
  



   Loc_f = Loc   
   while len(Loc_f) <= n - 1: 
         
         rlat = np.random.uniform(lat_min, lat_max)
         lat_temp = np.around(rlat, decimals=1)
         rlon = np.random.uniform(lon_min, lon_max)
         lon_temp = np.around(rlon, decimals=1)
         new_loc = [rlat,rlon]
         Loc=np.vstack((Loc,new_loc)) 
         Loc_f=distance_check(Loc, min_distance)

   return Loc_f[:, 0], Loc_f[:, 1] 
   
def distance_check(Loc, min_distance):
   n=len(Loc)
   dis=[]
   check = 0
   lat=Loc[:,0]
   lon=Loc[:,1]
   na = np.newaxis    
   dis = orthodrome.distance_accurate50m_numpy(lat[na,:], lat[:,na] , lon[na,:], lon[:,na])
   compar = np.ones([n,n])*min_distance
   C= (dis >= compar)
   d,d_=np.where(C==False)
   check = d.size
   
            
             
   if check==0:
 
       Loc_ = Loc
    
   else:
   
       Loc_ = Loc[0:-1]
             
   return Loc_

########################################
def ArrayBeamForming(sx,sy,traces, t_i_c ,t_e_c, nsta , center, lat, lon, t_o, deltat):
   PP = np.zeros([sx.size *sy.size])
   P = np.zeros([sx.size,sy.size])
   na = np.newaxis   
   dis_x = orthodrome.distance_accurate50m_numpy( lat[center], lon , lat[center], lon[center])
   dis_y = orthodrome.distance_accurate50m_numpy( lat, lon[center], lat[center], lon[center])
   dis_x[center]=0.0
   dis_y[center]=0.0
   dum_ = np.where( lat > lat[center] )
   dis_y[dum_]= -dis_y[dum_] 
            
   dum = np.where(lon  > lon[center])   
   dis_x[dum] = -1.0*dis_x[dum]
   X = np.array([dis_x , dis_y])
   S =np.transpose([np.tile(sx, len(sy)), np.repeat(sy, len(sx))])
   delay = np.dot(S,X)
               
   Q = Abeam(traces, delay, t_i_c , t_e_c)              
   PP = Q.trace_beam()  
  
   P=PP.reshape(sx.size, sy.size)
   P=P.transpose()
   return P, delay


###################


     
###################################
def sx_sy(lat, lon, center, t):  

     n = len(lat)    
     d_x = np.zeros([n])
     d_y = np.zeros([n])
     sy = np.zeros([n])
     sx = np.zeros([n])

     na = np.newaxis    
     dis_x = orthodrome.distance_accurate50m_numpy( lat[center], lon , lat[center], lon[center])
     dis_y = orthodrome.distance_accurate50m_numpy( lat, lon[center], lat[center], lon[center])
     dis_x[center]=0.0
     dis_y[center]=0.0
     dum_ = np.where( lat > lat[center] )
     dis_y[dum_]= -dis_y[dum_] 
            
     dum = np.where(lon  > lon[center])   
     dis_x[dum] = -1.0*dis_x[dum]
   
     GG = np.array([dis_x,dis_y])
     G = np.transpose(GG) 
     slowness_vector = la.lstsq(G, t)[0]   
   
     return slowness_vector


####################################################################     
def array_beam_s(center,lat,lon,Sx,Sy, traces, dum_tr):  

   nsta = len(traces)
   delay = np.zeros(nsta)
   traces_d = []
   for st in range(0,nsta,1):
      
      if st != center :
           dis_x = orthodrome.distance_accurate50m_numpy( lat[center], lon[st] , lat[center], lon[center])
           dis_y = orthodrome.distance_accurate50m_numpy( lat[st], lon[center], lat[center], lon[center])
           if lat[st] > lat[center] :
               dis_y = -1.0*dis_y
            
           if  lon[st] > lon[center]  :   
               dis_x = -1.0*dis_x
       
           delay[st] = dis_x * Sx + dis_y * Sy  
         
      else:
           delay[st] = 0.0
   
           
      dd = traces[st].copy() 
      dd.shift(delay[st])
      dum_tr.add(dd)
      traces_d.append(dd)
   traces_d.append(dum_tr) 
   return traces_d, delay
   
   
#####################################################################################   
def SectorMask(shape,centre,radius,angle_range,P):
    """
    Return a boolean mask for a circular sector. The start/stop angles in  
    `angle_range` should be given in clockwise order.
    """

    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = centre
    tmin,tmax = np.deg2rad(angle_range)

    # ensure stop angle > start angle
    if tmax < tmin:
            tmax += 2*np.pi

    # convert cartesian --> polar coordinates
    r2 = (P[x,cy]-P[cx,cy])**2 + (P[cx,y]-P[cx,cy])**2
    
    #print 'r2', r2
    theta = np.arctan2(x-cx,y-cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2*np.pi)

    # circular mask
    #print "radius*radius",radius*radius
    circmask = r2 <= radius*radius

    # angular mask
    anglemask = theta <= (tmax-tmin)
 
    return circmask*anglemask
    

# ############################################################################################    

class AbeamShow(object):
    def __init__(self, power, sx, sy):
        self.power = power
        self.sx = sx
        self.sy = sy
        
    def MaxPower(self):
        self.IndexMaxX, self.IndexMaxY = np.unravel_index(self.power.argmax(), self.power.shape)
        return   self.IndexMaxX,self.IndexMaxY 
          
    def PowerQuality(self, nradi, nphi):
        self.nradi = nradi
        self.nphi = nphi
        Quality = np.zeros([nradi, nphi])
        for j in range(0, nphi):
            for i in range(0, nradi):
                delt_t = 360.0/nphi
                Pdum = self.power.copy()
                mask = SectorMask(Pdum.shape, (self.IndexMaxX, self.IndexMaxY),  0.5*j , (0 + (i-1)*delt_t, i*delt_t) , Pdum )
                Pdum[~mask] = 0.0
                non_zeros = Pdum[np.nonzero(Pdum)]
                nnonzeros = non_zeros.size
                Quality [i][j] = Pdum.sum()#/nnonzeros
        self.Ql = Quality
          
    def Draw(self, i):
        self.i = i
        self.fig = plt.figure(self.i)
        self.ax = self.fig.add_subplot(1, 1, 1)
        levels = np.linspace(self.power.min(), self.power.max(), len(self.sx))
        plt.contourf(
            self.sx, self.sy, self.power, levels, ls='-', origin='lower',
            cmap=plt.cm.jet)
        plt.colorbar()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xticks(np.arange(min(self.sx), max(self.sx), 0.1*1.0e-4))
        plt.yticks(np.arange(min(self.sy), max(self.sx), 0.1*1.0e-4))
        plt.xlabel('sx')
        plt.ylabel('sy')
        
    def DrawCircle(self, DelRad): 
        self.rad = DelRad 
        self.rad2 = DelRad
        
        while self.rad < self.sx[-1]:
          #  print self.rad
            circle = Circle( (0, 0), radius= self.rad, linewidth=2, facecolor='none',
            edgecolor=(0, 0, 0), alpha=0.1)
            self.ax.add_patch(circle)
            self.rad = self.rad + self.rad2 
            
    def DrawArow(self):
        Arrow2 = Arrow(0., 0.0, self.sy[self.IndexMaxY], self.sx[self.IndexMaxX], width=0.1e-4, alpha=0.4)
        self.ax.add_patch(Arrow2)   
        
         
            

              
 





