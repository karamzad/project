import numpy as np
from pyrocko import gf, io, util, trace, pile, model, orthodrome
import numpy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arrow
import time
from pyrocko.parstack import parstack

class Abeam(trace.Trace):
    '''
    Abeam is an object defined based pyrocko trace object.
    this object added trace_beam method to pyrocko trace.
    trace_beam shifts the traces based on given delay time and allighn
    '''
    def __init__(self, trc, delay, tinit, tend):
        self.tr = trc
        self.td = delay
        self.start = tinit
        self.tend = tend
        self.x = np.size(self.tr)
        self.delt = trc[0].deltat

    def trace_beam(self):
        '''
        This method returns traces beam power. 
        '''
        
        c, b = self.td.shape
        P = np.zeros(c)
        for j in range(1,c): 
             dum = trace.Trace('', '', 'Beam', '', tmin= self.start, deltat=self.delt,
                 ydata=np.zeros(int((self.tend - self.start)/self.delt))) 
            
             for i in range(1, self.x):
                   tt = self.td[j][i]
                   dd = self.tr[i].copy()   
                   dd.shift(tt)
                   dum.add(dd)
                  
             dum.chop( self.start, self.tend)
             P[j]= sum((dum.get_ydata())**2)/self.x 

        P = P/P.max()
        return P

def abeam(traces, delays, tmin, tmax):
    narrays = len(traces)
    
   # for tr in traces:
    #    tr.lowpass(4, 0.1)
  #  trace.snuffle(traces)
        
    tr_tmins = np.array([tr.tmin for tr in traces], dtype=np.float)
    t0 = tr_tmins[0]
    deltat = traces[0].deltat

    offsets = np.round((tr_tmins - t0) / deltat).astype(np.int32)
    offsets = np.concatenate([offsets, offsets])
    
    nshifts = delays.shape[0]
    
    shifts = np.zeros((nshifts, narrays*2), dtype=np.int32)
    
    rshifts = delays / deltat
    fshifts = shifts[:,:narrays] = np.floor(rshifts)
    cshifts = shifts[:, narrays:] = np.ceil(rshifts)
    
    assert shifts.shape[1] == narrays*2

    arrays = [tr.get_ydata().astype(np.float) for tr in (traces + traces)]
    weights = np.zeros(shifts.shape, dtype=np.float)
    weights[:,:narrays] = 1.0 - (rshifts - fshifts)
    weights[:, narrays:] = (1.0 - (cshifts - rshifts)) * (cshifts - fshifts)
    
    p, _ = parstack(arrays, offsets, shifts, weights, 1)

    p = p/p.max()
    return p

            
#####################################################################################


def ArrayBeamForming(sx, sy, traces, t_i_c, t_e_c, nsta , center, lat, lon, t_o, deltat):
   '''

   '''
   PP = np.zeros([sx.size *sy.size])
   P = np.zeros([sx.size, sy.size])
      
   dis_x = orthodrome.distance_accurate50m_numpy( lat[center], lon, lat[center], lon[center])
   dis_y = orthodrome.distance_accurate50m_numpy( lat, lon[center], lat[center], lon[center])
   print "dis_x[center]",dis_x[center]
   dis_x[center]=0.0
   dis_y[center]=0.0
   dum_ = np.where( lat > lat[center])
   dis_y[dum_] = - 1.0* dis_y[dum_] 
            
   dum = np.where(lon  > lon[center])   
   dis_x[dum] = - 1.0*dis_x[dum]
   X = np.array([dis_x, dis_y])
   S = np.transpose([np.tile(sx, np.size(sy)), np.repeat(sy, np.size(sx))])
   delay = np.dot(S,X)
#   t0 = time.time()            
   PP = abeam(traces, delay, t_i_c, t_e_c)
#   t1 = time.time()
#   Q = Abeam(traces, delay, t_i_c, t_e_c)     
#   PP = Q.trace_beam()
 #  print PP.dtype
 #  t2 = time.time()  
#   print t1 - t0 
   P = PP.reshape(sx.size, sy.size)
   P = P.transpose()
   return P, delay


###################
   
   
#####################################################################################   
def SectorMask(shape, centre, radius, angle_range, P):
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
        return   self.IndexMaxX, self.IndexMaxY 
          
    def PowerQuality(self, nradi, nphi):
        self.nradi = nradi
        self.nphi = nphi
        Quality = np.zeros([nradi, nphi])
        for j in range(0, nphi):
            for i in range(0, nradi):
                delt_t = 360.0/nphi
                Pdum = self.power.copy()
                mask = SectorMask(Pdum.shape, (self.IndexMaxX, self.IndexMaxY), 0.5*j, (0 + (i-1)*delt_t, i*delt_t), Pdum )
                Pdum[~mask] = 0.0
                non_zeros = Pdum[np.nonzero(Pdum)]
                nnonzeros = non_zeros.size
                Quality [i][j] = Pdum.sum()#/nnonzeros
        self.Ql = Quality
          
    def Draw(self, i):
        self.i = i
        self.fig = plt.figure(self.i, figsize=(20., 14))
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.baz = 270.0-np.degrees(np.arctan2(self.sy[self.IndexMaxX], self.sx[self.IndexMaxY]))
        levels = np.linspace(self.power.min(), self.power.max(), np.size(self.sx))
        plt.contourf(
            self.sx, self.sy, self.power, levels, ls='-', origin='lower',
            cmap=plt.cm.jet)
           
        plt.colorbar()
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        plt.xticks(np.arange(min(self.sx), max(self.sx), 0.1*1.0e-4))
        plt.yticks(np.arange(min(self.sy), max(self.sx), 0.1*1.0e-4))
        self.ax.set_title('Array Beam Power')
        self.ax.annotate('sx: %3.2e         sy: %3.2e        BAZ: %s' %(self.sx[self.IndexMaxY],
                 self.sy[self.IndexMaxX], self.baz), xy=(1, 0), xycoords='axes fraction',
                 fontsize=13, xytext=(-20, -30), textcoords='offset points',
                ha='right', va='top')
        plt.xlabel('sx')
        plt.ylabel('sy')
        
    def DrawCircle(self, DelRad): 
        self.rad = DelRad 
        self.rad2 = DelRad
        
        while self.rad < self.sx[-1]:
          #  print self.rad
            circle = Circle((0, 0), radius= self.rad, linewidth=7, facecolor='none',
            edgecolor=(0, 0, 0), alpha=0.1)
            self.ax.text(0, self.rad, '%s' %self.rad, fontsize=14)
            self.ax.text(self.rad, 0,'%s' %self.rad, fontsize=14)
            self.ax.text(-self.rad, 0,'%s' %self.rad, fontsize=14)
            self.ax.text(0, -self.rad ,'%s' %self.rad, fontsize=14)
            self.ax.add_patch(circle)
            self.rad = self.rad + self.rad2 
            
    def DrawArow(self):
        Arrow2 = Arrow(0., 0.0, self.sy[self.IndexMaxY], self.sx[self.IndexMaxX], width=0.1e-4, alpha=0.4)
        self.ax.add_patch(Arrow2) 
        
    def saveplot(self,num):
        self.num = num
        plt.savefig('Power%s.png' %self.num)  
        
         
def window(lat, lon, lat_source, lon_source, depth_source, nsta, store): 

        center = Array_center(lat, lon)
        source_reciever_dis = orthodrome.distance_accurate50m_numpy(
            lat_source, lon_source, lat, lon)
    
        t_p = np.array([
            store.t("first(p|P)", (depth_source, int(source_reciever_dis[i])))
            for i in range(0, nsta)])
        t_s = np.array([
           store.t("first(s|S)", (depth_source, int(source_reciever_dis[i])))
         for i in range(0, nsta)])
       
        t_origin = util.str_to_time('2008-02-17 11:06:01.10')
       
        def win_(t_, t_l , t_r, center):
            wind_i = t_origin + t_ - t_l
            wind_e = t_origin + t_ + t_r
            t_o = - t_ + t_[center]
            return wind_i, wind_e, t_o
        
        P_wind_i, P_wind_e, t_op = win_(t_p, 5.0 , 20.0, center)
        S_wind_i, S_wind_e, t_os = win_(t_s, 2.0 , 18.0, center)
        
        return P_wind_i, P_wind_e, t_op , center#, S_wind_i, S_wind_e, t_os
        
def Array_center(lat, lon):
    '''
    This function returns centeral station of the randomly produced array. 
    Central station is defined so that it has the minimum distance of all other stations.
    '''
    na = np.newaxis
    dis = np.sqrt((lon[na, :]-lon[:, na])**2 + (lat[na, :] - lat[:, na])**2)
    all_dis = dis.sum(0)
    dd_s = all_dis.argsort()
    return dd_s[0]
 
def TrTimeMinMax(traces): 
    ''' 
    This Function returns minimum of initial time and maximum of ending time of all traces.
    The outputs are used to extend the traces considering the common initial and ending time
    '''
    ns = np.size(traces)  
    t_min = np.zeros([ns])
    t_max = np.zeros([ns])        
    for i in range(0, ns, 1):
        a = traces[i]
        t_min[i] = a.tmin
        t_max[i] = a.tmax
    return t_min.min(), t_max.max(), np.argsort(t_min)  
            


def array_extend(nsta, traces_min_time, traces_max_time, traces, P_wind_i, P_wind_e):
    '''
    Return the extended traces with common time interval (same start and end).
    '''
    traces_a = []
    for i in range(0, nsta, 1):
        a = traces[i]
        a.extend(traces_min_time, traces_max_time, fillmethod='zeros')
    
        a.lowpass(4, 0.1)
      
#   a.bandpass(4,1,8, demean= True)
        t_1 = (P_wind_i[i])
        t_2 = (P_wind_e[i])
        aa_ = traces[i].copy()
        aa_am = aa_.get_ydata()
       # aa_v = np.append(aa_am, 0)- np.append(0, aa_am)
        aa_.set_ydata(aa_am)
        aa_.chop(t_1, t_2)
        traces_a.append(aa_)
    return traces_a 

def array_normalized(traces_a_orig, nsta): 
    '''
    Return the normalized trace based on maximum amplitude in Phase's time window.
    '''  
    traces_a = []
    for i in range(0, nsta, 1):
        a = traces_a_orig[i].copy()
        a_am = a.get_ydata()
        a_v = np.append(a_am, 0)- np.append(0, a_am)
        aa_ = traces_a_orig[i].copy()
        aa_.set_ydata(a_v)
        traces_a.append(aa_)
        

    a_max = [traces_a[i].absmax() for i in range(0, nsta)]
    amp = np.zeros([nsta])

    for i in range(0, nsta, 1):
        (ind_, am_) = a_max[i]
        amp[i] = am_

    ind_max = amp.argmax()
    time, amp_max = a_max[ind_max]

    traces_b = []
    for i in range(0, nsta, 1):
        time, amp_max_ = a_max[i]
        b = traces_a[i]
        bb = traces_a[i].copy()
        bb_y = bb.get_ydata()
        amp_norm = bb_y/bb_y.max()
        bb.set_ydata(amp_norm)
        traces_b.append(bb)
    return traces_b


g_engine = gf.LocalEngine(store_dirs=['../global_2s'])

def GTrace(targets, s ):
       
        ev = model.Event(time=s.time, lat=s.lat, lon=s.lon)
#print ev
#for sta in stations:
#    print sta
#    print orthodrome.distance_accurate50m(ev, sta)

        store = g_engine.get_store('global_2s')
        response = g_engine.process(s, targets)
        traces = response.pyrocko_traces()
        return traces, store
        
        

def RMask(shape, centre, radius, P):
    """
    Return a boolean mask for a circular sector. 
    
    """
    x, y = np.ogrid[:shape[0], :shape[1]]
    cx, cy = centre
    # convert cartesian --> polar coordinates
    r2 = np.power((x - cx), 2) + np.power(( y - cy), 2)
    # circular mask
    circmask = r2 <= np.power(radius, 2) 
    return circmask

           

 
def ObjectFunction(Power, IndexMaxX, IndexMaxY):
    '''
    Return the value of object function.
    '''
    P1 = Power.copy()
    mask = RMask(Power.shape, (IndexMaxX, IndexMaxY), 3 ,Power)
    Power[~mask] = 0.0
    E0 = P1 - Power 
    Q = E0 <= 0.2
    E0[Q]=0
#    plt.figure(6)
#    levels = np.linspace(P.min(), P.max(), len(sx))
#    plt.contourf(
#            sx, sy, E0  , levels, ls='-', origin='lower',
#            cmap=plt.cm.jet)
#  plt.colorbar()
#plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
#plt.xticks(np.arange(min(sx), max(sx), 0.1*1.0e-4))
#plt.yticks(np.arange(min(sy), max(sx), 0.1*1.0e-4))
#plt.xlabel('sx')
#plt.ylabel('sy')
   
    #E = np.sqrt(np.sum(np.power(E0,2)))
    E = np.power(np.power(E0,3).sum(), 1/3.0) 
    return E
    
class RandomSource(object):
    '''
    Return the location of random sources. 
    After shocks.
    Input parameters:
    M : mainshocke magnitude.
    loc : lat, lon
    focal mechanism: Strike, Dip, Rake
    '''
    def __init__(self, lat, lon, depth, mag, strike, dip, rake, delta_lat, delta_lon, delta_dep, nsource):
        self.lat = lat
        self.lon = lon
        self.mg = mag
        self.dpth = depth
        self.strk = strike
        self.dp = dip
        self.rk = rake    
        self.dlat = delta_lat
        self.dlon = delta_lon
        self.ddep = delta_dep
        self.n = nsource 
 
        
    def RandomS(self):
        self.lats = RanRon(self.lat, self.dlat, 2, self.n)
        self.lons = RanRon(self.lon, self.dlon, 2, self.n)
        self.dpths = RanRon(self.dpth, self.ddep, 1, self.n)
        self.strks = RanRon(self.strk - 30.0, self.strk + 30.0, 0, self.n)
        self.dps = RanRon(self.dp - 30.0, self.dp + 30.0, 0, self.n)
        self.rks = RanRon(self.rk - 30.0, self.rk + 30.0, 0, self.n)
        self.mgs = RanRon(1, self.mg, 5, self.n)

        sources = [
            gf.DCSource(
            time = util.str_to_time('2008-02-17 11:06:01.10'),
            lat = self.lats[i],
            lon = self.lons[i],
            depth = self.dpths[i],
            strike = self.strks[i],
            dip =  self.dps[i],
            rake =  self.rks[i],
            magnitude = self.mgs[i] ) for i in range(self.n)]
        return sources 
        
        
def RanRon( init, delta, decim, n):
        init = init
        delta = delta
        decim = decim
        Rvalue = np.around(np.random.uniform(init-delta, init+delta, n))
        Rvalue = np.around(Rvalue, decimals=decim)
        return Rvalue       

def array_beam_s(center, lat, lon, Sx, Sy, traces, dum_tr):  

   nsta = np.size(traces)
   delay = np.zeros(nsta)
   traces_d = []
   for st in range(0,nsta,1):
      
      if st != center :
           dis_x = orthodrome.distance_accurate50m_numpy( lat[center], lon[st], lat[center], lon[center])
           dis_y = orthodrome.distance_accurate50m_numpy( lat[st], lon[center], lat[center], lon[center])
           if lat[st] > lat[center]:
               dis_y = -1.0*dis_y
            
           if  lon[st] > lon[center]:   
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

def AperMin(lat, lon):
    '''
    
    This function returns array aperture and minimum distance between stations.
    inputs are a set of lat and lon of stations.
     
    '''
    Ns = len(lat)                                         
    dis=np.zeros([Ns, Ns])                                                           
                                                                 
    for i in range(0,Ns):                                                           
        dis[i] = orthodrome.distance_accurate50m_numpy(lats, lons, lats[i], lons[i])                                                                   
    dis /= 1000.0                                                                   
    dis = np.around(dis, decimals=0)                                                                                                                            
                                                                                
    iu1 = np.triu_indices(Ns)                                                       
    dis[iu1]=0                                                                      
                                                                                
    ar_mindis = dis[np.nonzero(dis)].min()                                          
    ar_aperture = dis[np.nonzero(dis)].max()
    return ar_mindis, ar_aperture        
        
def Cost(lat, lon, nsta, nsource, sources, sx, sy, deltat):

    targets = [
        gf.Target(
            codes=('', 'STA%03i' % i, '%03i' % i, 'Z'), lat=lat[i], lon=lon[i],
            store_id='global_2s') for i in range(nsta)]
    ind=np.array([nsource])
    for it_source in ind:
        
        s = sources[it_source]
##        print s.lat, s.lon
        traces, store = GTrace(targets, s)
        traces_min_time, traces_max_time, order_trace = TrTimeMinMax(traces)
        P_wind_i, P_wind_e, t_op, center = window(lat, lon, s.lat, s.lon, s.depth, nsta, store)
            
        traces_a = array_extend(nsta, traces_min_time, traces_max_time, traces, P_wind_i, P_wind_e)
#traces_s = array_extend(nsta,traces_min_time, traces_max_time, traces, S_wind_i, S_wind_e)
# trace.snuffle(traces_a + traces)
        traces_b = array_normalized(traces_a, nsta)
         # trace.snuffle(traces_b + traces_a)
        P_, delay = ArrayBeamForming(
            sx, sy, traces_b, P_wind_i[center], P_wind_e[center], nsta, center, lat, lon, t_op, deltat)
##        print 'center = ', center
        P_ = P_.transpose()
        BeamShow = AbeamShow(power=P_, sx=sx, sy=sy)
        IndexMaxX, IndexMaxY = BeamShow.MaxPower()
     
  
 #       BeamShow.PowerQuality(5, 10) 
 #       Quality = BeamShow.Ql
 #       BeamShow.Draw(10)
 #       BeamShow.DrawCircle(1e-5)
  
#      BeamShow.DrawArow()
  #      plt.show()

      # To check if the traces are inline uncomment 
      # dum3 = trace.Trace(
      #    '', '', 'Beam', '', tmin=traces_min_time, deltat=traces[0].deltat,
      #    ydata=np.zeros(int((traces_max_time - traces_min_time)/traces[0].deltat)))
    
      #  traces_d, dela = array_beam_s(center, lat, lon, sy[IndexMaxY], sx[IndexMaxX], traces, dum3)
      #  trace.snuffle(traces_d )

    return ObjectFunction(P_, IndexMaxX, IndexMaxY)
              
 





