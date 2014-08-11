import numpy as np
from pyrocko import gf, io, util, trace, pile, model, orthodrome
import matplotlib.pyplot as plt
import math
from matplotlib.patches import Circle, Arrow
from utilities import *

def RMask(shape,centre,radius,P):
    """
    Return a boolean mask for a circular sector. 
    
    """

    x,y = np.ogrid[:shape[0],:shape[1]]
    cx,cy = centre
 
    # convert cartesian --> polar coordinates
    r2 = (x - cx)**2 + ( y - cy)**2
     
  

    # circular mask
    #print "radius*radius",radius*radius
    circmask = r2 <= radius*radius 


 
    return circmask


def TrTimeMinMax(traces): 
    ''' 
    This Function returns minimum of initial time and maximum of ending time of all traces.
    The outputs are used to extend the traces considering the common initial and ending time
    '''          
    for i in range(0, len(traces), 1):
        a = traces[i]
        t_min[i] = a.tmin
        t_max[i] = a.tmax
    return t_min.min(), t_max.max(), np.argsort(t_min)  
            
def Array_center(lat, lon):
    '''
    This function returns centeral station of the random array. 
    Central station is defined so that it has the minimum distance of all other stations.
    '''
    na = np.newaxis
    dis = np.sqrt((lon[na, :]-lon[:, na])**2+(lat[na, :]-lat[:, na])**2)
    all_dis = dis.sum(0)
    dd_s = all_dis.argsort()
    return dd_s[0]

def array_extend(nsta,traces_min_time, traces_max_time, traces, P_wind_i, P_wind_e):
    traces_a = []
    for i in range(0, nsta, 1):
        a = traces[i]
        a.extend(traces_min_time, traces_max_time, fillmethod='zeros')
#   a.bandpass(4,1,8, demean= True)
        t_1 = (P_wind_i[i])
        t_2 = (P_wind_e[i])
        aa_ = traces[i].copy()
        aa_am = aa_.get_ydata()
        aa_v = np.append(aa_am, 0)-np.append(0, aa_am)
        aa_.set_ydata(aa_v)
        aa_.chop(t_1, t_2)
        traces_a.append(aa_)
    return traces_a 

def array_normalized(traces_a, nsta):    
    a_max = [traces_a[i].absmax() for i in range(0, nsta)]
    amp = np.zeros([nsta])

    for i in range(0, nsta, 1):
        (ind_, am_) = a_max[i]
        amp[i] = am_

    ind_max = amp.argmax()
    time, amp_max = a_max[ind_max]
# ####################
    traces_b = []
    for i in range(0, nsta, 1):
        time, amp_max_ = a_max[i]
        b = traces_a[i]
        bb = traces_a[i].copy()
        bb_y = bb.get_ydata()
       # print "len(bb_y)", len(bb_y)
        amp_norm = bb_y/bb_y.max()
        bb.set_ydata(amp_norm)
        traces_b.append(bb)
    return traces_b

class RandLoc(object):
    def __init__(self, lmin, lmax, rsln, nsta):
        self.min = lmin
        self.max = lmax
        self.res = rsln
        self.nsta = nsta

    def RandomStation(self):
        self.InLoc = np.linspace(self.min, self.max, self.res)
        Pd_ = np.ones([ len(self.InLoc)]) 
        #n1 = len(self.InLoc)/4.0
       # Pd_[int(n1):-int(1*n1)] = 0.0
      #  Pd_[25:30]=0.0
      #  Pd_[50:60]=0.0     
        self.PbL = np.cumsum(Pd_)
        self.ul = self.PbL.max()
        RandSta = np.zeros(self.nsta)
        for h in range(0, self.nsta):
            Inx = np.random.randint(1, self.ul, 1)
           # print self.ul, Inx
            C = (self.PbL==Inx)
          
            d_ = np.where(C==True)
           
            RndLon = self.InLoc[d_[0]]
            RandSta[h]=RndLon[0]
        #plt.plot(self.InLoc, self.PbL)
        #plt.show()
        return RandSta 
                 

 
def ObjectFunction(Power, IndexMaxX, IndexMaxY):
    P = Power.copy()
    mask = RMask(Power.shape, (IndexMaxX, IndexMaxY), 1 ,Power)
    Power[~mask] = 0.0
    E0 = P - Power 
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
    E = np.power(np.power(E0,2).sum(), 1/2.0)
    
    
    return E
    



    
lat_min, lat_max, lon_min, lon_max, nsource, nsta = 10, 30, 20, 40, 1, 10
lat_source, lon_source, depth_source = 0., 80., 20000.0

a = np.array([1, 30000])
beam = np.array([1, 30000])
t_min = np.zeros([nsta])
t_max = np.zeros([nsta])

sources = [
    gf.DCSource(
        time=util.str_to_time('2008-02-17 11:06:01.10'), lat=lat_source,
        lon=lon_source, depth=depth_source, strike=20., dip=10., rake=90,
        magnitude=6.8) for i in range(nsource)]
#########################################################################

#lat, lon = array_config(nsta, 5, 100000., lat_min, lat_max, lon_min, lon_max)
#lat, lon, Prob, Lonsssss = RandomStation(lat_min, lat_max, lon_min, lon_max)
E = np.zeros([1000])
for it in range(0,100):
    XX = RandLoc(lat_min, lat_max, 100, nsta)
    lat = XX.RandomStation()

    YY = RandLoc(lon_min, lon_max, 100, nsta)
    lon = YY.RandomStation()
  

    targets = [
        gf.Target(
            codes=('', 'STA%03i' % i, '%03i' % i, 'Z'), lat=lat[i], lon=lon[i],
            store_id='global_2s') for i in range(nsta)]

    stations = [
        model.Station(t.codes[0], t.codes[1], t.codes[2], t.lat, t.lon)
        for t in targets]

    s = sources[0]
    ev = model.Event(time=s.time, lat=s.lat, lon=s.lon)

#print ev

#for sta in stations:
#    print sta
#    print orthodrome.distance_accurate50m(ev, sta)


    engine = gf.LocalEngine(store_dirs=['global_2s'])
    store = engine.get_store('global_2s')
    response = engine.process(sources, targets)
    traces = response.pyrocko_traces()

# ###FIRSTtrace.snuffle(traces, stations=stations, events=[ev])
# ## adjust the length of traces


    traces_min_time, traces_max_time, order_trace = TrTimeMinMax(traces)
    center = Array_center(lat, lon)

    source_reciever_dis = orthodrome.distance_accurate50m_numpy(
        lat_source, lon_source, lat, lon)
    
    t_p = np.array([
        store.t("first(p|P)", (depth_source, int(source_reciever_dis[i])))
        for i in range(0, nsta)])
    
    t_origin = util.str_to_time('2008-02-17 11:06:01.10')
    P_wind_i = t_origin + t_p - 2.
    P_wind_e = t_origin + t_p + 8.
    t_o = -t_p + t_p[center]
    
#t_s = np.array([
#    store.t("first(s|S)", (depth_source, int(source_reciever_dis[i])))
#    for i in range(0, nsta)])

#S_wind_i = t_origin + t_s - 2.
#S_wind_e = t_origin + t_s + 18.
#t_o_s = -t_s + t_s[center]



    deltat = 10.0/2.
    traces_a = array_extend(nsta,traces_min_time, traces_max_time, traces, P_wind_i, P_wind_e)
#traces_s = array_extend(nsta,traces_min_time, traces_max_time, traces, S_wind_i, S_wind_e)

# trace.snuffle(traces_a + traces)

    a_max = [traces_a[i].absmax() for i in range(0, nsta)]
    amp = np.zeros([nsta])

    for i in range(0, nsta, 1):
        (ind_, am_) = a_max[i]
        amp[i] = am_

    ind_max = amp.argmax()
    time, amp_max = a_max[ind_max]
# ####################
    traces_b = []
    for i in range(0, nsta, 1):
        time, amp_max_ = a_max[i]
        b = traces_a[i]
        bb = traces_a[i].copy()
        bb_y = bb.get_ydata()
   # print "L=",len(bb_y), bb_y.max()
        amp_norm = bb_y/bb_y.max()
        bb.set_ydata(amp_norm)
    #print "A=",len(bb_y), amp_norm.max()
        traces_b.append(bb)
# trace.snuffle(traces_b + traces_a)
    sy = np.linspace(-0.0001, 0.0001, 50)
    sx = np.linspace(-0.0001, 0.0001, 50)
    slowness = sx_sy(lat, lon, center, t_o)


#####################################################################
    P, delay = ArrayBeamForming(
        sx, sy, traces_b, P_wind_i[center], P_wind_e[center], len(traces), center, lat, lon, t_o, deltat)

    P = P.transpose()

    BeamShow = AbeamShow(power=P,sx=sx, sy=sy)
    IndexMaxX, IndexMaxY = BeamShow.MaxPower()
    BeamShow.PowerQuality(5, 10)
   

    Quality = BeamShow.Ql
    BeamShow.Draw(it)
    BeamShow.DrawCircle(1e-5)
    BeamShow.DrawArow()
    plt.savefig('power%s.png' %it)
    
    dum3 = trace.Trace(
        '', '', 'Beam', '', tmin=traces_min_time, deltat=traces[0].deltat,
        ydata=np.zeros(int((traces_max_time - traces_min_time)/traces[0].deltat)))
    
    traces_d, dela = array_beam_s(center, lat, lon, sy[IndexMaxY], sx[IndexMaxX], traces, dum3)
#trace.snuffle(traces_d )



    E[it] = ObjectFunction(P, IndexMaxX, IndexMaxY)
    #plt.savefig('beamforming3.png')


    plt.figure(it*500)
    plt.plot(lon, lat, '*', lon_source, lat_source, 'o')
    plt.text(lon[center], lat[center], 'center')
    plt.xlabel('lon(degree)')
    plt.ylabel('lat(degree)')
    
    
    #plt.savefig('station_event3.png')
    plt.xlim([min(lon_min, lon_source), max(lon_max, lon_source)])
    plt.ylim([min(lat_min, lat_source), max(lat_max, lat_source)])
    plt.savefig('station%s.png' %it)


#plt.plot(E)
#plt.show()
   # plt.show()
   # plt.show()

