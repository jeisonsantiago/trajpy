# -*- coding: utf-8 -*-
"""
Created on Sun Jul 31 23:21:25 2016

@author: Jeison Santiago jeison.santiago@gmail.com

This code is based on traj.gs 
a grads script used to calculate Horizontal 
Forward and Backward trajectories

ftp://cola.gmu.edu/grads/scripts/traj.gs

The algorithm is compatible with era-interim pressure levels dataset,
it's a crude work and I used only for testing
"""

#open file and plotimport numpy as np
from mpl_toolkits.basemap import Basemap
import  matplotlib.pyplot as plt
import numpy as np
import math
from netCDF4 import Dataset
import sys

#----------------------- -180 <-> 180 to 0 <-> 360 ----------------------------

# simple conversion
def to0_360(lon):
    return (lon + 180) % 360
    
def to180_180(lon):
    return lon - 180
#----------------------- -180 <-> 180 to 0 <-> 360 ----------------------------
# the net cdf file
netfile = 'filename'
ncfile = Dataset(netfile,'r') 

lat_top = ncfile.variables['latitude'][0]
lat_down = ncfile.variables['latitude'][-1]

lon_left = ncfile.variables['longitude'][0]
lon_right = ncfile.variables['longitude'][-1]

lat_size = ncfile.dimensions['latitude'].size
lon_size = ncfile.dimensions['longitude'].size

print lat_top,lat_down,lon_left, lon_right

print 'set basemap'
m = Basemap( projection='cyl',llcrnrlat=lat_down,urcrnrlat=lat_top,\
            llcrnrlon=lon_left,urcrnrlon=lon_right,resolution='c')
            
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
m.drawmapboundary(fill_color='aqua')
#plt.title("Equidistant Cylindrical Projection")

fig = plt.figure(1)

lat_data = [];
lon_data = [];

u = []
v = []

global _Xmax 
global _Xmin
global _Ymax
global _Ymin

global _At
global _PI
global _D2R
global _R2D
_At=6371229
_PI=3.141592654
_D2R=_PI/180
_R2D=180/_PI

#--------------------------get the angle ------------------------------
#wind_abs = sqrt(u_ms^2 + v_ms^2)
#wind_dir_trig_to = atan2(u_ms/wind_abs, v_ms/wind_abs) 
#wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi ## -111.6 degrees

def pointsTo(u,v,front=True):
    wind_abs = np.sqrt(u**2 + v**2)
    print 'absolute uv:%f'%(wind_abs)
    wind_dir_trig_to = math.atan2(u/wind_abs, v/wind_abs) 
    wind_dir_trig_to_degrees = wind_dir_trig_to * 180/np.pi ## -111.6 degrees
    wind_dir_trig_cardinal = 90 - wind_dir_trig_to_degrees
    #return wind_dir_trig_to
    return wind_dir_trig_cardinal
#--------------------------find nearest value ---------------------------------
def getnearpos(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx   

#-------------------------- fi function ---------------------------------------
def fi(i):
    
    global _Xmax 
    global _Xmin
    global _Ymax
    global _Ymin
    _Xmax = lonarray.size
    _Xmin = 0
    _Ymax = latarray.size
    _Ymin = 0    
    
    if(i == 1):
        print 'fi exit'
        sys.exit(1)

# --------------------------- distance function -------------------------------
def distance(lat1, lng1, lat2, lng2):
    #return distance as meter if you want km distance, remove "* 1000"
    radius = 6371 * 1000 

    dLat = (lat2-lat1) * math.pi / 180
    dLng = (lng2-lng1) * math.pi / 180

    lat1 = lat1 * math.pi / 180
    lat2 = lat2 * math.pi / 180

    val = math.sin(dLat/2) * math.sin(dLat/2) + math.sin(dLng/2) * math.sin(dLng/2) * math.cos(lat1) * math.cos(lat2)    
    ang = 2 * math.atan2(math.sqrt(val), math.sqrt(1-val))
    return radius * ang

def dist(lon1,lat1,lon2,lat2):

    phi=lon1*_D2R
    theta=(90-lat1)*_D2R
    x1=math.sin(theta)*math.cos(phi)
    y1=math.sin(theta)*math.sin(phi)
    z1=math.cos(theta)
    phi=lon2*_D2R
    theta=(90-lat2)*_D2R
    x2=math.sin(theta)*math.cos(phi)
    y2=math.sin(theta)*math.sin(phi)
    z2=math.cos(theta)
    x=y1*z2-y2*z1
    y=x2*z1-x1*z2
    z=x1*y2-x2*y1
    d2=x*x+y*y+z*z
    res = math.asin(math.sqrt(d2))
    return res * _At

#-------------------------- mou -----------------------------------------------
def mou(lon0, lat0, dist, alpha):
# Coordinates of a point located a distance dist (in meters) and angle
# alpha (in degrees and mathematical convention) from the point (lon0,lat0).
    if (dist<1):
        return lon0,lat0
        
    lon0=lon0*_D2R
    lat0=lat0*_D2R
    
    alpha=90-alpha
    
    if(alpha < 0):
        alpha=360+alpha

    if(alpha > 360):
        alpha=alpha-360
        
    A=alpha*_D2R
    b=dist/_At
    c=_PI/2-lat0
    
    #acos(cos('b')*cos('c')+sin('b')*sin('c')*cos('A'))'
    a = math.acos(math.cos(b)*math.cos(c)+math.sin(b)*math.sin(c)*math.cos(A))
    #lat1= ((_PI/2-a)*_R2D)-180
    lat1= (_PI/2-a)*_R2D
    #asin(sin('b')*sin('A')/sin('a'))'
        
    B = math.asin(math.sin(b)*math.sin(A)/math.sin(a))
    lon1 = (lon0+B)*_R2D
    
    return lon1,lat1

#-------------------------- interpolate ---------------------------------------
def interp(camp,LONt,LATt):
    
    #get points lonlat to xy
    xin = getnearpos(lonarray,LONt)
    yin = getnearpos(latarray,LATt)

    lon = np.zeros(4)
    lat = np.zeros(4)
    
    z = np.zeros(4)
    d = np.zeros(4)

    x = np.zeros(4,dtype=int)
    y = np.zeros(4,dtype=int)
    
    x[1] = xin
    x[2] = x[1]
    x[0] = x[1]+1
    x[3] = x[0]

    y[2] = yin
    y[3] = y[2]
    y[0] = y[2]+1
    y[1] = y[0]
    
    num = 0; den = 0
    
    #print 'x',x
    #print 'y',y
    
    i = 0
    while(i<4):
        if( x[i] < _Xmin or x[i] > _Xmax or y[i] < _Ymin or y[i] > _Ymax):
            print 'x:',x[i],' xmin:',_Xmin,' xmax:',_Xmax
            print 'y:',y[i],' ymin:',_Ymin,' ymax:',_Ymax
            #print 'lon:',lonarray[x[i]],to180_180(lonarray[x[i]])
            #print 'lat:',latarray[y[i]]            
            fi(1)

        #dev a IndexError try catch !!!  
        #print '_x:',x[i],' xmin:',_Xmin,' xmax:',_Xmax
        #print '_y:',y[i],' ymin:',_Ymin,' ymax:',_Ymax
         
        try:
            lon[i] = lonarray[x[i]]
            lat[i] = latarray[y[i]]                    
            z[i] = camp[y[i],[x[i]]]
        except IndexError:
            print 'out of data bounds'
            lon[i] = lonarray[-1]
            lat[i] = latarray[-1]
            z[i] = camp[-1,-1]
         
        
        #update camp???
        
        
        # distance
        #d[i] = distance(LATt,LONt,lat[i],lon[i])
        d[i] = dist(LONt,LATt,lon[i],lat[i])
        
        if(d[i] == 0):
            return z[i]
            
        num=num+1/d[i]*z[i]
        den=den+1/d[i]
        
        i=i+1
    
    return (num/den)
#--------------------------def what to do--------------------------------------
def afterLatLon(lon, lat):
    global ncfile
    print('xdata=%f, ydata=%f' %(lon, lat))
    #print ncfile.variables
    print '--------------------------------------'
    
       
    #convert lon 
    lon = to0_360(lon)    
    print 'pontos_init:',lon,lat    
    
    #get x and y
    global lonarray
    lonarray = ncfile.variables['longitude'][:]
    
    #test
    #lonarray = lonarray+180
    for i in range(0,lonarray.size):
        lonarray[i] = to0_360(lonarray[i])
    
    global latarray
    latarray = ncfile.variables['latitude'][:]
    
    #xlon = getnearpos(lonarray,lon)
    #ylat = getnearpos(latarray,lat)
    
    #print('xdata=%f, ydata=%f x=%d y=%d' %(lon, lat, xlon, ylat))
    #print('xdata_real=%f, ydata_real=%f' %( lonarray[xlon], latarray[ylat]))
    
    #the the desired level u and v first time
    
    #get U
    global u
    u = ncfile.variables['u'][:,0,:,:] #[time, pressure level, u, v]
    global v
    v = ncfile.variables['v'][:,0,:,:] #[time, pressure level, u, v]
    
    traj_lon = []
    traj_lat = []
        
    global _Xmax 
    global _Xmin
    global _Ymax
    global _Ymin
    _Xmax = lonarray.size
    _Xmin = 0
    _Ymax = latarray.size
    _Ymin = 0    
    
    t = 0  # begining of time index (backward and forward)
    TM = 10 # end of time index
    Tm = 0
    signe = 1 # sine -1 backward trajectory
    
    # time gap between datasets
    _DT = 21600 #( 3600 * hours ), 3600sec = 1h 
    
    Ut = interp(u[t,:,:],lon,lat)
    Vt = interp(v[t,:,:],lon,lat)
    
    print 'first interp',t,Ut, Vt, lon, lat
    
    Utm1 = Ut
    Vtm1 = Vt
    
    #print 'plot:',t,lon,lat
    print 'plot:',t,to180_180(lon),lat
    
    traj_lon.append(to180_180(lon))
    traj_lat.append(lat)
    
    while (t < TM):
    #while (t > Tm):
        
        Ut=Utm1
        Vt=Vtm1
        
        #t=t+signe # set t+1 or t-1 forward backward
        
        LON0tm1=lon
        LAT0tm1=lat
        DR=1000000
        iteration=0
        
        while(DR > 0.1 and iteration < 15):
            Umig=signe*(Utm1+Ut)/2
            Vmig=signe*(Vtm1+Vt)/2
            V = math.sqrt( Umig * Umig + Vmig * Vmig )
            D=V*_DT
            alpha = math.atan2(Vmig,Umig) * _R2D
            LON1tm1, LAT1tm1 = mou(lon,lat,D,alpha)
            
            # check lon lat if its inside perimeter
            
            #print LON1tm1,LAT1tm1
            #DR = distance(LAT1tm1,LON1tm1,LAT0tm1,LON0tm1)
            DR=dist(LON1tm1,LAT1tm1,LON0tm1,LAT0tm1)
            Utm1=interp(u[t,:,:],LON1tm1,LAT1tm1)
            Vtm1=interp(v[t,:,:],LON1tm1,LAT1tm1)
            LON0tm1=LON1tm1
            LAT0tm1=LAT1tm1
            iteration=iteration+1
        if(iteration > 15):
            print 'lack of convergence: DR:=',DR,'m'
        #print 'plot:',t+signe,LON1tm1,LAT1tm1
        print 'plot:',t+signe,to180_180(LON1tm1),LAT1tm1
        #.plot([point[0], point2[0]], [point[1], point2[1]])        
        #linia(LONt,LATt,LON1tm1,LAT1tm1)
        #plt.plot(to180_180(lon),lat,to180_180(LON1tm1),LAT1tm1)
        traj_lon.append(to180_180(LON1tm1))
        traj_lat.append(LAT1tm1)
        
        lon=LON1tm1
        lat=LAT1tm1
        t=t+signe
        
    fi(0)
            
    xp, yp = m(traj_lon, traj_lat)
    #plt.plot(xp,yp,'ro-')
    plt.plot(xp,yp,'bo-')
    plt.show()
    
#--------------------------def onclick-----------------------------------------
#after the image is open, click and get the desired lat lng
def onclick(event):
    afterLatLon(event.xdata, event.ydata)

cid = fig.canvas.mpl_connect('button_press_event', onclick)

fig.show()



