# -*- coding: utf-8 -*-
"""
WRFout files to site time-series with tendencies

@author: Javier Sanz Rodrigo (jsrodrigo@cener.com)
"""
import numpy as np
from scipy.io import netcdf
from scipy import interpolate
from scipy import spatial
import matplotlib.dates as mdates
import utm
import sys
import glob
import netCDF4 
import datetime

# Constants
g = 9.81    # [m s-2]tar -zxvf
P0 = 100000. # Reference pressure [Pa]
T0 = 300.    # Reference temperature for perturbation temperature [K]
kappa = 0.2854  # Poisson constant (R/Cp)
R_air = 287.058  # Specific gas constant for dry air [J kg-1 K-1]
Cp_air = 1005   # Specific heat of air [J kg-1 K-1]
omega = 7.2921159e-5    # angular speed of the Earth [rad/s]

# Station coordinates
lat_s = 51.971   # degrees N
lon_s = 4.927    # degrees E
siteID = 'GABLS3'# identifier
fc = 2*omega*np.sin(lat_s*np.pi/180)
utmX_s, utmY_s, utm_zonenumber_s, utm_zoneletter_s = utm.from_latlon(lat_s,lon_s)

# Evaluation period
datefrom = datetime.datetime(2006,6,30,12,0,0)
dateto = datetime.datetime(2006,7,2,12,0,0)

# Column uniform grid in (microscale) cartesian coordinates
intmicro = 0   # = 1 if you want to interpolate veritcally to a microscale column (uniform grid) defined with Vres and Ztop
               # otherwise it will ouput the meso vertical levels

Vres = 20.0    # Vertical resolution [m]
Ztop = 4000.0  # Column height [m]

xmicro = utmX_s
ymicro = utmY_s
zmicro = np.linspace(0,Ztop,1+Ztop/Vres)
Nz = zmicro.shape[0]
Zmicro, Ymicro, Xmicro = np.meshgrid(zmicro, ymicro, xmicro, indexing='ij')
XYZmicro = zip(np.ravel(Zmicro), np.ravel(Ymicro), np.ravel(Xmicro))
XYmicro = zip(np.ravel(Ymicro[0,:,:]), np.ravel(Xmicro[0,:,:]))

# ----------------------------------------------------------------------
# Load simulation data
# ----------------------------------------------------------------------
# Simulation filename
dirdata = 'meso'	# folder with wrfout* files
simID = 'YSU-LES'	
dom = 2			# domain to extract time series from wrfout_d0[dom]*
pbl = 'YSU'
dxmeso = 3000.0		# horizontal resolution of the mesoscale grid
L = 3000.0           # Length of microscale domain Lx = Ly for spatial avaraging
                     # L = 0            linear interpolation to site coordinates
                     # 0 < L < dxmeso   nearest grid point
                     # L > dxmeso       spatial average over Nav points
Nav = int(L/dxmeso) + 1	# number of points to include in the spatial averaging box
simdes = ': WRF3.8, ERA-Interim, YSU until 1km then LES(tke) until 100-m res. Data from d0' + str(dom) + ' at ' + str(int(dxmeso/1000)) + ' km res, Lav = ' + str(int(L/1000)) + ' km'
simdes = ': WRF3.8, ERA-Interim, YSU 9km > 3km >1km. Data from d0' + str(dom) + ' at ' + str(int(dxmeso/1000)) + ' km res, Lav = ' + str(int(L/1000)) + ' km'

wrfoutfiles = sorted(glob.glob('./'+dirdata+'/wrfout_d0'+str(dom)+'*'))
#auxfiles = wrfoutfiles
#for i in range(0,len(wrfoutfiles)):
#    auxfiles[i] = auxfiles[i].replace("wrfout","auxhist7")
    
Nfiles = len(wrfoutfiles)
Nblock = 10	# use this to partition the output files in blocks. Nblock determines how many wrfout files will be post-processed per output file

def int3D(Z,ixnear,iynear,points3D,it):
    zdim = Z.shape[1]
    Z = np.reshape(Z[it,:,iynear[0]:iynear[-1]+1,ixnear[0]:ixnear[-1]+1],(zdim,3,3))
    f = interpolate.LinearNDInterpolator(points3D, np.ravel(Z))        
    return f

def int2D(Z,ixnear,iynear,points2D,it):
    Z = np.reshape(Z[it,iynear[0]:iynear[-1]+1,ixnear[0]:ixnear[-1]+1],(3,3))
    f = interpolate.LinearNDInterpolator(points2D, np.ravel(Z))        
    return f

cnt = 1
ifile = 0
while ifile < Nfiles:
    f2 = netCDF4.Dataset(wrfoutfiles[ifile])    
    print wrfoutfiles[ifile]
    sys.stdout.flush()
    
    zdim = len(f2.dimensions['bottom_top'])
    xdim = len(f2.dimensions['west_east'])
    ydim = len(f2.dimensions['south_north'])

    datetime_sim = netCDF4.num2date(f2.variables['XTIME'][:].squeeze(), f2.variables['XTIME'].units)    
    time_s = mdates.date2num(datetime_sim)   # number of days since 001-01-01 00:00:00 UTC, plus one
    tdim = len(time_s)

    idates = np.where(np.logical_and(time_s >= mdates.date2num(datefrom), 
                                 time_s <= mdates.date2num(dateto)))[0] 
    
    if idates.shape[0] > 0:    
        # Cell-centered coordinates 
        XLAT = f2.variables['XLAT'][:]     # WRF indexing XLAT[time,lat,lon]
        XLONG = f2.variables['XLONG'][:]
                
        # Nearest grid points to the site        
        points = zip(np.ravel(XLAT[0,:,:]),np.ravel(XLONG[0,:,:]))    
        tree = spatial.KDTree(points)
        dist, index = tree.query(np.array([lat_s,lon_s]),1) # nearest grid point
        inear = index % xdim     # index to XLONG
        jnear = index / xdim     # index to XLAT        
        ixnear = np.linspace(inear-1,inear+1,3) # three nearest points for interpolation
        iynear = np.linspace(jnear-1,jnear+1,3)
        if Nav == 1:        # nearest grid point
            ixav = inear
            iyav = jnear
        elif Nav == 2:      # four nearest grid points
            ixav = np.array([inear,inear+1])
            iyav = np.array([jnear,jnear+1])
        else:               
            if Nav % 2 == 1: # Nav (odd) nearest points 
                ixav = np.arange(inear - 0.5*(Nav-1),inear + 0.5*(Nav-1) + 1)                
                iyav = np.arange(jnear - 0.5*(Nav-1),jnear + 0.5*(Nav-1) + 1)
            else:
                ixav = np.arange(inear - 0.5*Nav +1,inear + 0.5*Nav + 1)                
                iyav = np.arange(jnear - 0.5*Nav +1,jnear + 0.5*Nav + 1)                

        # 3D fields. WRF indexing Z[tdim,ydim,xdim]
        U10 = f2.variables['U10'][:]    # 10-m U velocity component [m s-1]
        V10 = f2.variables['V10'][:]    # 10-m V velocity component [m s-1]
        T2 = f2.variables['T2'][:]      # 2-m temperature [K]
        TSK = f2.variables['TSK'][:]    # skin temperature [K]
        UST = f2.variables['UST'][:]    # M-O friction velocity [m s-1]
        Psfc = f2.variables['PSFC'][:]  # Surface Pressure [Pa]    
        HFX = f2.variables['HFX'][:]    # sensible heat flux [W m-2]
        LH = f2.variables['LH'][:]      # latent heat flux [W m-2]
       
        # 4D fields. WRF indexing Z[tdim,zdim,ydim,xdim]    
        U = f2.variables['U'][:]
        U = 0.5*(U[:,:,:,0:-1] + U[:,:,:,1:]) 
        U = np.concatenate((np.zeros((tdim,1,ydim,xdim)),U), axis = 1)
        
        V = f2.variables['V'][:]
        V = 0.5*(V[:,:,0:-1,:] + V[:,:,1:,:])
        V = np.concatenate((np.zeros((tdim,1,ydim,xdim)),V), axis = 1)

        W = f2.variables['W'][:]
        W = 0.5*(W[:,0:-1,:,:] + W[:,1:,:,:])
        W = np.concatenate((np.zeros((tdim,1,ydim,xdim)),W), axis = 1)
    
        Th = f2.variables['T'][:]+T0    
        Th = np.concatenate((np.reshape(TSK,(tdim,1,ydim,xdim)),Th), axis = 1)    

        PH = f2.variables['PH'][:]      # Base geopotential [Pa]
        PHB = f2.variables['PHB'][:]    # Perturbation geopotential [Pa]
        HGT = f2.variables['HGT'][:]      # terrain height [m]
                
        HGT2 = np.zeros((tdim,zdim+1,ydim,xdim))
        for iz in range(0,zdim+1):
            HGT2[:,iz,:,:] = HGT
            
        Zagl = (PH + PHB)/g - HGT2
        Zagl = 0.5*(Zagl[:,0:-1,:,:] + Zagl[:,1:,:,:])
        Zagl = np.concatenate((np.zeros((tdim,1,ydim,xdim)),Zagl), axis = 1)
            
        # Momentum budget components
        MUU = f2.variables['MUU'][:]    # mu-coupled u [Pa m s-1]
        MUU = 0.5*(MUU[:,:,0:-1] + MUU[:,:,1:]) 
        MUV = f2.variables['MUV'][:]    # mu-coupled v [Pa m s-1]
        MUV = 0.5*(MUV[:,0:-1,:] + MUV[:,0:-1,:]) 
        MUT = f2.variables['MUT'][:]    # mu-coupled T [Pa m s-1]
        
        RU_TEND = f2.variables['RU_TEND'][:]        # dU/dt
        RU_TEND = 0.5*(RU_TEND[:,:,:,0:-1] + RU_TEND[:,:,:,1:]) 
        RV_TEND = f2.variables['RV_TEND'][:]        
        RV_TEND = 0.5*(RV_TEND[:,:,0:-1,:] + RV_TEND[:,:,0:-1,:]) 
        
        RU_TEND_ADV = f2.variables['RU_TEND_ADV'][:]        # advection
        RU_TEND_ADV = 0.5*(RU_TEND_ADV[:,:,:,0:-1] + RU_TEND_ADV[:,:,:,1:]) 
        RV_TEND_ADV = f2.variables['RV_TEND_ADV'][:]        
        RV_TEND_ADV = 0.5*(RV_TEND_ADV[:,:,0:-1,:] + RV_TEND_ADV[:,:,0:-1,:]) 
        
        RU_TEND_PGF = f2.variables['RU_TEND_PGF'][:]        # pressure gradient
        RU_TEND_PGF = 0.5*(RU_TEND_PGF[:,:,:,0:-1] + RU_TEND_PGF[:,:,:,1:]) 
        RV_TEND_PGF = f2.variables['RV_TEND_PGF'][:]        
        RV_TEND_PGF = 0.5*(RV_TEND_PGF[:,:,0:-1,:] + RV_TEND_PGF[:,:,0:-1,:]) 
         
        RU_TEND_COR = f2.variables['RU_TEND_COR'][:]        # coriolis
        RU_TEND_COR = 0.5*(RU_TEND_COR[:,:,:,0:-1] + RU_TEND_COR[:,:,:,1:]) 
        RV_TEND_COR = f2.variables['RV_TEND_COR'][:]        
        RV_TEND_COR = 0.5*(RV_TEND_COR[:,:,0:-1,:] + RV_TEND_COR[:,:,0:-1,:]) 
        
        RU_TEND_PHYS = f2.variables['RU_TEND_PHYS'][:]        # PBL scheme (vertical diffusion)
        RU_TEND_PHYS = 0.5*(RU_TEND_PHYS[:,:,:,0:-1] + RU_TEND_PHYS[:,:,:,1:]) 
        RV_TEND_PHYS = f2.variables['RV_TEND_PHYS'][:]        
        RV_TEND_PHYS = 0.5*(RV_TEND_PHYS[:,:,0:-1,:] + RV_TEND_PHYS[:,:,0:-1,:]) 
    
        T_TEND_ADV1 = f2.variables['T_TEND_ADV'][:]        # advection (potential temperature)
                # for some reason T_TEND_ADV is not writable so I'll create a new T_TEND_ADV    
        T_TEND_ADV = np.zeros((tdim,zdim,ydim,xdim)) 
        
        for iz in range(0,zdim):    # divide by mu to obtain [m s-2]
            RU_TEND[:,iz,:,:]       = RU_TEND[:,iz,:,:]/MUU
            RV_TEND[:,iz,:,:]       = RV_TEND[:,iz,:,:]/MUV
            RU_TEND_ADV[:,iz,:,:]   = RU_TEND_ADV[:,iz,:,:]/MUU
            RV_TEND_ADV[:,iz,:,:]   = RV_TEND_ADV[:,iz,:,:]/MUV
            RU_TEND_PGF[:,iz,:,:]   = RU_TEND_PGF[:,iz,:,:]/MUU
            RV_TEND_PGF[:,iz,:,:]   = RV_TEND_PGF[:,iz,:,:]/MUV
            RU_TEND_COR[:,iz,:,:]   = RU_TEND_COR[:,iz,:,:]/MUU
            RV_TEND_COR[:,iz,:,:]   = RV_TEND_COR[:,iz,:,:]/MUV
            RU_TEND_PHYS[:,iz,:,:]  = RU_TEND_PHYS[:,iz,:,:]/MUU
            RV_TEND_PHYS[:,iz,:,:]  = RV_TEND_PHYS[:,iz,:,:]/MUV        
            T_TEND_ADV[:,iz,:,:]    = T_TEND_ADV1[:,iz,:,:]/MUT
    
        RU_TEND = np.concatenate((np.reshape(RU_TEND[:,0,:,:],(tdim,1,ydim,xdim)),RU_TEND), axis = 1)
        RV_TEND = np.concatenate((np.reshape(RV_TEND[:,0,:,:],(tdim,1,ydim,xdim)),RV_TEND), axis = 1)
        RU_TEND_ADV = np.concatenate((np.reshape(RU_TEND_ADV[:,0,:,:],(tdim,1,ydim,xdim)),RU_TEND_ADV), axis = 1)
        RV_TEND_ADV = np.concatenate((np.reshape(RV_TEND_ADV[:,0,:,:],(tdim,1,ydim,xdim)),RV_TEND_ADV), axis = 1)
        RU_TEND_PGF = np.concatenate((np.reshape(RU_TEND_PGF[:,0,:,:],(tdim,1,ydim,xdim)),RU_TEND_PGF), axis = 1)
        RV_TEND_PGF = np.concatenate((np.reshape(RV_TEND_PGF[:,0,:,:],(tdim,1,ydim,xdim)),RV_TEND_PGF), axis = 1)
        RU_TEND_COR = np.concatenate((np.reshape(RU_TEND_COR[:,0,:,:],(tdim,1,ydim,xdim)),RU_TEND_COR), axis = 1)
        RV_TEND_COR = np.concatenate((np.reshape(RV_TEND_COR[:,0,:,:],(tdim,1,ydim,xdim)),RV_TEND_COR), axis = 1)
        RU_TEND_PHYS = np.concatenate((np.reshape(RU_TEND_PHYS[:,0,:,:],(tdim,1,ydim,xdim)),RU_TEND_PHYS), axis = 1)
        RV_TEND_PHYS = np.concatenate((np.reshape(RV_TEND_PHYS[:,0,:,:],(tdim,1,ydim,xdim)),RV_TEND_PHYS), axis = 1)
        T_TEND_ADV = np.concatenate((np.reshape(T_TEND_ADV[:,0,:,:],(tdim,1,ydim,xdim)),T_TEND_ADV), axis = 1)
        
        # LHS terms
        Utend = (1/fc)*RU_TEND
        Vtend = (1/fc)*RV_TEND
        
        # RHS terms
        Vg = -(1/fc)*RU_TEND_PGF    
        Ug = (1/fc)*RV_TEND_PGF    
        Sg = (Ug**2 + Vg**2)**0.5
        WDg = 180 + np.arctan2(Ug,Vg)*180/np.pi
        
        Uadv = (1/fc)*RU_TEND_ADV
        Vadv = (1/fc)*RV_TEND_ADV
        
        Ucor = (1/fc)*RU_TEND_COR
        Vcor = (1/fc)*RV_TEND_COR
        
        Uphys = (1/fc)*RU_TEND_PHYS
        Vphys = (1/fc)*RV_TEND_PHYS
    
        # Advective RHS term potential temperature equation
        Thadv = T_TEND_ADV
        
        # Coordinates of WRF grid in the UTM-proyected catersian system
        XUTM = np.zeros((zdim,3,3))
        YUTM = np.zeros((zdim,3,3))
        ixnear = np.linspace(inear-1,inear+1,3)
        iynear = np.linspace(jnear-1,jnear+1,3)
        for ix in range(0,3):
            for iy in range(0,3):
                for iz in range(0,zdim):
                    XUTM[iz,iy,ix], YUTM[iz,iy,ix], utm_zonenumber, utm_zoneletter = utm.from_latlon(XLAT[0,iynear[iy],ixnear[ix]],XLONG[0,iynear[iy],ixnear[ix]], force_zone_number = utm_zonenumber_s)
        XUTM = np.vstack((np.reshape(XUTM[0,:,:],(1,3,3)),XUTM))
        YUTM = np.vstack((np.reshape(YUTM[0,:,:],(1,3,3)),YUTM))
        points2D = np.array((np.ravel(YUTM[0,:,:]),np.ravel(XUTM[0,:,:]))).T  
        
        if intmicro == 1: # vertical interpolation to microscale grid
            U_s = np.zeros((tdim,Nz));     V_s = np.zeros((tdim,Nz));     W_s = np.zeros((tdim,Nz));
            Th_s = np.zeros((tdim,Nz))
            Utend_s = np.zeros((tdim,Nz));     Vtend_s = np.zeros((tdim,Nz));     
            Ug_s = np.zeros((tdim,Nz));     Vg_s = np.zeros((tdim,Nz));         
            Uadv_s = np.zeros((tdim,Nz));     Vadv_s = np.zeros((tdim,Nz))    
            Ucor_s = np.zeros((tdim,Nz));     Vcor_s = np.zeros((tdim,Nz))    
            Uphys_s = np.zeros((tdim,Nz));     Vphys_s = np.zeros((tdim,Nz))    
            Thadv_s = np.zeros((tdim,Nz));

            ust_s = np.zeros((tdim)); T2_s  = np.zeros((tdim)); TSK_s = np.zeros((tdim)); 
            HFX_s = np.zeros((tdim)); LH_s = np.zeros((tdim)); Psfc_s = np.zeros((tdim));
            if L == 0: # interpolate to site coordinates
                for it in range(0,tdim):
                    ZUTM = Zagl[it,:,iynear[0]:iynear[-1]+1,ixnear[0]:ixnear[-1]+1]
                    points3D = np.array((np.ravel(ZUTM), np.ravel(YUTM),np.ravel(XUTM))).T  
                    
                    U_s[it,:] = int3D(U,ixnear,iynear,points3D,it)(XYZmicro)
                    V_s[it,:] = int3D(V,ixnear,iynear,points3D,it)(XYZmicro)
                    W_s[it,:] = int3D(W,ixnear,iynear,points3D,it)(XYZmicro)
                    Th_s[it,:] = int3D(Th,ixnear,iynear,points3D,it)(XYZmicro)
                    Utend_s[it,:] = int3D(Utend,ixnear,iynear,points3D,it)(XYZmicro)
                    Vtend_s[it,:] = int3D(Vtend,ixnear,iynear,points3D,it)(XYZmicro)
                    Ug_s[it,:] = int3D(Ug,ixnear,iynear,points3D,it)(XYZmicro)
                    Vg_s[it,:] = int3D(Vg,ixnear,iynear,points3D,it)(XYZmicro)
                    Uadv_s[it,:] = int3D(Uadv,ixnear,iynear,points3D,it)(XYZmicro)
                    Vadv_s[it,:] = int3D(Vadv,ixnear,iynear,points3D,it)(XYZmicro)
                    Ucor_s[it,:] = int3D(Ucor,ixnear,iynear,points3D,it)(XYZmicro)
                    Vcor_s[it,:] = int3D(Vcor,ixnear,iynear,points3D,it)(XYZmicro)
                    Uphys_s[it,:] = int3D(Uphys,ixnear,iynear,points3D,it)(XYZmicro)
                    Vphys_s[it,:] = int3D(Vphys,ixnear,iynear,points3D,it)(XYZmicro)
                    Thadv_s[it,:] = int3D(Thadv,ixnear,iynear,points3D,it)(XYZmicro)
                    ust_s[it] = int2D(UST,ixnear,iynear,points2D,it)(XYmicro)
                    T2_s[it] = int2D(T2,ixnear,iynear,points2D,it)(XYmicro)
                    TSK_s[it] = int2D(TSK,ixnear,iynear,points2D,it)(XYmicro)
                    HFX_s[it] = int2D(HFX,ixnear,iynear,points2D,it)(XYmicro)
                    LH_s[it] = int2D(LH,ixnear,iynear,points2D,it)(XYmicro)
                    Psfc_s[it] = int2D(Psfc,ixnear,iynear,points2D,it)(XYmicro)
            else: # spatial averaging based on L
                for it in range(0,tdim):
                    Znear = Zagl[it,:,jnear,inear]
                    U_s[it,:] = np.interp(zmicro,Znear,np.mean(U[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    V_s[it,:] = np.interp(zmicro,Znear,np.mean(V[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    W_s[it,:] = np.interp(zmicro,Znear,np.mean(W[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Th_s[it,:] = np.interp(zmicro,Znear,np.mean(Th[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Utend_s[it,:] = np.interp(zmicro,Znear,np.mean(Utend[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Vtend_s[it,:] = np.interp(zmicro,Znear,np.mean(Vtend[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Ug_s[it,:] = np.interp(zmicro,Znear,np.mean(Ug[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Vg_s[it,:] = np.interp(zmicro,Znear,np.mean(Vg[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Uadv_s[it,:] = np.interp(zmicro,Znear,np.mean(Uadv[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Vadv_s[it,:] = np.interp(zmicro,Znear,np.mean(Vadv[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Ucor_s[it,:] = np.interp(zmicro,Znear,np.mean(Ucor[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Vcor_s[it,:] = np.interp(zmicro,Znear,np.mean(Vcor[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Uphys_s[it,:] = np.interp(zmicro,Znear,np.mean(Uphys[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Vphys_s[it,:] = np.interp(zmicro,Znear,np.mean(Vphys[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    Thadv_s[it,:] = np.interp(zmicro,Znear,np.mean(Thadv[it,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2)))
                    ust_s[it] = np.mean(UST[it,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (0,1))
                    T2_s[it] = np.mean(T2[it,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (0,1))
                    TSK_s[it] = np.mean(TSK[it,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (0,1))
                    HFX_s[it] = np.mean(HFX[it,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (0,1))
                    LH_s[it] = np.mean(LH[it,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (0,1))
                    Psfc_s[it] = np.mean(Psfc[it,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (0,1))
    
        else:   # no vertical interpolation
            if L == 0:  # nearest grid point
                z_s = Zagl[:,:,jnear,inear]
                U_s = U[:,:,jnear,inear]
                V_s = V[:,:,jnear,inear]
                W_s = W[:,:,jnear,inear]
                Th_s = Th[:,:,jnear,inear]
                Utend_s = Utend[:,:,jnear,inear]
                Vtend_s = Vtend[:,:,jnear,inear]
                Ug_s = Ug[:,:,jnear,inear]
                Vg_s = Vg[:,:,jnear,inear]
                Uadv_s = Uadv[:,:,jnear,inear]
                Vadv_s = Vadv[:,:,jnear,inear]
                Ucor_s = Ucor[:,:,jnear,inear]
                Vcor_s = Vcor[:,:,jnear,inear]
                Uphys_s = Uphys[:,:,jnear,inear]
                Vphys_s = Vphys[:,:,jnear,inear]
                Thadv_s = Thadv[:,:,jnear,inear]                
                ust_s = UST[:,jnear,inear]
                T2_s = T2[:,jnear,inear]
                TSK_s = TSK[:,jnear,inear]
                HFX_s = HFX[:,jnear,inear]
                LH_s = LH[:,jnear,inear]
                Psfc_s = Psfc[:,jnear,inear]
            else:   
                z_s = np.mean(Zagl[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                U_s = np.mean(U[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                V_s = np.mean(V[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                W_s = np.mean(W[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Th_s = np.mean(Th[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Utend_s = np.mean(Utend[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Vtend_s = np.mean(Vtend[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Ug_s = np.mean(Ug[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Vg_s = np.mean(Vg[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Uadv_s = np.mean(Uadv[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Vadv_s = np.mean(Vadv[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Ucor_s = np.mean(Ucor[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Vcor_s = np.mean(Vcor[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Uphys_s = np.mean(Uphys[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Vphys_s = np.mean(Vphys[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))
                Thadv_s = np.mean(Thadv[:,:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (2,3))          
                ust_s = np.mean(UST[:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2))          
                T2_s = np.mean(T2[:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2))
                TSK_s = np.mean(TSK[:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2))
                HFX_s = np.mean(HFX[:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2))
                LH_s = np.mean(LH[:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2))
                Psfc_s = np.mean(Psfc[:,iyav[0]:iyav[-1]+1,ixav[0]:ixav[-1]+1], axis = (1,2))
            z_s = np.mean(z_s, axis = 0)    
    
        if cnt == 1:
            U_out = U_s; V_out = V_s; W_out = W_s; Th_out = Th_s
            Ug_out = Ug_s; Vg_out = Vg_s; Uadv_out = Uadv_s; Vadv_out = Vadv_s;
            Utend_out = Utend_s; Vtend_out = Vtend_s; Thadv_out = Thadv_s;
            Ucor_out = Ucor_s; Vcor_out = Vcor_s; 
            Uphys_out = Uphys_s; Vphys_out = Vphys_s
            ust_out = ust_s; T2_out = T2_s; TSK_out = TSK_s; HFX_out = HFX_s; LH_out = LH_s; Psfc_out = Psfc_s
            time_out = time_s; 
        else:
            time_out = np.hstack((time_out,time_s))
            U_out = np.vstack((U_out,U_s))
            V_out = np.vstack((V_out,V_s))
            W_out = np.vstack((W_out,W_s))
            Th_out = np.vstack((Th_out,Th_s))
            Ug_out = np.vstack((Ug_out,Ug_s))
            Vg_out = np.vstack((Vg_out,Vg_s))
            Utend_out = np.vstack((Utend_out,Utend_s))
            Vtend_out = np.vstack((Vtend_out,Vtend_s))
            Uadv_out = np.vstack((Uadv_out,Uadv_s))
            Vadv_out = np.vstack((Vadv_out,Vadv_s))
            Ucor_out = np.vstack((Ucor_out,Ucor_s))
            Vcor_out = np.vstack((Vcor_out,Vcor_s))
            Uphys_out = np.vstack((Uphys_out,Uphys_s))
            Vphys_out = np.vstack((Vphys_out,Vphys_s))
            Thadv_out = np.vstack((Thadv_out,Thadv_s))
            ust_out = np.hstack((ust_out,ust_s))
            T2_out = np.hstack((T2_out,T2_s))
            TSK_out = np.hstack((TSK_out,TSK_s))
            HFX_out = np.hstack((HFX_out,HFX_s))
            LH_out = np.hstack((LH_out,LH_s))
            Psfc_out = np.hstack((Psfc_out,Psfc_s))
        
        if (cnt == Nblock) or (ifile == Nfiles-1):
            isort = np.argsort(time_out)      
            time_out = time_out[isort]
            if intmicro == 1:
                z_out = zmicro
            else:
                z_out = z_s                
            U_out = U_out[isort,:]
            V_out = V_out[isort,:]
            W_out = W_out[isort,:]
            Th_out = Th_out[isort,:]
            Utend_out = Utend_out[isort,:]
            Vtend_out = Vtend_out[isort,:]
            Ug_out = Ug_out[isort,:]
            Vg_out = Vg_out[isort,:]
            Uadv_out = Uadv_out[isort,:]
            Vadv_out = Vadv_out[isort,:]
            Ucor_out = Ucor_out[isort,:]
            Vcor_out = Vcor_out[isort,:]
            Uphys_out = Uphys_out[isort,:]
            Vphys_out = Vphys_out[isort,:]
            Thadv_out = Thadv_out[isort,:]
            ust_out = ust_out[isort]
            T2_out = T2_out[isort]
            TSK_out = TSK_out[isort]
            HFX_out = HFX_out[isort]
            LH_out = LH_out[isort]
            Psfc_out = Psfc_out[isort]
                
            #datefrom = mdates.num2date(time_out[0]).strftime("%Y-%m-%d")
            #dateto = mdates.num2date(time_out[-1]).strftime("%Y-%m-%d")
            
            # Write to netcdf file
            fileout= siteID + '_' + simID + '_tend_L' + str(int(L)) + '_' + str(ifile) + '.nc'
            f = netcdf.netcdf_file(fileout, 'w')
            f.history = siteID +': '+ simdes
            
            f.createDimension('time', np.shape(time_out)[0])
            f.createDimension('site', 1)
            f.createDimension('z', len(z_out))
            
            lats = f.createVariable('lat', 'float', ('site',))
            lats.name = 'Site latitude'
            lats.units = 'degrees North'
            lats[:] = lat_s
            
            lons = f.createVariable('lon', 'float', ('site',))
            lons.name = 'Site longitude'
            lons.units = 'degrees East'
            lons[:] = lon_s
            
            fcs = f.createVariable('fc', 'float', ('site',))
            fcs.name = 'Coriolis parameter'
            fcs.units = 's-1'
            fcs[:] = fc
            
            times = f.createVariable('time', 'float', ('time',))
            times.name = 'Time'
            times.units = 'Days since 001-01-01 00:00:00 UTC, plus one'
            times[:] = time_out
            
            heights = f.createVariable('z', 'float', ('z',))
            heights.name = 'Height above ground level'
            heights.units = 'm'
            heights[:] = z_out
    
            Us = f.createVariable('U', 'float', ('time','z',))
            Us.name = 'U velocity component' 
            Us.units = 'm s-1'
            Us[:] = U_out
            
            Vs = f.createVariable('V', 'float', ('time','z',))
            Vs.name = 'V velocity component' 
            Vs.units = 'm s-1'
            Vs[:] = V_out
            
            Ws = f.createVariable('W', 'float', ('time','z',))
            Ws.name = 'W velocity component'
            Ws.units = 'm s-1'
            Ws[:] = W_out

            Ths = f.createVariable('Th', 'float', ('time','z',))
            Ths.name = 'Potential temperature' 
            Ths.units = 'K'
            Ths[:] = Th_out
            
            Utends = f.createVariable('Utend', 'float', ('time','z',))
            Utends.name = 'Tendency momentum LHS term (divided by fc) U-component' 
            Utends.units = 'm s-1'
            Utends[:] = Utend_out
            
            Vtends = f.createVariable('Vtend', 'float', ('time','z',))
            Vtends.name = 'Tendency momentum LHS term (divided by fc) V-component' 
            Vtends.units = 'm s-1'
            Vtends[:] = Vtend_out
            
            Ugs = f.createVariable('Ug', 'float', ('time','z',))
            Ugs.name = 'Geostrophic wind U-component = Pressure gradient RHS term (divided by fc) V-component' 
            Ugs.units = 'm s-1'
            Ugs[:] = Ug_out
            
            Vgs = f.createVariable('Vg', 'float', ('time','z',))
            Vgs.name = 'Geostrophic wind V-component = - Pressure gradient RHS term (divided by fc) U-component' 
            Vgs.units = 'm s-1'
            Vgs[:] = Vg_out
            
            Uadvs = f.createVariable('Uadv', 'float', ('time','z',))
            Uadvs.name = 'Advective momentum RHS term (divided by fc) U-component' 
            Uadvs.units = 'm s-1'
            Uadvs[:] = Uadv_out
            
            Vadvs = f.createVariable('Vadv', 'float', ('time','z',))
            Vadvs.name = 'Advective momentum RHS term (divided by fc) V-component' 
            Vadvs.units = 'm s-1'
            Vadvs[:] = Vadv_out
            
            Ucors = f.createVariable('Ucor', 'float', ('time','z',))
            Ucors.name = 'Coriolis momentum RHS term (divided by fc) U-component = V velocity component' 
            Ucors.units = 'm s-1'
            Ucors[:] = Ucor_out
            
            Vcors = f.createVariable('Vcor', 'float', ('time','z',))
            Vcors.name = 'Coriolis momentum RHS term (divided by fc) V-component = -U velocity component' 
            Vcors.units = 'm s-1'
            Vcors[:] = Vcor_out
            
            Uphyss = f.createVariable('Uphys', 'float', ('time','z',))
            Uphyss.name = 'Vertical diffusion momentum RHS term (divided by fc) U-component' 
            Uphyss.units = 'm s-1'
            Uphyss[:] = Uphys_out
            
            Vphyss = f.createVariable('Vphys', 'float', ('time','z',))
            Vphyss.name = 'Vertical diffusion momentum RHS term (divided by fc) V-component' 
            Vphyss.units = 'm s-1'
            Vphyss[:] = Vphys_out
            
            Thadvs = f.createVariable('Thadv', 'float', ('time','z',))
            Thadvs.name = 'Potential temperature advective RHS term' 
            Thadvs.units = 'K s-1'
            Thadvs[:] = Thadv_out

            usts = f.createVariable('ust', 'float', ('time',))
            usts.name = 'Monin Obukhov friction velocity' 
            usts.units = 'm s-1'
            usts[:] = ust_out
            
            T2s = f.createVariable('T2', 'float', ('time',))
            T2s.name = '2-m temperature' 
            T2s.units = 'K'
            T2s[:] = T2_out
            
            TSKs = f.createVariable('TSK', 'float', ('time',))
            TSKs.name = 'skin temperature' 
            TSKs.units = 'K'
            TSKs[:] = TSK_out
            
            HFXs = f.createVariable('HFX', 'float', ('time',))
            HFXs.name = 'Upward sensible heat flux at surface' 
            HFXs.units = 'W m-2'
            HFXs[:] = HFX_out
            
            LHs = f.createVariable('LH', 'float', ('time',))
            LHs.name = 'Upward latent heat flux at surface' 
            LHs.units = 'W m-2'
            LHs[:] = LH_out
            
            Psfcs = f.createVariable('Psfc', 'float', ('time',))
            Psfcs.name = 'Surface pressure' 
            Psfcs.units = 'Pa'
            Psfcs[:] = Psfc_out

            f.close()
            cnt = 1        
            ifile = ifile + 1
        else:
            cnt = cnt + 1
            ifile = ifile + 1
    else:
        ifile = ifile +1
    
    #Save to text files
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_Ug.txt', Ug_out, delimiter=',', fmt='%1.3f')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_Vg.txt', Vg_out, delimiter=',', fmt='%1.3f')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_W.txt', W_out, delimiter=',', fmt='%1.3f')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_Thadv.txt', Thadv_out, delimiter=',', fmt='%1.2e')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_Uadv.txt', Uadv_out, delimiter=',', fmt='%1.2e')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_Vadv.txt', Vadv_out, delimiter=',', fmt='%1.2e')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_t.txt', time_out, delimiter=',', fmt='%1.2e')
    #np.savetxt(siteID+'_'+datefrom+'_'+dateto+'_z.txt', z_out, delimiter=',', fmt='%1.2e')
