#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script to write forcing, initial and boundary conditions of the CFD based
on a procesed mesoscale (WRF) output file.  
National Renewable Energy Centre of Spain (CENER) 

Created:   Feb/2018 Roberto Chavez Arroyo <rchavez@cener.com>

"""
import sys
import numpy as np
import os
import csv
import matplotlib.pyplot as plt
import datetime
import pickle

sys.path.append(os.environ.get("PYTHON_COMMON_LIBRARIES"))

from libOFpreProcessing import read_nc_wrf, plot_nc_tzData

#%%
def main():
	#%%  ######################## INPUT DATA ###################################
	openFoamCasePath= str(os.environ.get("OFPATH","OpenFOAM path does not exist"))
	tendenciesFile  = str(os.environ.get("tendenciesFile","OpenFOAM path does not exist"))
	datefromSys     = str(os.environ.get("datefrom"," "))
	datetoSys       = str(os.environ.get("dateto"," "))
	plotData        = float(os.environ.get("plotncData",0))
	#%%
	
	# Initialization/forcing period
	#datefrom = datetime.datetime(2006,7,1,0,0,0)   # 1/2 day spin up
	#datefrom = datetime.datetime(2006,7,1,12,0,0)   # 0 day spin up 
	#datefrom = datetime.datetime(2006,7,1,8,0,0)   # 4 hrs spin up 
	#dateto   = datetime.datetime(2006,7,2,12,0,0)
	
	datefrom = datetime.datetime.strptime(datefromSys,'%Y-%m-%d %H:%M')
	dateto   = datetime.datetime.strptime(datetoSys,'%Y-%m-%d %H:%M')	

	print(("Initializing "+openFoamCasePath))
	print(("\nForcing File: "+tendenciesFile))
	print(("\nPeriod to consider: from "+str(datefrom)+" until "+str(dateto) +"\n"))
	
	# List with: [ writeoutput data?, where (path) to write it?]
	writeInitialValues = [True, openFoamCasePath+'/forcing/']
	writeForcings      = [True, openFoamCasePath+'/forcing/']
	writeSurfaceValues = [True, openFoamCasePath+'/forcing/']
	
	# Mesoscale tendencies averaging settings
	mesoDataSetup = {
			'ncfile': tendenciesFile,  # tendencies file
			'datesOffset': 1.0,     # time offset in ncfile to get timestamp in "days since 0001-01-01 00:00:00 UTC".
			'datesSlope':  1.0,     # time slope in ncfile to get timestamp in "days since 0001-01-01 00:00:00 UTC".
			'tav':         60.0,    # Time averaging time used in simulations [min]
			'Lav':         9000.0,  # Spatial averaging [m]
			'ts':          10       # sampling frequency to evaluate [min]     
				}                                     
	
	# Constants
	P0      = 100000         # standard reference pressure for computing potential temperature [Pa]
	g       = 9.81           # [m s-2]
	R_air   = 287.058        # Specific gas constant for dry air [J kg-1 K-1]
	Cp_air  = 1005           # Specific heat of air [J kg-1 K-1]
	K       = 0.4            # von Karman constant
	kappa   =R_air/Cp_air    # Poisson constant (R/Cp)
			
	#%%     CODE###############################################

	# Extract data from mesoscale. The structure is [name in nc file, should they by multiplied by fc?]
	vars2extract =[ ['U',    False],
				['V',    False],
				['W',    False],
				['Th',   False],
				['Ug',   True],
				['Vg',   True],
				['Uadv', True],
				['Vadv', True],
				['Thadv',False],
				['T2',   False],
				['TSK',  False],
				['Psfc', False]]
	
	mesoData = read_nc_wrf(mesoDataSetup['ncfile'], mesoDataSetup['datesOffset'], mesoDataSetup['datesSlope'],
													vars2extract, datefrom, dateto, plotData)

	#  define if the values are provided by WRF outputs
	z    = mesoData['z']
	t    = mesoData['hrs_since_t0'] * 3600
	U    = mesoData['U']
	V    = mesoData['V']
	Th   = mesoData['Th']

	# set-up forcing: Momentum and energy budget from WRF used as driving force of the micro
	# Note that, the Gestrophic wind refers to the Pressure gradient componet 
	# which by convention is derived as dP/dx=-Vg  & dP/dy=Ug

	Uf = mesoData['Uadv'] - mesoData['Vg']    # mometum forcing east-west component.  
	Vf = mesoData['Vadv'] + mesoData['Ug']    # mometum forcing north-south component
	Wf = np.zeros(np.shape(Uf))               # mometum forcing vertical component

	Thf= mesoData['Thadv']                    # heat flux tendency
	

	T0  = mesoData['TSK']	
	T2  = mesoData['T2']
	Th0 = T0 * ( (P0/mesoData['Psfc'])**kappa )      # wall potential temperature [K]
	Th2 = T2 * ( (P0/mesoData['Psfc'])**kappa )     # 2-m potential temperature [K]	
	
#	us0  = mesoData['ust']
	rho0 = mesoData['Psfc']/(R_air*mesoData['T2'])
#	wTh0 = mesoData['wt']                 # wall heat flux [K m s-1]
		
#	beta0 = 1.0/Th0                                  # expansion coefficient [K-1]
#	L0 = -us0**3/(K*g*beta0*wTh0)   # Surface-layer Obukhov length [m]
	
	nt   = len(t)
	nz   = len(z)
	
	#%%  dump output files
	def writeForcing1D(fid,headerName,variable,n):
		fid.write(headerName+'\n (\n')    
		for j in range(n):
			fid.write('    ' + str(variable[j]) + '\n')
		fid.write(');\n\n')
	
	def writeForcing2D(fid,headerName,t,variable,nz,nt):
		fid.write(headerName + '\n (\n')   
		for n in range(nt):
			fid.write('    (' + str(t[n]) + ' ')
			for j in range(nz):
				fid.write(str(variable[n,j]) + ' ')
			fid.write(')\n')
		fid.write(');\n\n')
		
		
	if (writeForcings[0]):
		print ("writing forcing")
		
		if plotData:
			plt.figure(1); plot_nc_tzData('Uforcing', mesoData['hrs_since_t0'], mesoData['z'], Uf, datefrom,dateto, True)
			plt.figure(2); plot_nc_tzData('Vforcing', mesoData['hrs_since_t0'], mesoData['z'], Vf, datefrom,dateto, True)
			plt.figure(3); plot_nc_tzData('Thforcing',mesoData['hrs_since_t0'], mesoData['z'], Thf,datefrom,dateto, True)	
		
		fid = open(writeForcings[1]+'forcingTable','w') 
		writeForcing1D(fid,'sourceHeightsMomentum', z, nz) # Write the height list for the momentum forcing
		writeForcing2D(fid,'sourceTableMomentumX', t, Uf, nz, nt) # Write the x-dir momentum forcing
		writeForcing2D(fid,'sourceTableMomentumY', t, Vf, nz, nt) # Write the y-dir momentum forcing
		writeForcing2D(fid,'sourceTableMomentumZ', t, Wf, nz, nt) # Write the z-dir momentum forcing
		writeForcing1D(fid,'sourceHeightsTemperature', z, nz) # Write the height list for the temperature forcing
		writeForcing2D(fid,'sourceTableTemperature', t, Thf, nz, nt) # Write the temperature forcing
		fid.close()
		
	if (writeSurfaceValues[0]):
		print ("writing surface values")
		# Skin real temperature   
		with open(writeSurfaceValues[1]+'surfaceSkinTemperatureTable', 'w') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerows( np.column_stack((t,T0)) )
		# 2-m real temperature
		with open(writeSurfaceValues[1]+'surface2mTemperatureTable', 'w') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerows( np.column_stack((t,T2)) )
		# Skin potential temperature       
		with open(writeSurfaceValues[1]+'surfaceSkinPotentialTemperatureTable', 'w') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerows( np.column_stack((t,Th0)) )
		# 2-m potential temperature
		with open(writeSurfaceValues[1]+'surface2mPotentialTemperatureTable', 'w') as csvfile:
			writer = csv.writer(csvfile)
			writer.writerows( np.column_stack((t,Th2)) )
				
	
	if (writeInitialValues[0]):
		print ("writing initial values")
		fid = open(writeInitialValues[1]+'initialValues','w')   
		for j in range(nz):
			fid.write(' (' + str(z[j]) + ' ' + str(U[0,j]) + ' ' + str(V[0,j]) + ' ' + str(Th[0,j]) + ')\n')
		fid.close()
				
	# save forcing data for possibel further usage
	with open(openFoamCasePath+'/forcing.pckl', 'wb') as f:
		pickle.dump(mesoData, f, protocol=-1)
			
	# Save inputs of the case initialization for further usage in the postprocessing stage
	inputs = {'lat':mesoData['lat'][0], 'lon': mesoData['lon'][0], 'fc':mesoData['fc'][0]}
	outFilePy = {'datefrom':datefrom, 'dateto':dateto, 'inputs': inputs}
	with open(openFoamCasePath+'/initialization_inputs.pckl', 'wb') as f:
		pickle.dump(outFilePy, f, protocol=-1)


	return 0
#%%
if __name__ == '__main__':
	main()
