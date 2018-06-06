#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Main script to read openfoam set lines and convert them to both pickl and netcdf outputs
'''
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.dates as mdates
import os
sys.path.append(os.environ.get("PYTHON_COMMON_LIBRARIES"))
from libOFpostProcessing import readSetFile
from write_nc_OF import OF2nc_headers, OF2nc_data

def main():
	#%%
#  INPUTS
	caseName = os.environ.get("OFCASE","OF case does not exist")
	OFpath = os.environ.get("OFPATH","Path do not exist")
		
	outputFile=OFpath+'/'+caseName
	
	# so far the names of the set files are hardcoded to the followin:g	
	scalarsFileName='lineZ_alphaB_C1ast_C3_k_epsilon_nut_lm_lMY_p_rgh_p_Prt_T_Rig_RiG_gradTz_wTz_uStar_L.xy'
	UFileName='lineZ_U.xy'

	init_Inputs=pd.read_pickle(OFpath+'/initialization_inputs.pckl')
	
	inputs = init_Inputs['inputs']
	datefrom = init_Inputs['datefrom']

	omega   = 7.2921159e-5    # angular speed of the Earth [rad/s]
	fc  = 2*omega*np.sin(inputs[0]['lat']*np.pi/180)     # Coriolis parameter [s-1]
	
	# SAVE AS PICKL FORMAT
	varList = ['alphaB','C1ast','C3','k','epsilon','nut','lm','lMY',
				'p_rgh','p','Prt','T','RiG','gradTz','wTz','uStar','L']		
		
	U, t, coords = readSetFile(OFpath,UFileName,'U',True)
	
	fdataPy = dict(t=np.float64(t)/(3600*24.0) + mdates.date2num(datefrom), 
							z=coords['z'], ux=U['ux'], uy=U['uy'], uz=U['uz'])
		
	fdataPy['S'] = np.sqrt(U['ux']**2 + U['uy']**2 + U['uz']**2) 
	fdataPy['D'] = 180.0 + np.arctan2(U['ux'],U['uy'])*180.0/np.pi
																
	for ii in varList:
		if (ii is not 'U'):
			varData = readSetFile(OFpath,scalarsFileName,ii)
			fdataPy[ii] = varData['data']

			
	with open(outputFile+'.pckl', 'wb') as f:
		pickle.dump(fdataPy, f, protocol=-1)
			
	# SAVE IN NETCDF FORMAT
	#inputname in OF, outputname in ncfile, long_name, units, interpolate to 1D, heigh above ground (in case of interpolation)
	vars2extract =[
			['T',    'Th',   'Potential Temperature' ,   'K',      0, 0],
			['T',    'T2',   '2m Potential Temperature', 'K',      1, 2.0],
			['k',    'TKE',  'Turbulent kinetic energy', 'm2 s-2', 0, 0],
			['nut',  'NU_T', 'Eddy Viscosity',           'm2 s-1', 0, 0],
			['uStar','ust',  'Friction velocity',        'm s-1' , 0, 0],
			['wTz',  'wt',   'Kinematic heat flux',      'K m s-1',0, 0],
			['lm',   'lm',   'Mixing length'  ,          'm',      0, 0],
			['L',    'L',    'Obukhov length'  ,         'm',      0, 0]
			]
		
	# read common data from the first variable defined
	U, t, coords = readSetFile(OFpath,UFileName,'U',True)
	siteData = {'fc': fc}
	zData = coords['z'].as_matrix()
	timeData = t/(3600*24.0) + mdates.date2num(datefrom)
		
	OF2nc_headers(outputFile+'.nc', inputs, siteData, timeData, zData)
		
	OF2nc_data(outputFile+'.nc',2 ,'U','U velocity component','m s-1', U['ux'].as_matrix().T)
	OF2nc_data(outputFile+'.nc',2 ,'V','V velocity component','m s-1', U['uy'].as_matrix().T)
			
	for ii in range(np.shape(vars2extract)[0]):			
		
		varNameInSetsFile=vars2extract[ii][0]
		varName     = vars2extract[ii][1]
		varLongName = vars2extract[ii][2]
		varUnits    = vars2extract[ii][3]
		print("proccesing variable: "+ varName)
	
		varData = readSetFile(OFpath,scalarsFileName,varNameInSetsFile)
		if (vars2extract[ii][4]==0):
			OF2nc_data(outputFile+'.nc',2 ,varName,varLongName,varUnits, varData['data'].as_matrix().T)
		else:
			zI = vars2extract[ii][5]   # height to extract
			index = (np.abs(zData-zI)).argmin()
			OF2nc_data(outputFile+'.nc',1 ,varName,varLongName,varUnits,varData['data'].as_matrix().T[:,index])
	
	
#%%	

	return 0
#%%
if __name__ == '__main__':
	main()

