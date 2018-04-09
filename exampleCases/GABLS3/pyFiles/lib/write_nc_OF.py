# -*- coding: utf-8 -*-
"""
Auxiliary function to create netcdf data. 
So far only works with z and t-dependent information

Created  Feb 08 2018.  Roberto Chavez

"""
import sys
import netCDF4 
from os.path import isfile

def OF2nc_headers(fileout, inputs, siteData, timeData, zData):
	n = len(timeData)
	m = len(zData)
	if not isfile(fileout):
		print 'Creating netcdf file!!'
		f = netCDF4.Dataset(fileout,'w',format= "NETCDF4")
		# create dimensions, variables and attributes:
		#f.history = inputs[0]['siteID'] + ': Meso forcing = ' + filemeso + ','.join('{}{}'.format(key, val) for key, val in model[0].items())	
		f.history = 'openfoam simulation'	
		f.createDimension('time', n)
		f.createDimension('z', m)
		f.createDimension('site', 1)
		lats = f.createVariable('lat', 'float', ('site',))
		lats.long_name = 'Site latitude'
		lats.units = 'degrees North'
		lats[:] = inputs[0]['lat']
		    
		lons = f.createVariable('lon', 'float', ('site',))
		lons.long_name = 'Site longitude'
		lons.units = 'degrees East'
		lons[:] = inputs[0]['lat']
				
		fcs = f.createVariable('fc', 'float', ('site',))
		fcs.long_name = 'Coriolis parameter'
		fcs.units = 's-1'
		fcs[:] = siteData['fc']
		   
		times = f.createVariable('time', 'float', ('time',))
		times.long_name = 'Time'
		times.units = 'Days since 001-01-01 00:00:00 UTC, plus one'
		times[:] = timeData
		    
		heights = f.createVariable('z', 'float', ('z',))
		heights.long_name = 'Height above ground level'
		heights.units = 'm'
		heights[:] = zData
		
		f.close()	
	
	else:
		print "File already exist"
		
		
def OF2nc_data(fileout, nDimensions, varName, varLongName, varUnits, varData):

	print "appending variable "+ varName 
	f = netCDF4.Dataset(fileout,'a',format= "NETCDF4")
	
	if (nDimensions==1):
		var = f.createVariable(varName, 'float', ('time',))
	elif (nDimensions==2):
		var = f.createVariable(varName, 'float', ('time','z',))
	else:
		print "only 1(time) or 2 (time,z) dimensions allowed so far"
		
	var.long_name = varLongName
	var.units = varUnits
	var[:] = varData
		
	f.close()

	print 'Saving ' + fileout
	sys.stdout.flush()
		