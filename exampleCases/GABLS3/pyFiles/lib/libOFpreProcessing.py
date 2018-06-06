#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Auxiliary functions for the preprocessing of the NEWAFoam
Feb2018 Roberto Chavez rchavez@cener.com
'''
import sys
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import datetime
import matplotlib.dates as mdates

def read_nc_wrf(ncfile, ncDatesOffset, ncDatesSlope, vars2extract, datefrom, dateto, plotData=False):
	# Extract data from mesoscale. 
	# vars2extract list should contain  [name of variable in nc file, multiply fc?]
     # so far only time, z-varying data is accepted in the netcfd file
	#       time requires an offset and slope to get ncTime data to "Days since 001-01-01 00:00:00 UTC".

	mesoData = {}
	f = netCDF4.Dataset(ncfile, 'r')
	fc = f.variables['fc'][:]    #Coriolis Factor
	date_f = ncDatesOffset + ncDatesSlope * f.variables['time'][:]
	idates = np.logical_and(date_f >= our_date2num(datefrom), date_f <= our_date2num(dateto))
	print(("ncfile dates coverage: " + str(our_num2date(date_f[0])) + " to " + str(our_num2date(date_f[-1]))))
	
	mesoData['t'] = date_f[idates]
	mesoData['z'] = f.variables['z'][:]
	mesoData['hrs_since_t0'] = 24.0*(mesoData['t'] - mesoData['t'][0])  # time in hours since datefrom
	
	for ii,vv in enumerate(vars2extract):
		ncvar = vv[0]
		if (ncvar in f.variables):
			if (len(f.variables[ncvar].shape) == 1):
				mesoData[ncvar] = f.variables[ncvar][:][idates]
			elif (len(f.variables[ncvar].shape) == 2):
				mesoData[ncvar] = f.variables[ncvar][:][idates,:]
				if plotData:
					plt.figure(ii); 
					plot_nc_tzData(ncvar, mesoData['t'], mesoData['z'], mesoData[ncvar], datefrom, dateto)
			else:
				print ("so far it is only prepared for 1d or 2d nc variables")
			
			if vv[1]:
				mesoData[ncvar] = mesoData[ncvar] * fc
				
		else:
			print((ncvar + " does not exist in the netcdf file"))
	
	if ( ('Ug' in mesoData) & ('Vg' in mesoData) ):
		# Add Geostrophic wind velocity [m s-1] and direction [deg]
		mesoData['Sg'] = (mesoData['Ug']**2 + mesoData['Vg']**2)**0.5
		mesoData['Dg'] = 180 + np.arctan2(mesoData['Ug'],mesoData['Vg'])*180.0/np.pi
	
	f.close()
	
	return mesoData


def plot_nc_tzData(varname, t, z, tzData, datefrom, dateto, relativeTime=False):
	# plots generic t,z data
	zlim = 4000     # maximum height in z axis
	Nm = 6     # every 3 hrs

	hoursFmt = mdates.DateFormatter('%H')
	Nticks = int((our_date2num(dateto)-our_date2num(datefrom))*24/Nm + 1)
	tticks = [datefrom + datetime.timedelta(hours = Nm*x) for x in range(Nticks)]
	
	if not relativeTime:
		t = our_num2date(t)
 
	cs = plt.contourf(t, z, tzData.T, len(z), cmap=plt.cm.Spectral_r)
	#cs = plt.contourf(t, z, tzData.T, len(z), cmap=plt.cm.jet)
	plt.ylabel('$z$ [$m$]');	plt.title(varname)
	ax = plt.gca();	ax.set_ylim([0,zlim])
	plt.xlabel('hrs since '+ datefrom.strftime('%Y-%m-%d %H:%M'));
	plt.grid(which='major',color='grey',linestyle=':')
	
	if not relativeTime:
		ax.xaxis.set_major_formatter(hoursFmt)
		plt.xticks(tticks)
		plt.colorbar(cs, orientation='vertical')
	
	plt.show()
	sys.stdout.flush()


def our_date2num(dateObject):
	# This function defines our common date2num format to use in all functions
	return netCDF4.date2num(dateObject, units='days since 0001-01-01 00:00:00.0', calendar='gregorian')

def our_num2date(dateNumObject):
	# This function defines our common date2num format to use in all functions
	return netCDF4.num2date(dateNumObject, units='days since 0001-01-01 00:00:00.0', calendar='gregorian')
