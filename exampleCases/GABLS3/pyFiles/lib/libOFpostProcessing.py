#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Auxiliary functions for the postprocessing of the NEWAFoam sets
Feb2018 Roberto Chavez rchavez@cener.com
'''
import numpy as np
import re
import pandas as pd
import os

def findLastIteration(casePath):
	intFolders = list()
	for o in os.listdir(casePath):
		if os.path.isdir(os.path.join(casePath,o)):
			try:
				intFolders.append(int(o))
			except:
				continue	
	lastIter = max(intFolders)
	return (lastIter,intFolders)

#%
def readSetFile(ofPath,setFileName,variab,outputTimeCoords=False):
	
	setsPath = ofPath+'/postProcessing/sets/'
	[lastIt,iteras] = findLastIteration(setsPath)
	
	t=np.sort(np.array(iteras))
		
	# the issue with sets is that the string is divided with "_" to we need to hack p_rgh
	is_prgh = setFileName.find('p_rgh')
	
	setFileNameNew = setFileName
	if is_prgh > 0:
		setFileNameNew = setFileName[0:is_prgh] + 'prgh' + setFileName[is_prgh+5:]
	
	if (variab=='p_rgh'):
		variab='prgh'
		
	temp=re.split(r'[_.]\s*',setFileNameNew)

	# IM DISCARDING FIRST TIME (0) because it doesn't contain many of the computed variables (only initial conditions)
	#t=np.delete(t,0)

	fileName=setsPath+'/'+str(t[0])+'/'+setFileName
	ndata=np.genfromtxt(fileName)		
	coords = pd.DataFrame( dict( x=ndata[:,0], y=ndata[:,1], z=ndata[:,2]) )
	
	if (variab=='U'):
		ux=pd.DataFrame()
		uy=pd.DataFrame()
		uz=pd.DataFrame()

		for ii in t:	
			print 'Reading t = '+str(ii)+'s'		
			fileName=setsPath+'/'+str(ii)+'/'+setFileName

			ndata=np.genfromtxt(fileName)
			ux['t'+str(ii)]=ndata[:,3]
			uy['t'+str(ii)]=ndata[:,4]
			uz['t'+str(ii)]=ndata[:,5]					
		
		res=dict(ux=ux, uy=uy, uz=uz)
	else:
		params=temp[1:-1]
		print "Available parameters: "
		print params

		ipar = 0
		sflag = False
		while sflag==False:
			sflag = params[ipar]==variab
			ipar = ipar +1
		
		ipar=ipar-1
		
		data = pd.DataFrame()
		
		for ii in t:
			print 'Reading t = '+str(ii)+'s'
			fileName=setsPath+'/'+str(ii)+'/'+setFileName
			
			ndata=np.genfromtxt(fileName)		
			data['t'+str(ii)]=ndata[:,ipar+3]
		
		res=dict(data=data)
	if (outputTimeCoords==True):
		return res, t, coords
	else:
		return res
