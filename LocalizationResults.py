#!/usr/bin/env python

import numpy
import PublicTableXMLTools
import fileinput
import math
import matplotlib.pylab as plt
import code
import IDL
##########################################################################################

def AngularDistance(RA1, DEC1, RA2, DEC2):

	try:
		RA1 = float(RA1)
		DEC1 = float(DEC1)
		RA2 = float(RA2)
		DEC2 = float(DEC2)

		# Converting everything to radians
		ra1rad = RA1 * math.pi/180.
		ra2rad = RA2 * math.pi/180.
		dec1rad = DEC1 * math.pi/180.
		dec2rad = DEC2 * math.pi/180.

		# Calculate scalar product for determination of angular separation      
		x=math.cos(ra1rad)*math.cos(dec1rad)*math.cos(ra2rad)*math.cos(dec2rad)
		y=math.sin(ra1rad)*math.cos(dec1rad)*math.sin(ra2rad)*math.cos(dec2rad)
		z=math.sin(dec1rad)*math.sin(dec2rad)
		rad=math.acos(x+y+z)

		# Use Pythargoras approximation if rad < 1 arcsec
		if rad<0.000004848:
			rad=math.sqrt((math.cos(dec1rad)*(ra1rad-ra2rad))**2+(dec1rad-dec2rad)**2)
			pass

		# Angular separation in degrees
		Angle = rad*180/math.pi

		return Angle

	except Exception, message:
		print message

		return float('nan')


##########################################################################################

def Parse(xmlfile, ProjectDirectory='/Users/kocevski/Research/Analysis/LocalizationStudies/', version='ver1'):
	"""
	Usage Example:
	function(Variable=True)
	"""

	import numpy
	import PublicTableXMLTools
	import fileinput
	import math

	# Define the default project directory
	P7_P203_Directory = "%s/Likelihood/%s/P7_P203/" % (ProjectDirectory, version)
	P8_P301_Directory = "%s/Likelihood/%s/P8_P301/" % (ProjectDirectory, version)

	# Load the xml file
	xml = PublicTableXMLTools.xml(xmlfile)

	# Extract a dictionary containing all of the GRB information
	xmlData = xml.ExtractData()

	# Extract the GRB information from the xml file
	GRBs = xmlData['GRBNAME']
	METs = xmlData['MET']
	RAs = xmlData['LATRA']
	DECs = xmlData['LATDEC']
	LIKESTARTs = xmlData['LIKESTART']	
	LIKESTOPs = xmlData['LIKESTOP']
	TSs = xmlData['TS']
	THETAs = xmlData['THETA']
	ZENITHs = xmlData['ZENITH']
	LLE = xmlData['LLESIGMA']
	TS = xmlData['TS']
	BestRAs = xmlData['RA']
	BestDECs = xmlData['DEC']
	BestErrors = xmlData['ERROR']
	BestSource = xmlData['POSITIONSOURCE']
	LIKE = xmlData['ABOVE75MEVDETECTION']

	# Create empty arrays to store the localization information
	gtfindsrc_P7_P203_RA = []
	gtfindsrc_P7_P203_Dec =[]
	gtfindsrc_P7_P203_Error =[]

	gtfindsrc_P8_P301_RA = []
	gtfindsrc_P8_P301_Dec = []
	gtfindsrc_P8_P301_Error = []

	tsmap_P7_P203_MaxTS = numpy.array([],dtype='float')
	tsmap_P7_P203_RA = numpy.array([],dtype='float')
	tsmap_P7_P203_Dec = numpy.array([],dtype='float')
	tsmap_P7_P203_Error68 = numpy.array([],dtype='float')
	tsmap_P7_P203_Error90 = numpy.array([],dtype='float')
	tsmap_P7_P203_Error95 = numpy.array([],dtype='float')

	tsmap_P8_P301_MaxTS	= numpy.array([],dtype='float')
	tsmap_P8_P301_RA = numpy.array([],dtype='float')
	tsmap_P8_P301_Dec = numpy.array([],dtype='float')
	tsmap_P8_P301_Error68 =numpy.array([],dtype='float')
	tsmap_P8_P301_Error90 =numpy.array([],dtype='float')
	tsmap_P8_P301_Error95 =numpy.array([],dtype='float')


	# Loop through each GRB and extract the localizaton information.  If localization informaton isn't found, return nan values 
	counter = 1
	failedGRBs =[]
	passedGRBs = []
	for grb, met, ra, dec, tmin, tmax, ts, theta, zenith in zip(GRBs, METs, RAs, DECs, LIKESTARTs, LIKESTOPs, TSs, THETAs, ZENITHs):
		
		failed = False

		gtfindsrc_P7_P203_File = "%s/%s/gtfindsrc_%s.txt" % (P7_P203_Directory, grb, grb)
		try:
			Ra, Dec, AngularSeperation, Error = numpy.loadtxt(gtfindsrc_P7_P203_File, unpack=True, comments='#', dtype = str)
			firstBadRow = numpy.where(Ra == 'initial')[0][0]
			gtfindsrc_P7_P203_RA.append(Ra[firstBadRow-1])
			gtfindsrc_P7_P203_Dec.append(Dec[firstBadRow-1])
			gtfindsrc_P7_P203_Error.append(Error[firstBadRow-1])
		except Exception, message:
			#print message
			failed = True
			gtfindsrc_P7_P203_RA.append(float('nan'))
			gtfindsrc_P7_P203_Dec.append(float('nan'))
			gtfindsrc_P7_P203_Error.append(float('nan'))		

		gtfindsrc_P8_P301_File = "%s/%s/gtfindsrc_%s.txt" % (P8_P301_Directory, grb, grb)
		try:
			Ra, Dec, AngularSeperation, Error = numpy.loadtxt(gtfindsrc_P8_P301_File, unpack=True, comments='#', dtype = str)
			firstBadRow = numpy.where(Ra == 'initial')[0][0]
			gtfindsrc_P8_P301_RA.append(Ra[firstBadRow-1])
			gtfindsrc_P8_P301_Dec.append(Dec[firstBadRow-1])
			gtfindsrc_P8_P301_Error.append(Error[firstBadRow-1])
		except Exception, message:
			#print message
			failed = True			
			gtfindsrc_P8_P301_RA.append(float('nan'))
			gtfindsrc_P8_P301_Dec.append(float('nan'))
			gtfindsrc_P8_P301_Error.append(float('nan'))


		tsmap_P7_P203_File = "%s/%s/tsmap_%s.txt" % (P7_P203_Directory, grb, grb)
		try:
			for line in fileinput.input([tsmap_P7_P203_File]):
				if 'Max TS:' in line:
					tsmap_P7_P203_MaxTS = numpy.append(tsmap_P7_P203_MaxTS, line.split()[2])
				if 'Ra Dec:' in line:
					tsmap_P7_P203_RA = numpy.append(tsmap_P7_P203_RA, line.split()[2])
					tsmap_P7_P203_Dec = numpy.append(tsmap_P7_P203_Dec, line.split()[3])
				if '68 percent' in line:
					tsmap_P7_P203_Error68 = numpy.append(tsmap_P7_P203_Error68, line.split()[4])
				if '90 percent' in line:
					tsmap_P7_P203_Error90 = numpy.append(tsmap_P7_P203_Error90, line.split()[4])
				if '95 percent' in line:
					tsmap_P7_P203_Error95 = numpy.append(tsmap_P7_P203_Error95, line.split()[4])
			fileinput.close()		
		except Exception, message:
			#print message
			failed = True			
			tsmap_P7_P203_MaxTS = numpy.append(tsmap_P7_P203_MaxTS, float('nan'))
			tsmap_P7_P203_RA = numpy.append(tsmap_P7_P203_RA, float('nan'))
			tsmap_P7_P203_Dec = numpy.append(tsmap_P7_P203_Dec, float('nan'))
			tsmap_P7_P203_Error68 = numpy.append(tsmap_P7_P203_Error68, float('nan'))
			tsmap_P7_P203_Error90 = numpy.append(tsmap_P7_P203_Error90, float('nan'))
			tsmap_P7_P203_Error95 = numpy.append(tsmap_P7_P203_Error95, float('nan'))


		tsmap_P8_P301_File = "%s/%s/tsmap_%s.txt" % (P8_P301_Directory, grb, grb)
		try:
			for line in fileinput.input([tsmap_P8_P301_File]):
				if 'Max TS:' in line:
					tsmap_P8_P301_MaxTS = numpy.append(tsmap_P8_P301_MaxTS, line.split()[2])
				if 'Ra Dec:' in line:
					tsmap_P8_P301_RA = numpy.append(tsmap_P8_P301_RA, line.split()[2])
					tsmap_P8_P301_Dec = numpy.append(tsmap_P8_P301_Dec, line.split()[3])
				if '68 percent' in line:
					tsmap_P8_P301_Error68 = numpy.append(tsmap_P8_P301_Error68, line.split()[4])
				if '90 percent' in line:
					tsmap_P8_P301_Error90 = numpy.append(tsmap_P8_P301_Error90, line.split()[4])
				if '95 percent' in line:
					tsmap_P8_P301_Error95 = numpy.append(tsmap_P8_P301_Error95, line.split()[4])
			fileinput.close()	
		except Exception, message:
			#print message
			failed = True			
			tsmap_P8_P301_MaxTS = numpy.append(tsmap_P8_P301_MaxTS, float('nan'))
			tsmap_P8_P301_RA = numpy.append(tsmap_P8_P301_RA, float('nan'))
			tsmap_P8_P301_Dec = numpy.append(tsmap_P8_P301_Dec, float('nan'))
			tsmap_P8_P301_Error68 = numpy.append(tsmap_P8_P301_Error68, float('nan'))
			tsmap_P8_P301_Error90 = numpy.append(tsmap_P8_P301_Error90, float('nan'))
			tsmap_P8_P301_Error95 = numpy.append(tsmap_P8_P301_Error95, float('nan'))

		if failed == True:
			failedGRBs.append(grb)
		else:
			passedGRBs.append(grb)

	# Print the bursts that passed and failed
	print "\nGood Bursts:"
	for grb in passedGRBs:
		print grb

	print "\nBad Bursts:"
	for grb in failedGRBs:
		print grb

		
	print "\nLikelihood Results: %s" % version
	print "GRBs		LLE	CatTS	P7MaxTS	P8MaxTS	P8RA	P8Dec	P8Err90"
	for grb, lle, ts, P7ts, p8ts, ra, dec, error in zip(GRBs, LLE, TS, tsmap_P7_P203_MaxTS, tsmap_P8_P301_MaxTS, tsmap_P8_P301_RA, tsmap_P8_P301_Dec, tsmap_P8_P301_Error90):
		print "%s	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	+/- %.3f" % (grb, float(lle), float(ts), float(P7ts), float(p8ts), float(ra), float(dec), float(error))

	#Make sure the arrays are floats
	tsmap_P7_P203_MaxTS = numpy.array(tsmap_P7_P203_MaxTS).astype(float)
	tsmap_P7_P203_Error90 = numpy.array(tsmap_P7_P203_Error90).astype(float)
	tsmap_P7_P203_Error95 = numpy.array(tsmap_P7_P203_Error95).astype(float)

	tsmap_P8_P301_MaxTS = numpy.array(tsmap_P8_P301_MaxTS).astype(float)
	tsmap_P8_P301_Error90 = numpy.array(tsmap_P8_P301_Error90).astype(float)
	tsmap_P8_P301_Error95 = numpy.array(tsmap_P8_P301_Error95).astype(float)


	# Set the plotting format
	try:
		IDL.plotformat()

 
	# Plot P7_P203 TS vs P8_P301 TS
	GRB130427A = numpy.where(GRBs == '130427324')
	good = numpy.where((tsmap_P7_P203_MaxTS > 0) & (tsmap_P8_P301_MaxTS > 0))[0]
	plt.scatter(tsmap_P7_P203_MaxTS[good], tsmap_P8_P301_MaxTS[good])
	plt.annotate('GRB130427A', xy=(tsmap_P7_P203_MaxTS[GRB130427A],tsmap_P8_P301_MaxTS[GRB130427A]), xytext=(-35,10), textcoords='offset points', ha='center', va='bottom')
	plt.plot([1,10000],[1,10000], '--')
	plt.xscale('log')
	plt.yscale('log')
	plt.xlim(1,10000)
	plt.ylim(1,10000)
	plt.xlabel('TS (P7_203)')
	plt.ylabel('TS (P8_301)')	
	plt.show()

	# Plot P7_P203 90% Error vs P8_P301 90% Error
	good = numpy.where((tsmap_P7_P203_MaxTS > 0) & (tsmap_P8_P301_MaxTS > 0))[0]
	plt.scatter(tsmap_P7_P203_Error90[good], tsmap_P8_P301_Error90[good], c=numpy.log10(tsmap_P8_P301_MaxTS[good]))
	plt.annotate('GRB130427A', xy=(tsmap_P7_P203_Error90[GRB130427A],tsmap_P8_P301_Error90[GRB130427A]), xytext=(40,-10), textcoords='offset points', ha='center', va='bottom')
	cbar = plt.colorbar(pad = 0.02)
	cbar.set_label(r'log TS$_{\rm P8}$')
	plt.plot([0.001,10],[0.001,10], '--')
	plt.xlim(0.01,10)
	plt.ylim(0.01,10)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'$\sigma_{90\%}$ (P7_203)')
	plt.ylabel(r'$\sigma_{90\%}$ (P8_301)')
	plt.show()

	# Calculate the angular seperation between P7 and P8
	import AngularSeparation
	angularSeperation_P7toP8 = numpy.array([])
	for ra1, dec1, ra2, dec2 in zip(tsmap_P7_P203_RA, tsmap_P7_P203_Dec, tsmap_P8_P301_RA, tsmap_P8_P301_Dec):
		angle = AngularSeparation.degrees(ra1, dec1, ra2, dec2)
		angularSeperation_P7toP8 = numpy.append(angularSeperation_P7toP8, angle)
		pass

	# Create a subselection of all bursts and bursts with TS > 25
	good_All = numpy.where(tsmap_P8_P301_Error90 > 0)
	good_HighTS = numpy.where(tsmap_P8_P301_MaxTS > 25)

	# Obtain the number of bursts that survived the cuts
	numberOfBurstsP8_90_All = len(angularSeperation_P7toP8[good_All]/tsmap_P8_P301_Error90[good_All])
	numberOfBurstsP8_95_All = len(angularSeperation_P7toP8[good_All]/tsmap_P8_P301_Error95[good_All])
	numberOfBurstsP8_90_HighTS = len(angularSeperation_P7toP8[good_HighTS]/tsmap_P8_P301_Error90[good_HighTS])
	numberOfBurstsP8_95_HighTS = len(angularSeperation_P7toP8[good_HighTS]/tsmap_P8_P301_Error95[good_HighTS])

	# Generate a cumulative sum array
	cumulativeSumP8_90_All = numpy.arange(numberOfBurstsP8_90_All)/(float(numberOfBurstsP8_90_All)-1)
	cumulativeSumP8_95_All = numpy.arange(numberOfBurstsP8_95_All)/(float(numberOfBurstsP8_95_All)-1)

	# Plot the angular seperation between P7 and P8 normalized by the P8_P301 90% Error
	plt.step(numpy.sort(angularSeperation_P7toP8[good_All]/tsmap_P8_P301_Error90[good_All]), cumulativeSumP8_90_All)
	plt.step(numpy.sort(angularSeperation_P7toP8[good_All]/tsmap_P8_P301_Error95[good_All]), cumulativeSumP8_95_All)
	i = numpy.where(numpy.sort(angularSeperation_P7toP8[good_All]/tsmap_P8_P301_Error90[good_All]) == angularSeperation_P7toP8[GRB130427A]/tsmap_P8_P301_Error90[GRB130427A])
	plt.scatter(angularSeperation_P7toP8[GRB130427A]/tsmap_P8_P301_Error90[GRB130427A], cumulativeSumP8_90_All[i])
	plt.annotate('GRB130427A', xy=(angularSeperation_P7toP8[GRB130427A]/tsmap_P8_P301_Error90[GRB130427A], cumulativeSumP8_90_All[i] ), xytext=(-20,-20), textcoords='offset points', ha='center', va='bottom')
	plt.legend(('P8_301 90% C.L.', 'P8_301 95% C.L.'), frameon=False, scatterpoints=1, loc=4)
	plt.plot([1,1],[0,1.05], '--')
	plt.ylim(0,1.05)
	plt.xlim(0,5)
	plt.xlabel(r'$\theta_{\rm P7 to P8} / \sigma$')	
	plt.show()

	# Select a subset of bursts with known x-ray, optical, or radio localizations
	#good = numpy.where((LIKE == 1) & (BestSource != 'Fermi-LAT') & (BestSource != 'Fermi-GBM') & (BestSource != 'IPN'))[0]
	good = numpy.where((tsmap_P8_P301_Error95 > 0) & (BestSource != 'Fermi-LAT') & (BestSource != 'Fermi-GBM') & (BestSource != 'IPN'))[0]
	GRB130427A = numpy.where(GRBs[good] == '130427324')

	# Find the angular seperation between P7 and the best localizaton
	angularSeperation_BestToP7 = numpy.array([])
	for ra1, dec1, ra2, dec2 in zip(BestRAs[good], BestDECs[good], tsmap_P7_P203_RA[good], tsmap_P7_P203_Dec[good]):
		angle = AngularSeparation.degrees(ra1, dec1, ra2, dec2)
		angularSeperation_BestToP7 = numpy.append(angularSeperation_BestToP7,angle)
		pass

	# Normalized the angular seperation between P7 and the best localizaton by the P7_P203 90% and 95% Error
	normalizedAngularSeperationP7_90 = angularSeperation_BestToP7/tsmap_P7_P203_Error90[good]
	normalizedAngularSeperationP7_95 = angularSeperation_BestToP7/tsmap_P7_P203_Error95[good]

	# Find the angular seperation between P8 and the best localizaton
	angularSeperation_BestToP8 = numpy.array([])
	for ra1, dec1, ra2, dec2 in zip(BestRAs[good], BestDECs[good], tsmap_P8_P301_RA[good], tsmap_P8_P301_Dec[good]):
		angle = AngularSeparation.degrees(ra1, dec1, ra2, dec2)
		angularSeperation_BestToP8 = numpy.append(angularSeperation_BestToP8,angle)
		pass

	# Normalized the angular seperation between P78 and the best localizaton by the P8_P301 90% and 95% Error
	normalizedAngularSeperationP8_90 = angularSeperation_BestToP8/tsmap_P8_P301_Error90[good]
	normalizedAngularSeperationP8_95 = angularSeperation_BestToP8/tsmap_P8_P301_Error95[good]

	print '\nTrue Position to LAT Localization Comparison (P7_203)'
	print "GRB		BestRa	BestDec	BestError  		angle	LatRa	LatDec	LatErr90 "
	for grb, bestra, bestdec, besterror, angle, latra, latdec, laterror in zip(GRBs[good], BestRAs[good], BestDECs[good], BestErrors[good], angularSeperation_BestToP7, tsmap_P7_P203_RA[good], tsmap_P7_P203_Dec[good], tsmap_P7_P203_Error90[good]):
		print "%s	%.2f	%.2f	+/- %.3f 		%.2f	%.2f 	%.2f 	+/- %.3f" % (grb, float(bestra), float(bestdec), float(besterror), float(angle), float(latra), float(latdec), float(laterror))

	print '\nTrue Position to LAT Localization Comparison (P8_P301)'
	print "GRB		BestRa	BestDec	BestError  		angle	LatRa	LatDec	LatErr90 "	
	for grb, bestra, bestdec, besterror, angle, latra, latdec, laterror in zip(GRBs[good], BestRAs[good], BestDECs[good], BestErrors[good], angularSeperation_BestToP8, tsmap_P8_P301_RA[good], tsmap_P8_P301_Dec[good], tsmap_P8_P301_Error90[good]):
		print "%s	%.2f	%.2f	+/- %.3f 		%.2f	%.2f 	%.2f 	+/- %.3f" % (grb, float(bestra), float(bestdec), float(besterror), float(angle), float(latra), float(latdec), float(laterror))

	# Make sure all values are finite (no inf!)
	finiteP7_90 = numpy.isfinite(numpy.array(normalizedAngularSeperationP7_90).astype(float))
	finiteP7_95 = numpy.isfinite(numpy.array(normalizedAngularSeperationP7_95).astype(float))
	finiteP8_90 = numpy.isfinite(numpy.array(normalizedAngularSeperationP8_90).astype(float))
	finiteP8_95 = numpy.isfinite(numpy.array(normalizedAngularSeperationP8_95).astype(float))

	# Find the number of bursts that survived the cuts
	numberOfBurstsP7_90 = len(normalizedAngularSeperationP7_90[finiteP7_90])
	numberOfBurstsP7_95 = len(normalizedAngularSeperationP7_95[finiteP7_95])
	numberOfBurstsP8_90 = len(normalizedAngularSeperationP8_90[finiteP8_90])
	numberOfBurstsP8_95 = len(normalizedAngularSeperationP8_95[finiteP8_95])

	# Generate a cumulative sum arrays
	cumulativeSumP7_90 = numpy.arange(numberOfBurstsP7_90)/(float(numberOfBurstsP7_90)-1)
	cumulativeSumP7_95 = numpy.arange(numberOfBurstsP7_95)/(float(numberOfBurstsP7_95)-1)
	cumulativeSumP8_90 = numpy.arange(numberOfBurstsP8_90)/(float(numberOfBurstsP8_90)-1)
	cumulativeSumP8_95 = numpy.arange(numberOfBurstsP8_95)/(float(numberOfBurstsP8_95)-1)
	i = numpy.where(numpy.sort(normalizedAngularSeperationP8_90[finiteP8_90] == normalizedAngularSeperationP8_90[GRB130427A]))

	# Plot the angular seperation between P7 and P8 localizations and the best known position normalized by their respective errors
	plt.step(numpy.sort(normalizedAngularSeperationP7_90[finiteP7_90]), cumulativeSumP7_90)
	plt.step(numpy.sort(normalizedAngularSeperationP7_95[finiteP7_95]), cumulativeSumP7_95)
	plt.step(numpy.sort(normalizedAngularSeperationP8_90[finiteP8_90]), cumulativeSumP8_90)
	plt.step(numpy.sort(normalizedAngularSeperationP8_95[finiteP8_95]), cumulativeSumP8_95)
	plt.scatter(numpy.sort(normalizedAngularSeperationP8_90[GRB130427A]), cumulativeSumP8_90[i])
	plt.annotate('GRB130427A', xy=(numpy.sort(normalizedAngularSeperationP8_90[GRB130427A]), cumulativeSumP8_90[i]), xytext=(0,-50), textcoords='offset points', ha='center', va='bottom')

	# Setup the plot
	plt.legend(('P7_203 90% C.L.', 'P7_203 95% C.L.', 'P8_301 90% C.L.', 'P8_301 95% C.L.'), frameon=False, scatterpoints=1, loc=4)
	plt.plot([1,1],[0,1.05], '--')
	plt.ylim(0,1.05)
	plt.xlim(0,4)
	plt.xlabel(r'$\Delta \theta_{\rm True} / \sigma$')	
	plt.show()

	#code.interact(local=locals())

	# Plot the normalized angular seperation between the P7 and P8 localizations and the best known position vs the P8 TS
	plt.scatter(normalizedAngularSeperationP7_90[finiteP8_90], normalizedAngularSeperationP8_90[finiteP8_90], c=numpy.log(tsmap_P8_P301_MaxTS[good][finiteP8_90]))
	plt.plot([0.001,100],[0.001,100], '--')
	cbar = plt.colorbar(pad = 0.02)
	cbar.set_label(r'log TS$_{\rm P8}$')
	plt.ylim(0,5)
	plt.xlim(0,5)
	plt.xlabel(r'$\Delta \theta_{\rm True->P7} / \sigma_{\rm P7}$')	
	plt.ylabel(r'$\Delta \theta_{\rm True->P8} / \sigma_{\rm P8}$')		
	plt.show()

	#return tsmap_P7_P203_MaxTS, tsmap_P8_P301_MaxTS
	return

##########################################################################################
		
if __name__=='__main__':

	# Check to see if any arguments were passed
	if(len(sys.argv) > 1):

		# Extract the first arugument
		Argument = sys.argv[1]
							
				
	# If no arguments were passed, tell the user how to use the script		
	else:
		print 'Usage: MyProject.py Argument'
		sys.exit(0)


	function()		
