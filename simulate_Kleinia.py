#!/usr/bin/python3

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

import random
import os
import math
import shlex, subprocess
import numpy as np

#Function to transform simulations from the ms format into NumPy arrays
def ms2nparray(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	f = []
	for i in k:
		L = g[i+4:i+nDNANsam+4]
		q = []
		for i in L:
			i = i = [int(j) for j in list(i.decode('utf-8'))]
			i = np.array(i)
			q.append(i)
		q = np.array(q)
		q = q.astype("int8")
		f.append(np.array(q))
	return f

### variable declarations

#define the number of simulations
Priorsize = 10000

# Sample sizes (times 2 for diploids)
## Sample size of Fuerteventura+Lanzarote (25 samples)
nFL = 25*2
## Sample size of Gran Canaria (8 samples)
nGC = 8*2
## Sample size of Oriental Tenerife  (10 samples)
nTor = 10*2
## Sample size of Ocidental Tenerife  (112 samples)
nToc = 12*2
## Sample size of La Gomera (10 samples)
nLG = 10*2
## Sample size of La Palma (9 samples)
nLP = 9*2
## Sample size of El Hierro (10 samples)
nEH = 10*2

## Sample sample size of all lineages.
nDNANsam = nFL + nGC + nTor + nToc + nLG + nLP + nEH

# Create empty lists to store the simulations from each scenario
Mod_SS = []
Mod_SSHclim = []
Mod_CIH1 = []
Mod_CIH2 = []

## create files to store parameters from each scenario
parSS = open("parSS.txt","w")
parSSHclim = open("parSSHclim.txt","w")
parCIH1 = open("parCIH1.txt","w")
parCIH2 = open("parCIH2.txt","w")

## SCENARIO 1: STEPPING STONES with Continent as colonizer 
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	
	## migration prior set to 0 in this model.
	m=0

	## divergence time prior in years, following uniform distributions.
	#Tm is not used in this model (assign default value).
	Tm1=0
	Tm2=0
	T6=random.uniform(1190000,4000000)
	T5=random.uniform(0,T6)
	T4=random.uniform(0,T5)
	T3=0
	T2=random.uniform(0,T4)
	T1=random.uniform(0,T2)

	## Transform to coalescent units -> Time in years/(number of years per generation*4Ne)
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=0
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate according to the formula in the ms manual
	GrowthT1 = -(1/T6)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T5)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T4)*math.log((1/NTor)/(1/FoundedSizeRatio))
	#Divergence time for both parts of Tenerife is 0
	GrowthT4 = 0
	GrowthT5 = -(1/T2)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
	
	## simulate SNPs by plugging the parameter values into a ms command
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 7 %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -g 1 %f -g 2 %f -g 5 %f -g 6 %f -g 7 %f -ej 0 4 3 -eg 0 3 %f -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 3 %f -ej %f 5 3 -en %f 2 %f -ej %f 3 2 -en %f 1 %f -ej %f 2 1" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLP, nLG, nEH, NGC, NTor, NToc, NLP, NLG, NEH, GrowthT1, GrowthT2, GrowthT5, GrowthT6, GrowthT7, GrowthT3, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6), shell=True, stdout=subprocess.PIPE).stdout
	#read the ms output
	output = com.read().splitlines()
	#Convert the ms output into a NumPy array
	Mod_SS.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parSS.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm1, Tm2, T1, T2, T3, T4, T5, T6, Ne, NGC, NTor, NToc, NLP, NLG, NEH, m, FoundedSizeRatio))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))

#Save NumPy array containing the simulated SNPs
Mod_SS=np.array(Mod_SS)
np.save('Mod_SS.npy', Mod_SS)
del(Mod_SS)

## MODEL 2: Surfing Syngaemon Hypothesis with migration associated with climate
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	
	## migration prior from 0 to 10 in this model.
	m=random.uniform(0,10)

	## divergence time prior in years, following uniform distributions.
	#T8 is not used in this model (assign default value).
	Tm2=random.uniform(14000,1250000)
	Tm1=random.uniform(14000,Tm2)
	T6=random.uniform(1190000,4000000)
	#make sure T6 is not younger than Tm2
	if T6<Tm2:
		T6=Tm2
	T5=random.uniform(Tm2,T6)
	T4=random.uniform(Tm2,T5)
	T3=0
	T2=random.uniform(Tm2,T4)
	T1=random.uniform(Tm2,T2)

	## Transform to coalescent units
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=0
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)
	coalTm2=Tm2/(genlen*4.0*Ne)
	coalTm1=Tm1/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT1 = -(1/T6)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T5)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T4)*math.log((1/NTor)/(1/FoundedSizeRatio))
	#Divergence time for both parts of Tenerife is 0
	GrowthT4 = 0
	GrowthT5 = -(1/T2)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))
	
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 7 %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -g 1 %f -g 2 %f -g 5 %f -g 6 %f -g 7 %f -ej 0 4 3 -eg 0 3 %f -em %f 1 2 %f -em %f 2 3 %f -em %f 1 2 0 -em %f 2 3 0 -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 3 %f -ej %f 5 3 -en %f 2 %f -ej %f 3 2 -en %f 1 %f -ej %f 2 1" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLP, nLG, nEH, NGC, NTor, NToc, NLP, NLG, NEH, GrowthT1, GrowthT2, GrowthT5, GrowthT6, GrowthT7, GrowthT3, coalTm1, m, coalTm1, m, coalTm2, coalTm2, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_SSHclim.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parSSHclim.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm1, Tm2, T1, T2, T3, T4, T5, T6, Ne, NGC, NTor, NToc, NLP, NLG, NEH, m, FoundedSizeRatio))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))

Mod_SSHclim=np.array(Mod_SSHclim)
np.save('Mod_SSHclim.npy', Mod_SSHclim)
del(Mod_SSHclim)

## MODEL 3: Central Island Hypothesis (CIH) with Tenerife populations combined
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	NC = random.uniform(0.2,2)

	## migration prior set to 0 in this model.
	m=0
	
	## divergence time prior in years, following uniform distributions.
	Tm1=0
	Tm2=0
	T6=random.uniform(1190000,4000000)
	T5=random.uniform(0,T6)
	T4=random.uniform(0,T5)
	T3=0
	T2=random.uniform(0,T4)
	T1=random.uniform(0,T2)

	## Transform to coalescent units
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=0
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT1 = -(1/T4)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T4)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T6)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = 0
	GrowthT5 = -(1/T2)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
		
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 7 %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -g 1 %f -g 2 %f -g 5 %f -g 6 %f -g 7 %f -ej 0 4 3 -eg 0 3 %f -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 2 %f -ej %f 2 1 -ej %f 5 1 -en %f 3 %f -en %f 1 %f -ej %f 1 3" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLG, nLP, nEH, NGC, NTor, NToc, NLG, NLP, NEH, GrowthT1, GrowthT2, GrowthT5, GrowthT6, GrowthT7, GrowthT3, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT4, FoundedSizeRatio, coalT4, coalT5, coalT6, FoundedSizeRatio, coalT6, FoundedSizeRatio, coalT6), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_CIH1.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parCIH1.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm1, Tm2, T1, T2, T3, T4, T5, T6, Ne, NGC, NTor, NToc, NLP, NLG, NEH, m, FoundedSizeRatio))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

Mod_CIH1=np.array(Mod_CIH1)
np.save('Mod_CIH1.npy', Mod_CIH1)
del(Mod_CIH1)

## MODEL 4: Central Island Hypothesis (CIH) with Tenerife populations separated
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)

	## migration prior set to 0 in this model.
	m=0
	
	## divergence time prior in years, following uniform distributions.
	Tm1=0
	Tm2=0
	T6=random.uniform(1190000,4000000)
	T5=random.uniform(0,T6)
	T4=random.uniform(0,T5)
	T3=random.uniform(0,T4)
	T2=random.uniform(0,T3)
	T1=random.uniform(0,T2)

	## Transform to coalescent units
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=T3/(genlen*4.0*Ne)
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT1 = -(1/T6)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T3)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T5)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = -(1/T4)*math.log((1/NToc)/(1/FoundedSizeRatio))
	GrowthT5 = -(1/T2)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
			
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 7 %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -g 1 %f -g 2 %f -g 3 %f -g 4 %f -g 5 %f -g 6 %f -g 7 %f -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 2 %f -ej %f 2 1 -en %f 4 %f -ej %f 5 4 -en %f 3 %f -ej %f 1 3 -en %f 1 %f -ej %f 4 3" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLG, nLP, nEH, NGC, NTor, NToc, NLG, NLP, NEH, GrowthT1, GrowthT2, GrowthT3, GrowthT4, GrowthT5, GrowthT6, GrowthT7, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT3, FoundedSizeRatio, coalT3, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_CIH2.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values
	parCIH2.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm1, Tm2, T1, T2, T3, T4, T5, T6, Ne, NGC, NTor, NToc, NLP, NLG, NEH, m, FoundedSizeRatio))
	print("Completed %d %% of Model 4 simulations" % (float(i)/Priorsize*100))

Mod_CIH2=np.array(Mod_CIH2)
np.save('Mod_CIH2.npy', Mod_CIH2)
del(Mod_CIH2)