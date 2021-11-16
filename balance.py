import glob 
import numpy as np
import math
import os
import time
from decimal import *

##Number of electrons you want at the end
ENDNUM=2e5
## Rate of changing weight, i.e. in one step WEIGHT is multiplied by either 0.95 or 1.05
PERCENT=5


##get electron number from last line of conv.dat
def get_enumber():
	a=np.loadtxt('conv.dat')
	return a[-1,1]

##bash command to change weight in params.h
def command(newweight,oldweight):
	return r"sed -i 's/constexpr double      WEIGHT_COM.*=.*{}; /constexpr double      WEIGHT_COM= {}; /g' params.h".format(int(oldweight),int(newweight))

##get current weight from params.h
def get_weight():
	x=0
	with open('params.h') as f:
		line=f.readline()
		while line:
			line=f.readline()
			if "WEIGHT_COM" in line and x==0:
				x+=1
				s=line.split('=')[1][:-1]
				g=s.split(';')[0][:]
				print("before splitting is {}".format(g))
				os.system("sed -i 's/constexpr double      WEIGHT_COM.*=.*{}; /constexpr double      WEIGHT_COM={}; /g' params.h".format(g,int(float(g))))
				weight=float(g)	
	print(weight)
	return weight

##decide whether weight should be smaller or larger, then change weight by PERCENT, then wait 1 hour (can be changed)
def do_iteration(weight, enumber):
	new_weight=weight/(ENDNUM/enumber)
	ratio=enumber/ENDNUM	

	print(ratio)
	intermed_weight=weight
	oldweight=weight
	while(new_weight<intermed_weight):	
		if ratio>(1-PERCENT/100):	
			intermed_weight*=ratio
		else:
			intermed_weight*=(1-PERCENT/100)
		os.system(command(intermed_weight,oldweight))
		os.system("bash comp_PICit")
		time.sleep(3600)
		oldweight=intermed_weight
		ratio=new_weight/intermed_weight

	while(new_weight>intermed_weight):	
		if ratio<(1+PERCENT/100):	
			intermed_weight*=ratio
		else:
			intermed_weight*=(1+PERCENT/100)
		os.system(command(intermed_weight,oldweight))
		os.system("bash comp_PICit")
		time.sleep(3600)
		oldweight=intermed_weight
		ratio=new_weight/intermed_weight


## Loop until electron number is close to ENDNUM, wait between iterations if necessare (time.sleep)
s=get_enumber()
while(np.abs(s-ENDNUM)/ENDNUM>0.02):
	do_iteration(get_weight(),get_enumber())
	time.sleep(1500)
	s=get_enumber()

















	
