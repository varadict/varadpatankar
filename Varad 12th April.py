# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 21:54:55 2016

@author: hp

"""
from scipy import exp,linspace
from scipy.integrate import odeint
from pylab import plot,show
from scipy.optimize import fsolve
from matplotlib import pyplot as plt    


Q=2.0           #Flow Rate of the Reactant
V=1.0           #Volume of the Reactor
U=110.0         #Overall Heat Tr Coeff for the HX
a=150.0         #HT Area
R=8.314         
T0=1035.0       #Inlet Temperature
Ta=1150.0       #Temperature of the surrounding
P0=162000.0     #Inlet Pressure
X0=0            #Conversion at inlet
Ca0=18.8        #COncentration at inlet
Fa0=Q*Ca0
n=25          #no of CSTRs in Series
'Increase n to get a smoother curve.After a point we get the same curve of CSTRS and the PFR'


def ideal(x,v):                          #For an ideal PFR
    X=x[0]
    T=x[1]
    k=3.58*exp(34222.0*(1.0/T0-1.0/T))   #Reaction Constant
    Cpa=26.63+0.183*T-45.86*10**-6*T**2  #Cp values for ketone
    Cpb=20.04+0.0954*T-30.95*10**-6*T**2 #Cp values for Ketene
    Cpc=13.39+0.077*T-18.71*10**-6*T**2  #Cp values for Methane
    dHr=80770.0+6.8*(T-298.0)-5.75*10**-3*(T**2-298.0**2)-1.27*10**-6*(T**3-298.0**3) #Reaction Heat
    ra=k*Ca0*(1.0-X)/(1.0+X)*T0/T
    e1=ra/Fa0
    e2=(U*a*(Ta-T)-ra*dHr)/(Fa0*(Cpa+X*(Cpb+Cpc-Cpa)))
    return [e1,e2]
    
def nonideal(x):                            #For a non ideal PFR i.e Number of CSTRS in series
    X=x[0]
    T=x[1]
    k=3.58*exp(34222.0*(1.0/T0-1.0/T))      #Reaction Constant
    Cpa=26.63+0.183*T-45.86*10**-6*T**2     #Cp values for Ketone
    Cpb=20.04+0.0954*T-30.95*10**-6*T**2    #Cp values for Ketene
    Cpc=13.39+0.077*T-18.71*10**-6*T**2     #Cp values for Methane
    dHr=80770.0+6.8*(T-298.0)-5.75*10**-3*(T**2-298.0**2)-1.27*10**-6*(T**3-298.0**3) #Reaction Heat
    ra=k*Ca0*(1.0-X)/(1.0+X)*T0/T
    e1=(X-X0)/V*n-ra/Fa0
    e2=(T-T0)/V*n-(U*a*(Ta-T)-ra*dHr)/(Fa0*(Cpa+X*(Cpb+Cpc-Cpa)))
    return [e1,e2]
    
X=[X0]                                  #Initial Conversion
T=[T0]                                  #Initial Temperature
vv=[0]                                  
v=linspace(0,V,101)                     #Dividing the entire length into infinitismal units
x=odeint(ideal,[X0,T0],v)               #Solving  for n CSTRS in series 
for i in range(n):
    X.append(X0)
    T.append(T0)
    vv.append((i+1)*V/float(n))         
    "The volume of one  cstr is taken by dividing the volume of pfr by no. of CSTSRs"

    
    [X0,T0]=fsolve(nonideal,[0.3,1150])
    X.append(X0)
    T.append(T0)
    vv.append((i+1)*V/float(n))
    
print "Conversion and Final Temperature for a PFR is :", X[-1],T[-1]
print "Conversion and Final Temperature for n CSTRS is :", x[-1,0],x[-1,1]   

plt.plot(vv,X)
plt.title('Conversion Profile')
plt.ylabel('Conversion')
plt.xlabel('CSTR No.')
plt.show()

plt.plot(vv,T)
plt.title('Temperature Profile')
plt.ylabel('Temperature')
plt.xlabel('CSTR No.')
plt.show()

plt.plot(v,x[:,0])
plt.title('Conversion Profile')
plt.ylabel('Conversion')
plt.xlabel('Length')
plt.show()

plt.plot(v,x[:,1])
plt.title('Temperature Profile')
plt.ylabel('Temperature')
plt.xlabel('CSTR No.')
plt.show()

