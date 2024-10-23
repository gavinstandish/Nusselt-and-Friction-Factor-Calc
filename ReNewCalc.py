# -*- coding: utf-8 -*-
"""
@author: gavin
"""

import numpy as np
import matplotlib.pyplot as plt
import iapws.iapws97
from iapws.iapws97 import IAPWS97

'Liquid Water Properties'
t1 = 293 # K (Standard Temperature)
p1 = 0.1013 # MPa (Standard Pressure)
LiqWater = IAPWS97(T=t1,P=p1,x=0)
Cp = LiqWater.cp # kJ/kg*K , Specific Heat
p = LiqWater.rho # kg/m^3 , Density
Mu = LiqWater.mu # Pa*s , Dynamic Viscosity
k = LiqWater.k # W/m*K , Thermal Conductivity
Pr = LiqWater.Prandt # Prandtl Number
m_dot = 1 # kg/s
m = 1 # kg

'Twisting Tube Properties'
d = 0.05 # m , Tube Diameter
y = 3.0 # H/d , twist ratio as suggested by Manglik & Bergles papers pt.I and pt.II
L = 0.217 # m , axial length based on L=sqrt((P^2)+(pi*d)^2) and y=P/d
BT = 0.02 # Coeff of Isobaric Thermal Expansion * Change in Temp (pg.891 pt.II for water)
deld = 0.05 # Dimensionless ratio between thickness of twisted tape and tube inside diameter (pg.890 pt.II)
n = 0.18 # Liquid Heating value (pg. 894 pt.II)

'Other Properties'
pi = np.pi
MuBMuW = 0.99    # Correlation Coefficient, MuB/MuW, described to be equal to 0.99 in pg.887 of pt.I
g = 9.8    # m/s^2 , gravity constant

'Reynolds Numbers'
Re = np.linspace(1,50000,1000) # These values just function as the x-values
ReL = np.linspace(1,2300,1000) # Laminar Re Range
ReT = np.linspace(10000,50000,1000) # Turbulent Re Range

'Defining Secondary Functions'
def Ra(g,p,d,BT,Mu): # Rayleigh's Number
    Ra = (g)*(p**2)*(d**3)*(BT) / (Mu**2)
    return Ra
def Gz(m_dot,Cp,k,L): # Graetz Number
    Gz = (m_dot*Cp) / (k*L)
    return Gz
def SwL(ReL,y,deld,pi): # Laminar Swirl Parameter
    SwL = (ReL/np.sqrt(y)) * (pi/(pi-(4*deld))) * (1+(pi/((2*y)**2)))**(1/2)
    return SwL
def SwT(ReT,y,deld,pi): # Turbulent Swirl Parameter
    SwT = (ReT/np.sqrt(y)) * (pi/(pi-(4*deld))) * (1+(pi/((2*y)**2)))**(1/2)
    return SwT

'Manglik and Bergles (M&B) Twisting Tape Nusselt Numbers and Friction Factors'
def NuL(Gz,SwL,Pr,ReL,Ra,MuBMuW): # Laminar Nusselt Numbers - Eq. 17 M&B pt.I
    NuL = 4.612 * ((((((1+0.0951*(Gz(m_dot,Cp,k,L))**0.894)**2.5)+(6.413*10**-9)*(((SwL(ReL,y,deld,pi))*Pr**0.391))**3.835)**2.0)+(2.132*10**-14)*(ReL*(Ra(g,p,d,BT,Mu)))**2.23)**0.1) * (MuBMuW)**0.14
    return NuL
def NuT(ReT,Pr,deld,MuBMuW,n,pi): # Turbulent Nusselt Numbers - Eq. 9 M&B pt.II
    NuT = 0.023 * (ReT**0.8) * (Pr**0.4) * ((pi/(pi-(4*deld)))**0.8) * (((pi+2-(2*deld))/(pi-(4*deld)))**0.2) * (MuBMuW)**n
    return NuT
def fL(ReL,SwL,pi,y,deld): # Laminar Friction Factor - Eq. 10 M&B pt.II
    fL = (15.767/ReL)*((1+(10**(-6))*((SwL(ReL,y,deld,pi))**2.55))**(1/6))*(1+((pi/(2*y))**2))*(((pi+2-(2*deld))/(pi-(4*deld)))**2)*(pi/(pi-(4*deld)))
    return fL
def fT(ReT,pi,deld,y): # Turbulent Friction Factor - Eq. 7 M&B pt.II
    fT = (0.0791/(ReT**0.25))*((pi/(pi-(4*deld)))**1.75)*(((pi+2-(2*deld))/(pi-(4*deld)))**1.25)*(1+(2.752/(y**1.29)))
    return fT

'Dittus-Boelter Nusselt Number Correlation'
# Only valid when 0.6<Pr<160 and Re>10000
def NuDh(ReT,Pr):
    NuDh = 0.023 * (ReT**0.8) * (Pr**0.4)
    return NuDh

'Blasius Friction-Factor Correlation'
# Only valid for Re>=TransitionalFlow, using Re>=10,000 for comparison
def fBla(ReT):
    fBla = 0.3164 * (ReT**(-1/4))
    return fBla

'Creation of Plots'
# M&B Twisting Tube Friction-Factor Laminar vs. Turbulent including Blasisus correlation.
plt.figure() 
plt.plot(ReL, fL(ReL,SwL,pi,y,deld), label='Laminar')
plt.plot(ReT, fT(ReT,pi,deld,y), label='Turbulent')
plt.plot(ReT, fBla(ReT), label='Blasius Comparison')
plt.ylim(0.001, 100)
plt.xlim(1, 50000)
plt.xscale('log')    # Setting plot axes to logarithmic scale. Just used log for this one because it looked similar to results from paper.
plt.yscale('log')
plt.xlabel('Re')
plt.ylabel('f')

plt.legend()
plt.grid()

# Comparison of Turbulent M&B Friction-Factor to Blasius Straight Tube comparison.
plt.figure()
plt.plot(ReT, fT(ReT,pi,deld,y), label='Turbulent')
plt.plot(ReT, fBla(ReT), label='Blasius Comparison')
# plt.ylim(0.02, 0.04)
# plt.xlim(10000, 50000)
# plt.xscale('log')    # Setting plot axes to logarithmic scale. Just used log for this one because it looked similar to results from paper.
# plt.yscale('log')
plt.xlabel('Re')
plt.ylabel('f')

plt.legend()
plt.grid()

# M&B Twisting Tube Nusselt Number Laminar vs. Turbulent including Dittus-Boelter correlation.
plt.figure()
plt.plot(ReL, NuL(Gz,SwL,Pr,ReL,Ra,MuBMuW), label='Laminar')
plt.plot(ReT, NuT(ReT,Pr,deld,MuBMuW,n,pi), label='Turbulent')
plt.plot(ReT, NuDh(ReT,Pr), label='Dittus-Boelter Comparison')
# plt.ylim(1, 100)
plt.xlim(1, 50000)
# plt.xscale('log')    # Setting plot axes to logarithmic scale.
# plt.yscale('log')
plt.xlabel('Re')
plt.ylabel('Nu')

plt.legend()
plt.grid()
plt.show()