# -*- coding: utf-8 -*-
"""
@author: Gavin Standish
"""

import numpy as np
import matplotlib.pyplot as plt

'Defining of Variables'
MuB = 0.99    # Correlation Coefficient, MuB/MuW, described to be equal to 0.99 in pg.887 of pt.I
MuW = 1
pi = np.pi
g = 9.8    # m/s^2 , gravity

'Properties of Fluid'
Cp = 4.186    # J/g*K , specific heat constant (water) 
p = 1000    # kg/m^3 , fluid density (water)
Mu = 1.0016    # kg/m*s , fluid dynamic viscosity (water @ 293K)
n = 0.18    # liquid Heating value (pg. 894 pt.II)
m_dot = 1    # kg/s , mass flow rate (Assumed)
m = 1    # kg , mass of fluid (Assumed)
K = 293    # K, temperature of fluid (Based on Mu value)

'Properties of Twisting Tube'
W = 400    # K , at tube wall temperature (Assumed)
d = 0.05    # m , tube inside diameter (Assumed)
BdelT = 0.02    # , Coeff of Isobaric Thermal Expansion * Change in Temp (pg.891 pt.II for water)
deltad = 0.05    # Dimensionless ratio between thickness of twisted tape and tube inside diameter (pg.890 pt.II)
y = 3    # H/d , twist ratio as suggested by Manglik & Bergles papers pt.I and pt.II
L = 0.217    # m , axial length based on L=sqrt((P^2)+(pi*d)^2) and y=P/d

'Various Number Calculations'
Re = np.linspace(1,50000,1000)    # Assumption that Reynolds Number is Laminar from <2300
ReL = np.linspace(1,2300,1000)    # Assumption that Reynolds Number is Turbulent from <4000
ReT = np.linspace(4000,50000,1000)    # It should be noted that from 2300-4000 the fluid is in a transitional period.
Ra = (g)*(p**2)*(d**3)*(BdelT) / (Mu**2)    # Rayleigh Number
Pr = (Mu*Cp) / (W/((m**2)*K))    # Prandtl Number
# Pr = 7
Gz = (m_dot*Cp) / ((W/((m**2)*K))*L)    # Graetz Number
SwL = (ReL/np.sqrt(y)) * (pi/(pi-(4*deltad))) * (1+(pi/((2*y)**2)))**(1/2)    # Dimensionless Swirl Parameter for Laminar Re values
SwT = (ReT/np.sqrt(y)) * (pi/(pi-(4*deltad))) * (1+(pi/((2*y)**2)))**(1/2)    # Same but for Turbulent Re values

'Nusselt Number Formula'
NuL = 4.612 * ((((((1+0.0951*Gz**0.894)**2.5)+(6.413*10**-9)*((SwL*Pr**0.391))**3.835)**2.0)+(2.132*10**-14)*(ReL*Ra)**2.23)**0.1) * (MuB/MuW)**0.14
# Mean Nusselt Number Equation 17 of Manglik & Bergles pt.I - Laminar Flow
NuT = 0.023 * (ReT**0.8) * (Pr**0.4) * ((pi/(pi-(4*deltad)))**0.8) * (((pi+2-(2*deltad))/(pi-(4*deltad)))**0.2) * (MuB/MuW)**n
# Mean Nusselt Number Equation 9 of Manglik & Bergles pt.II - Transition and Turbulent Flow

fL = (15.767/ReL)*((1+(10**-6)*(SwL**2.55))**(1/6))*(1+((pi/(2*y))**2))*(((pi+2-(2*deltad))/(pi-(4*deltad)))**2)*(pi/(pi-(4*deltad)))
# Mean Nusselt Number Equation 10 of Manglik & Bergles pt.II - Transition and Turbulent Flow
fT = (0.0791/(ReT**0.25))*((pi/(pi-(4*deltad)))**1.75)*(((pi+2-(2*deltad))/(pi-(4*deltad)))**1.25)*(1+(2.752/(y**1.29)))
# Mean Nusselt Number Equation 7 of Manglik & Bergles pt.II - Transition and Turbulent Flow


'Creation of Plots'
plt.figure()
plt.plot(ReL, fL, label='Laminar')
plt.plot(ReT, fT, label='Turbulent')
plt.ylim(0.001, 1)
plt.xlim(1, 100000)
plt.xscale('log')    # Setting plot axes to logarithmic scale. Just used log for this one because it looked similar to results from paper.
plt.yscale('log')
plt.xlabel('Re')
plt.ylabel('f')

plt.legend()
plt.grid()

plt.figure()
plt.plot(ReL, NuL, label='Laminar')
plt.plot(ReT, NuT, label='Turbulent')
# plt.ylim(1, 100)
plt.xlim(1, 50000)
# plt.xscale('log')    # Setting plot axes to logarithmic scale.
# plt.yscale('log')
plt.xlabel('Re')
plt.ylabel('Nu')

plt.legend()
plt.grid()
plt.show()


