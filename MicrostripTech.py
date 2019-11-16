# Functions to calculate the fundamental characteristics of microstrip technology based on certain 
# classical approximations.

import numpy as np
import math

#################################################################################################
# Function to calculate the effective dielectric constant, only depends of
# the geometric values and dielectric constant of the substrate
# The approximation is due to Wheeler and Schneider, Book: Microstrip Lines
# and Slotlines, Gupta, 2nd Edition, 1996
def Effective(Er, W, H):
    # Condition for W/H ratio given by function F
    F = 0.0
    if((W/H) <= 1.0):
        F = (1.0/np.sqrt(1.0+(12.0*H/W)))+(0.04*((1-(W/H))**2))
    else:
        F = (1.0/np.sqrt(1.0+(12.0*H/W)))
    Eeff = ((Er+1.0)/2.0)+(((Er-1.0)/2.0)*F)
    # Return value of effective constant
    return Eeff


# Function to calculate the impendance of the line, which only depends of
# the geometric values and Effective/realtive dielectric constant
# The approximation was founded in Pozar, Microwave Engineering (3rd Edition)
# Additionally the coorection of finite strip thickness
def Impedance(Eff, W, H):
    # Condition for (W/H ratio)
    if((W/H) <= 1):
        Z = (60.0/np.sqrt(Eff))*np.log((8.0*H/W)+(W/(4.0*H)))
    else:
        Z = (120.0*np.pi)/((np.sqrt(Eff))*((W/H)+1.393+(0.667*np.log((W/H)+1.444))))
    # return the value of impedance calculated
    return Z


# Function to calculate the capacitance of Microstrip transmision line, the
# final value is in Farads, the value has units of output F/m, capacitance
# per length.
# The approximation is based in M. V. Schneider, 1969, "Microstrip Lines
# for Microwave Integrated Circuits"
def Capacitance(Er, W, H):
    c = 299792548                                                               # Velocity of light in vacuum
    e0 = (8.8541878176e-12)                                                     # Farads/inches
    # It is another value for effective dielectric constant due to
    # capacitance calculus and approximation
    Eff = Effective(Er, W, H)
    # Ratio W/H classification
    if ((W/H) <= 1):
        CL = Eff/(60.0*c*np.log(((8.0*H)/(W)) + ((W)/(4.0*H))))
    else:
        CL = (Eff*e0)*((W/H)+2.42-(0.44*H/W)+np.power((1.0-(H/W)),6.0))
    # return the capacitance per length value
    return CL


# Function to calculate the inductance of Microstrip transmision line, the
# final value is in Farads, the value has units of output F/m, capacitance
# per length.
# The approximation is based in M. V. Schneider, 1969, "Microstrip Lines
# for Microwave Integrated Circuits"
def Inductance(W, H):
    c = 299792548;                                                          # Vacuum light velocity
    Iv0 = 120*np.pi;                                                        # Impedance of vacuum
    if((W/H) <= 1):
        Ln = (Iv0/(2*np.pi*c))*(np.log((8.0*H/W)+(W/(4*H))))
    else:
        Ln = (Iv0/c)*(np.power(((W/H)+1.393+(0.667*(np.log((W/H)+1.444)))),-1.0))
    # return the inductance per length of the transmission line
    return Ln


# Function to calculate the frequency depend line resistance using the
# classical approximation of R. A. Pucel (1968) "Losses in Microstrip"
# The function calculate the resistance per length (Ohm/m) and to calculate
# the losses for conducting effect only required the impedance value.
# R = sqrt((Rdc^2)+(Rac^2))
# And each value of resistance corresponds to 
# Rdc = RdcS + RdsG (Strip + Ground)
# Rac = RacS + RacG (Strip + Ground)
# The 8.68 factor is the case pro dB transformation presented in the article, so
# this case only evaluated the losses in Np per length
def Resistance(W, H, T, Freq, mu, sigma):
    RD = 1.0/(np.sqrt(np.pi*Freq*mu*sigma))                                        # Skin depth (assuming same material for strip and ground)
    Rdc = 2.0/(W*T*sigma)                                                          # DC line resistance, assuming equal W reflects in Ground part
    W1 = W + ((T/np.pi)*(np.log(1.0+((2.0*H)/(T)))))                               # Width modified for the calculation of resistance
    # Ratio classification using W/H cocient
    if((W/H) <= (1.0/(2.0*pi))):
        RacS = (RD/(2*np.pi*H))*(1.0-((W1/(4.0*H))**2.0)*(1.0+(((2.0*H)/W1)*(1.0+(T/(np.pi*W))+((1.0/np.pi)*(log((4.0*np.pi*W)/(T)))))))
        RacG = (RD/(2.0*np.pi*H))*(1.0-(((W1)/(4.0*H))**2.0))
    elif (((W/H) > (1.0/(2.0*np.pi))) and ((W/H) <= 2.0)):
        RacS = (RD/(2.0*np.pi*H))*(1.0-(((W1)/(4.0*H))**2.0))*(1.0-(T/(np.pi*W1))+(((2.0*H)/W1)*(1.0+((1.0/np.pi)*(np.log((4.0*H)/(T)))))))
        RacG = (RD/(2.0*np.pi*H))*(1.0-(((W1)/(4.0*H))**2.0))*(1.0-(T/(np.pi*W1)))
    else:
        RacSaux =(RD/H)*(((W1/H)+((2.0/np.pi)*np.log(2.0*np.pi*((W1/(2.0*H))+(0.94*np.exp(1))))))**(-2.0))*((W1/H)+((W1/(np.pi*H))/(((W1)/(2.0*H))+0.94)))
        RacS = RacSaux*(1.0-((T)/(np.pi*W1))+((2.0*H)/(W1))+((2.0*H)/(np.pi*W1*np.log((2.0*H)/(T)))))
        RacGaux =(RD/(np.pi*H))*(((W1/H)+((2.0/np.pi)*np.log(2.0*np.pi*((W1/(2.0*H))+(0.94*np.exp(1))))))**(-2.0))*((W1/H)+((W1/(np.pi*H))/(((W1)/(2.0*H))+0.94)))
        RacG = RacGaux*(1.0-((T)/(np.pi*W1)))
    # Final value of resistance per length
    Rac = RacS + RacG
    R = np.sqrt((Rac**2.0) + (Rdc**2.0))
    return R


# Function to total conductor losses, the final value has units of dB/m
# The approximation is presented in Pozar Book Microwave Engineering
# The final element correspond to rugosity profile additional model.
# 0.4um for Roger 3540B material rugosity profile
def ConductorLoss2(W, Freq, Z0, mu, sigma):
    RD = np.sqrt(np.pi*Freq*mu/sigma);
    # Losses in Np/m
    aCT = RD/(Z0*W);
    # Losses in dB/m
    aCT = 8.686*aCT;
    aCT = aCT*(1+((2/np.pi)*np.arctan(1.4*(0.4e-6)/RD)));
    # return the total aproximation 
    return aCT


# Function to calculate the dielectric losses in the material, the final
# value has units of dB/m
# The apprximation is based in Gupta, Hammerstad and Bekkadal. The book
# referenced is Foundations for Microstrip Circuit Design, Edwards 2016
# L is the length of the section analized (different to L0 that is the 
# wavelength in the vacuum)
# TaD is the losses constant for the dielectric
def DielectricLoss(L0, Er, Ef, TaD):
    # return the standard value of losses in Np/m (important for Q value)
    aD = (((2*np.pi/L0)*(Er))/(2*np.sqrt(Ef)))*((Ef-1)/(Er-1))*(TaD);
    # Convert to dB/m
    aD = 8.686*aD;
    # return the total dielectric losses
    return aD


# Function to calculate the two capacitances of the gap between adjacent
# microstrips lines. The Gap has a length of G and the final values of
# capacitance are:
# Cp = Capacitances between the edge and ground (like a open ended microstrip)
# Cg = Capacitance between the two edges.
# Ck = (Cp, Cg) [F/m]
# The model has restrictions for 0.5 <= (W/H) <= 2 however we will use this
# approximation for any near value of limits
# The model is based in R. Gang (1978) "Microstrip Discontinuities" and
# used for Reference book of Gupta
def CapGap(Er, G, W, H):                                                         
    # Classification for odd capacitance of Er = 9.6
    m0 = (W/H)*((0.619*(np.log10(W/H)))-0.3853)
    k0 = 4.26-(1.453*(np.log10(W/H)))
    Codd = W*((G/W)**(m0))*(np.exp(k0))*((Er/9.6)**0.8)
    # Classification for even capacitance of Er = 9.6
    if((G/W) <= 0.3):
        me = 0.8675
        ke = 2.043*((W/H)**(0.12))
    else((G/W) > 0.3):
        me = ((1.565)/((W/H)**(0.16)))-1.0
        ke = 1.97-((0.03)/(W/H))
    Ceven = 12.0*W*((G/W)**(me))*(np.exp(ke))*((Er/9.6)**0.9)
    # Final calculation for Cp and Cg (C1 and C2 in the thesis document)
    Cp = 0.5*Ceven
    Cg = 0.5*(Codd-(0.5*Ceven))
    Cap = np.array([(Cp*(1e-12)), (Cg*(1e-12))])                                                       
    # return an array of the two capacitances, Cp and Cg
    return Cap