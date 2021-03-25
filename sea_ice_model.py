import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


#definition des constantes :

#Physical constants
Lfus = 3.35e5        #Latent heat of fusion for water [J/kg]
rhoi = 917           #Sea ice density [kg/m3]
ki = 2.2             #Sea ice thermal conductivity [W/m/K]
ks = 0.31            #Snow thermal conductivity [W/m/K]
sec_per_day = 86400  #Second in one day [s/day]

#Bottom boundary condition
T_bo = -1.8          #Bottom temperature [C]



#conditions initiales
thick_0 = 0.1      #épaisseur initiale [m]
T_0 = -10           #température initiale [C]
Q_w = 2             #flux de chaleur océanique [W/m2]
snow_0 = 0.15       #initiale snow thickness [m]



#Evolution de l'épaisseur de glace-------------------------------------------------------------------
thick = np.zeros(30)  #tableau de l'épaisseur de glace
thick[0] = thick_0    #set de l'épaisseur initiale dans notre tableau    

Q_c = np.zeros(30)    #tableau du flux de chaleur

def get_thick(T_0,T_bo, ki, rhoi, Lfus):
    
    for j in range(0,29):
        Q_c[j] = ki * (T_0-T_bo)/thick[j]                         #Calcul du flux de chaleur
        thick[j+1] =  thick[j] - Q_c[j]*sec_per_day/(rhoi*Lfus)   #Calcul de l'épaisseur de glace
        
    return[thick]   

#Avec le flux de chaleur océanique-------------------------------------------------------------------
thick_w = np.zeros(30)  #tableau de l'épaisseur de glace
thick_w[0] = thick_0    #set de l'épaisseur initiale dans notre tableau

Q_cw = np.zeros(30)     #tableau du flux de chaleur 

def get_thick_w(T_0,T_bo, ki, rhoi, Lfus, Q_w):

    for j in range(0,29):
        Q_cw[j] = ki * (T_0-T_bo)/thick_w[j]                                    #Calcul du flux de chaleur
        thick_w[j+1] =  thick_w[j] - (Q_cw[j]+Q_w)*sec_per_day/(rhoi*Lfus)      #Calcul de l'épaisseur de glace avec le Q_w
        
    return[thick_w]   

#Avec de la neige------------------------------------------------------------------------------------
thick_s = np.zeros(30)  #tableau de l'épaisseur de glace
thick_s[0] = thick_0

Q_cs = np.zeros(30)     #tableau du flux de chaleur


def get_thick_snow(T_0,T_bo, ki, ks, rhoi, Lfus, snow_0, Q_w):

    for j in range(0,29):
        Q_cs[j] = (ki*ks*(thick_s[j]+snow_0))/(ki*snow_0+ks*thick_s[j]) * (T_0-T_bo)/(thick_s[j]+snow_0)    #Calcul du flux de chaleur avec neige 
        thick_s[j+1] =  thick_s[j] - (Q_cs[j] + Q_w)*sec_per_day/(rhoi*Lfus)                                        #Calcul de l'épaisseur de glace avec Q_cs
        
    return[thick_s]


#Surface heat fluxes-----------------------------------------------------------------------------------
#Other physical constants
epsilon = 0.99          #surface emissivity
sigma = 5.67e-8         #Stefan-Boltzmann constantes
Kelvin = 273.15         #Conversion from Celsius to Kelvin
alb = 0.8               #surface albedo


Q_SOL = np.zeros(365)
Q_NSOL = np.zeros(365) 

def solar_flux(day):
    Q_sol= 314 * np.exp((-(day-164)**2)/4608)
    return(Q_sol)

def non_solar_flux(day):
    Q_nsol = 118 * np.exp((-0.5*(day-206)**2)/(53**2))+179
    return(Q_nsol)
    

for j in range(1,366):
    Q_sol = solar_flux(j-1)
    Q_SOL[j-1] = Q_sol
    
for j in range(1,366):
    Q_nsol = non_solar_flux(j-1)
    Q_NSOL[j-1] = Q_nsol
    
#Surface temperature----------------------------------------------------------------------------------
T_0K = T_0 + Kelvin        #convertion en Kelvin
T_boK = T_bo + Kelvin      #convertion en Kelvin
h = 0.5
x0 = 200.15

def f(T_su,h,ki,epsilon,sigma,alb,day):
    return (-(h*epsilon*sigma)/ki*T_su**4 - T_su + h/ki*((solar_flux(day))*(1-alb)+non_solar_flux(day))+T_boK)
 
def df(T_su,h,ki,epsilon,sigma,alb,day):
    return (-(4*h*epsilon*sigma)/ki*T_su**3-1)


def get_root(day,h,ki,epsilon,sigma,alb):
    root = optimize.newton(f,x0,fprime = df, args= (h,ki,epsilon,sigma,alb,day))
    return(root)
    
ROOT = np.zeros(365)

for j in range(365):
    ROOT[j] = get_root(j,h,ki,epsilon,sigma,alb)

       
#Couple temperature and thickness---------------------------------------------------------------------
thick_temp = np.zeros(365) 
thick_temp[0] = thick_0
Q_csurf = np.zeros(365)
ROOT_surf = np.zeros(365)

def get_thick_temp(T_boK, ki, ks, rhoi, Lfus, snow_0, Q_w, epsilon,sigma, alb):

    for j in range(1,365):
        #func = f(T_su,thick_temp[j],ki,epsilon,sigma,alb,j)
        #dfunc = df(T_su,thick_temp[j],ki,epsilon,sigma,alb,j)
        ROOT_surf[j-1] = get_root(j,thick_temp[j-1],ki,epsilon,sigma,alb)
        if ROOT_surf[j-1] > 273.15:
            ROOT_surf[j-1] = 273.15
        Q_csurf[j-1] = (ki*ks*(thick_temp[j-1]+snow_0))/(ki*thick_temp[j-1]+ks*snow_0) * (ROOT_surf[j-1]-T_boK)/(thick_temp[j-1]+snow_0)    #Calcul du flux de chaleur avec neige
        thick_temp[j] =  thick_temp[j-1] - (Q_csurf[j-1]+Q_w)*sec_per_day/(rhoi*Lfus)                                        #Calcul de l'épaisseur de glace avec Q_cs
        
    return(thick_temp, ROOT_surf) 





#Sortie des tableaux et figures-----------------------------------------------------------------------
day30 = np.arange(int(30)) #vecteur temps (30 jours)
[thick] = get_thick(T_0,T_bo, ki, rhoi, Lfus)                       #vecteur épaisseur de glace 
[thick_w] = get_thick_w(T_0,T_bo, ki, rhoi, Lfus, Q_w)              #vecteur épaisseur de glace avec flux de chaleur océanique
[thick_s] = get_thick_snow(T_0,T_bo, ki, ks, rhoi, Lfus, snow_0, Q_w)    #vecteur épaisseur de glace avec neige
thick_temp = get_thick_temp(T_boK, ki, ks, rhoi, Lfus, snow_0, Q_w, epsilon,sigma, alb)[0]  #vecteur épaisseur de la glace en prenant en compte la variation de température en surface

ROOT_surf = get_thick_temp(T_boK, ki, ks, rhoi, Lfus, snow_0, Q_w, epsilon,sigma, alb)[1]

ROOT_surf[-1] = get_root(365,thick_temp[-1], ki, epsilon,sigma,alb)   
Q_csurf[-1] = (ki*ks*(thick_temp[-1]+snow_0))/(ki*snow_0+ks*thick_temp[-1]) * (ROOT_surf[-1]-T_boK)/(thick_temp[-1]+snow_0)



year = np.arange(1,366)


fig1 = plt.figure()
plt.plot(day30, thick)



fig2 = plt.figure()
plt.plot(day30, thick_w)



fig3 = plt.figure()
plt.plot(day30, thick_s)
plt.close()


fig6 = plt.figure()
plt.plot(year, ROOT_surf)
plt.show()
plt.close()



















