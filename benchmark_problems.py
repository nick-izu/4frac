#Sandia National Laboratories is a multimission laboratory managed and operated 
#by National Technology & Engineering Solutions of Sandia, LLC, a wholly owned 
#subsidiary of Honeywell International Inc., for the U.S. Department of Energy’s 
#National Nuclear Security Administration under contract DE-NA0003525.

###########################
##DECOVALEX Benchmark Cases

# Based on:
# Kolditz, et al. (2015) Thermo-Hydro-Mechanical-Chemical Processes in
# Fractured Porous Media: Modelling and Benchmarking, Closed Form Solutions,
# Springer International Publishing, Switzerland.

# Author: Rosie Leone
# SAND2020-7449 O

from numpy import *
from scipy import special
import matplotlib.pyplot as plt
import pandas as pd

filename1 = "Flow_Benchmarks.xlsx"
filename2 = "Transport_Benchmarks.xlsx"

#1D steady flow boundary conditions of first kind
def steady_1d_flow():

    Lx = 100 #model domain [m]
    nx = 10 #num of grid cells
    dx = Lx/nx #length of grid cells [m]
    P0 = 2 #Pressure at x=0 [MPa]
    
    x_soln = linspace(0.+(dx/2.),Lx-(dx/2.),nx) 
    p_soln = zeros(nx+2)
    
    p_soln = array((-x_soln/Lx) + P0)*1e6 #[Pa]
    soln = column_stack((x_soln,p_soln))
    
    return soln

#2D steady flow boundary conditions of first kind
def steady_2d_flow():
    Lx = 1 #Model domain in x direction [m]
    nx = 20 #num of grid cells in x direction
    Ly = 1 #Model domain in y direction [m]
    ny = 20 #num of grid cells in y direction
    dx = Lx/nx #length of grid cells in x direction [m]
    dy = Ly/ny #length of grid cells in y direction [m]
    p0 = 2e6 #[Pa]
    p_offset = 1e6 #[Pa]
    
    p_soln = zeros((nx,ny))                      
    
    x_soln = linspace(0.+(dx/2.),Lx-(dx/2.),nx)  # [m]
    y_soln = linspace(0.+(dy/2.),Ly-(dy/2.),ny)  # [m]
    
    # calc the analytical solution
    i = -1
    for x in x_soln:
      i = i + 1
      j = -1
      for y in y_soln:
        j = j + 1
        p_soln[i,j] = (p0-p_offset)*(x/Lx)*(y/Ly) + p_offset #[Pa]
        
    return x_soln, y_soln, p_soln

#3D steady flow boundary conditions of first kind
def steady_3d_flow():
    Lx = 1 #Model domain in x direction [m]
    Ly = 1 #Model domain in y direction [m]
    Lz = 1 #Model domain in z direction [m]
    nx = 10 #num of grid cells in x direction
    ny = 10 #num of grid cells in y direction
    nz = 10 #num of grid cells in z direction
    dx = Lx/nx #length of grid cells in x direction [m]
    dy = Ly/ny #length of grid cells in y direction [m]
    dz = Lz/nz #length of grid cells in z direction [m]
    
    p0 = 1e6 #[Pa]
    
    p_soln = zeros((nx,ny,nz))                   
    
    x_soln = linspace(0.+(dx/2.),Lx-(dx/2.),nx)  # [m]
    y_soln = linspace(0.+(dy/2.),Ly-(dx/2.),ny)  # [m]
    z_soln = linspace(0.+(dz/2.),Lz-(dx/2.),nz)  # [m]
    
    i = -1
    for x in x_soln:
      i = i + 1
      j = -1
      for y in y_soln:
          j = j + 1
          k = -1
          for z in z_soln:
              k = k + 1
              p_soln[i,j,k] = p0*((x/Lx)+(y/Ly)+(z/Lz)) #[Pa]
              
    return x_soln,y_soln,z_soln,p_soln


#Transient-Conservative Tracer
def F_inert(x,t,c0,phi,D,q,R,T):

    return c0*0.5*((exp((q*x)/(phi*D))*special.erfc((x+(q/phi)*t)/(2*(D*t)**0.5)))+special.erfc((x-(q/phi)*t)/(2*(D*t)**0.5)))

#Transient-Adsorbing Tracer
def F_sorb(x,t,c0,phi,D,q,R,T):

   return c0*0.5*((exp((q*x)/(phi*D))*special.erfc((R*x+(q/phi)*t)/(2*(D*R*t)**0.5)))+special.erfc((R*x-(q/phi)*t)/(2*(R*D*t)**0.5)))

#Transient-Decaying Tracer
def F_decay(x,t,c0,phi,D,q,R,T):

   A = q/(2*phi*D)
   B = log(2)/(D*T) +A**2
   
   return c0*0.5*((exp(x*(A+B**0.5))*special.erfc((x+2*t*(B*D**2)**0.5)/(2*(D*t)**0.5)))+(exp(x*(A-B**0.5))*special.erfc((x-2*t*(B*D**2)**0.5)/(2*(D*t)**0.5))))

#Transient Transport    
def transient_transport_1d(tracer):
    
    L= 10 #Model domain in x direction [x]
    k = 1e-11 #permeability [m^2]
    phi = 0.4 #porosity
    mu = 1e-3 #viscosity [mPa *s]
    D = 1e-4 #Diffusion coefficient [m^2/s]
    Kd = 6.8e-4 #Distribution coefficient [m^3/kg]
    ps = 2000. #Solid grain density [kg/m^3]
    p0 = 1e5 #Pressure at inlets [Pa]
    R = 1+(1-phi)/phi*Kd*ps #Retardation Factor
    T=log(2)/0.5e-4 #half life [s]
   
    c0 = 1.0 #Initial concentration [mol/L]
    q = k/mu * (p0/L) #Specific discharge [m/s]
    
    if tracer == "inert":
        F = F_inert
        
    elif tracer == "sorb":
        F = F_sorb
        
    elif tracer == "decay":
        F = F_decay
  

    nx = 200 #num of grid cells in x direction
    t0 = 15000 #time when concentration goes to zero at x=0

    dx = L/nx #discretization [m]

    x_soln = linspace(0+(dx/2),L-(dx/2),nx) #[m]

    t_soln = ([20000]) #[s]
    
    soln = zeros(len(x_soln))
       
    for t in t_soln:

        for i in range(len(x_soln)):
            if t>t0:
              soln[i] = F(x_soln[i],t,c0,phi,D,q,R,T)-F(x_soln[i],t-t0,c0,phi,D,q,R,T)
    
            else:
              soln[i] = F(x_soln[i],t,c0,phi,D,q,R,T)
              
    soln = column_stack((x_soln,soln))         
    return soln

###############################################################################
######Benchmark Steady Flow 1D
soln = steady_1d_flow()
##plot
plt.figure()
plt.plot(soln[:,0],soln[:,1])
plt.title('1D Flow')

#Write to Excel
df1 = pd.DataFrame(soln)
writer1 = pd.ExcelWriter(filename1,engine='xlsxwriter')
df1.to_excel(writer1,index=False,header=['X Cord Cell Center [m]','Analytical Pressure Solution [Pa]'],sheet_name='1D Flow')


#######Benchmark Steady Flow 2D
x,y,p_soln = steady_2d_flow()


##plot 
plt.figure()
X,Y = meshgrid(x,y) 
levels = linspace(amin(p_soln),amax(p_soln),11)   
plt.contourf(X,Y,p_soln,levels)
plt.colorbar()
plt.title('2D Flow')

#Get solution at y = 0.475 m
soln = column_stack((x,p_soln[:,9])) 

#Write to Excel
df2 = pd.DataFrame(soln)
df2.to_excel(writer1,index=False,header=['X Cord Cell Center [m]','Analytical Pressure Solution [Pa] at y = 0.475 m'],sheet_name='2D Flow')


######Benchmark Steady Flow 3D
x,y,z,p_soln = steady_3d_flow()

##plot at z = 0.45
plt.figure()
X,Y = meshgrid(x,y) 
levels = linspace(amin(p_soln[:,:,4]),amax(p_soln[:,:,4]),11)   
plt.contourf(X,Y,p_soln[:,:,4],levels)
plt.colorbar()
plt.title('3D Flow at z=0.45 m')

##write to excel
soln = column_stack((x,p_soln[:,4,4]))
df3 = pd.DataFrame(soln) 
df3.to_excel(writer1,index=False,header=['X Cord Cell Center [m]','Analytical Pressure Solution [Pa] at y=0.45 m, z=0.45'],sheet_name='3D Flow')

#Save Flow excel spreadsheet
writer1.save()   


#####Benchmark Transient Transport
#Conservative Tracer
soln = transient_transport_1d("inert")
##plot
plt.figure()
plt.plot(soln[:,0],soln[:,1],label='Conservative Tracer')

#save to excel       
df4 = pd.DataFrame(soln)
writer2 = pd.ExcelWriter(filename2,engine='xlsxwriter')
df4.to_excel(writer2,index=False,header=['X Cord Cell Center [m]','Analytical Concentration Solution [M]'],sheet_name='Conservative')

#Decaying Tracer
soln = transient_transport_1d("decay")
#plot
plt.plot(soln[:,0],soln[:,1],label='Decaying Tracer')

#save to excel  
df5 = pd.DataFrame(soln)
df5.to_excel(writer2,index=False,header=['X Cord Cell Center [m]','Analytical Concentration Solution [M]'],sheet_name='Decay')

#Adsorbing Tracer
soln = transient_transport_1d("sorb")
#plot
plt.plot(soln[:,0],soln[:,1],label='Adsorbing Tracer')
plt.legend()
plt.title('Transient Transport')

#save to excel  
df6 = pd.DataFrame(soln)
df6.to_excel(writer2,index=False,header=['X Cord Cell Center [m]','Analytical Concentration Solution [M]'],sheet_name='Sorb')

#Save Transport spread sheet
writer2.save()




