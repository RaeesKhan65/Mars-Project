#Packages used
import numpy as np
import matplotlib.pyplot as plt
from numpy import cos,sin,arccos,sqrt
from scipy.optimize import fsolve
from datetime import date


file = open("jup.txt", "r")
data = file.readlines() #Reading in the file
file.close()

ra_obs = [] #List to store all observed Right Ascension Values
de_obs = [] #List to store all observed Declination Values
dates = [] #List to store dates at which values were observed

for val in data:
    
    time = (val[1:13]).replace(" ","") #Extracting the date
    
    ra_h = float(val[23:25]) #Extracting the RA by hours,minutes,seconds
    ra_m = float(val[26:28])
    ra_s = float(val[29:34])
    
    de_sign = float(str(val[35])+str(1))
    de_de = float(val[36:38]) #Extracting the DE by degrees,minutes,seconds
    de_m = float(val[39:41])
    de_s = float(val[42:46])
    
    phi = (ra_h + ra_m/60 + ra_s/3600)*15 #Converting RA to azimuthal angle between 0 and 360

    
    if(de_sign != -1):
        theta = 90.0 - (de_de + (60*de_m+de_s)/3600) #Converting DE to polar angle between 0 and 180
    else:
        theta = 90.0 - (-de_de - (60*de_m+de_s)/3600)
    
    ra_obs.append(phi) #Adding values to my list for each date
    de_obs.append(theta)
    dates.append(time)



def deg(rad): #Function to convert radians to degrees
    return rad*180/np.pi

def rad(deg):#Function to convert degrees to radians
    return deg*np.pi/180

def sub_vec(r1,r2): #Function to subtract two vectors
    x = r1[0]-r2[0]
    y = r1[1]-r2[1]
    z = r1[2]-r2[2]
    
    r=np.array([x,y,z])
    
    return r


def day(val): #Function to return number of days 'd' since Epoch
    Y,M,D = val
    d0 = date(2000, 1, 1)
    d1 = date(Y,M,D)
    delta = d1 - d0
    return (delta.days)


def mean_anomaly_mars(d): #Function to compute mean anomaly of Mars
    return rad((19.41248 + (360*0.0843170)*d)%360)

def mean_anomaly_jupiter(d): #Function to compute mean anomaly of Mars
    return rad(((34.40438-14.75385) + (1/12)*d)%360)

def mean_anomaly_earth(d): #Function to compute mean anomaly of Earth
    return rad((-2.48284 + 0.9856002585*d)%360)

def E_A(data): #Function to calculate Eccentric anomaly
    root = fsolve(func, 1,args = data)
    return root[0]

def func(x,*data): #Function used in calculation of Eccentric anomaly
    e,M = data
    return (x-e*np.sin(x)-M) 

def position(a,e,E): #Function to calculate position vector in orbit
    x = (a)*(cos(E)-e)
    y = (a)*(sin(E)) * (sqrt(1 - e**2))
    z = 0
    
    r = np.array([x,y,z])
    return r


def Rx(angle): #Rotation matrix around x
    return np.matrix([[ 1, 0           , 0           ],
                   [ 0, np.cos(angle),-np.sin(angle)],
                   [ 0, np.sin(angle), np.cos(angle)]])

  
def Rz(angle): #Rotation matrix around z
    return np.matrix([[ np.cos(angle), -np.sin(angle), 0 ],
                   [ np.sin(angle), np.cos(angle) , 0 ],
                   [ 0           , 0            , 1 ]])


def euler_rotation(phi,theta,psi): #Euler rotation matrix 
    
    R = np.matmul(np.matmul(Rz(phi),Rx(theta)),Rz(psi))
    return R

def equator_rotation(): #Rotation matrix to bring vector to equator 
    obliquity = rad(23.4393)
    return Rx(obliquity)
    

def calculate_RA_DE(r): #Function which calculates RA and DE 
    x = r[0]
    y = r[1]
    z = r[2]
    
    RA = deg(np.arctan2(y,x))
    
    if(RA<0):
        RA = RA+360
    DEC = 90 - deg(np.arctan2(z,sqrt(x**2+y**2)))
    return (RA,DEC)

def matrix_prod(R,r):
 	prod = np.asarray(np.matmul(R,r)).reshape(-1)
 	return prod



e_earth = 0.01671 #Earth Parameters
a_earth = 1
earth_phi = rad(18.272)
earth_theta = rad(0)
earth_psi = rad(85.901)

e_jupiter = 0.04839266  #Jupiter parameters
a_jupiter = 5.20336301
jupiter_phi = rad(100.398)
jupiter_theta = rad(1.30530)
jupiter_psi = rad(273.442)


#Dictionary helpful when calculating days since Epoch
months = {"Jan": 1, "Feb":2, "Mar":3, "Apr":4, 
          "May": 5, "Jun":6,"Jul": 7, "Aug":8, 
          "Sep":9, "Oct":10, "Nov": 11, "Dec":12}

ra_alg = [] #List to store what my algorithm predicts for RA
de_alg = [] #List to store what my algorithm predicts for DE


for ele in dates: #Looping over dates
    
    Y = int(ele[0:4]) #Extracting Year,Month,Date
    M = months[ele[5:8]]
    D = int(ele[9:11])
    
    d = day((Y,M,D)) #Calculating days since Epoch

    M_jupiter = mean_anomaly_jupiter(d) #Calculating Mean anomaly for Mars
    E_jupiter = E_A((e_jupiter,M_jupiter)) #Calculating Eccentric anomaly for Mars
    pos_jupiter = position(a_jupiter,e_jupiter,E_jupiter) #Calculating orbit position vector
    R_euler_jupiter = euler_rotation(jupiter_phi,jupiter_theta,jupiter_psi) #Defining euler rotation
    eclipitic_pos_jupiter = matrix_prod(R_euler_jupiter,pos_jupiter) #Applying euler rotation
      
    
    M_earth = mean_anomaly_earth(d) #Calculating Mean anomaly for Earth
    E_earth = E_A((e_earth,M_earth)) #Calculating Eccentric anomaly for Earth
    pos_earth = position(a_earth,e_earth,E_earth) #Calculating orbit position vector
    R_euler_earth = euler_rotation(earth_phi,earth_theta,earth_psi) #Defining euler rotation
    eclipitic_pos_earth = matrix_prod(R_euler_earth,pos_earth) #Applying euler rotation
      

    r_eclip = sub_vec(eclipitic_pos_jupiter,eclipitic_pos_earth) #FInding geo-centric position of Mars 
    R_equat = equator_rotation() #Defining rotation from ellipitic to equatorial plane
    r_equat = matrix_prod(R_equat,r_eclip) #Applying the rotation
    

    RA,DE = calculate_RA_DE(r_equat) #Calculating RA and DE
    ra_alg.append(RA) #Appending RA and DE to respective lists
    de_alg.append(DE)


#Plotting results

time_scale = list(range(len(data)))

plt.style.use('dark_background')
plt.scatter(time_scale,ra_obs,color="r", label = "Actual data given")
plt.plot(time_scale[0:50],ra_alg[0:50],color = 'b', label = 'Prediction from algorithm')
plt.plot(time_scale[50:],ra_alg[50:],color = 'b')
plt.grid()
plt.legend(loc = 'upper center',fontsize='xx-large')
plt.title("Right Ascension vs Time", fontsize="xx-large")
plt.xlabel("Time",fontsize= 'xx-large')
plt.ylabel("Right Ascension",fontsize= 'xx-large')
plt.show()

plt.style.use('dark_background')
plt.scatter(time_scale,de_obs,color="r", label = "Actual data given")
plt.plot(time_scale,de_alg,color = 'b', label = 'Prediction from algorithm')
plt.grid()
plt.legend(fontsize= 'xx-large')
plt.xlabel("Time",fontsize= 'xx-large')
plt.ylabel("Declination",fontsize= 'xx-large')
plt.title("Delicnation vs Time", fontsize="xx-large")
plt.show()

#Lists to append percent error to
RA_perc_error_list = [] 
DE_perc_error_list = []

for i in range(len(ra_alg)):
    #calculating percent error
    RA_perc_error = 100*(np.abs(ra_obs[i]-ra_alg[i]))/ra_alg[i] 
    DE_perc_error = 100*(np.abs(de_obs[i]-de_alg[i]))/de_alg[i]
    
    #appending percent error to lists for all points
    RA_perc_error_list.append(RA_perc_error)
    DE_perc_error_list.append(DE_perc_error)

#computing avg percent error
mean_RA_perc_error = np.mean(RA_perc_error_list) 
mean_DE_perc_error = np.mean(DE_perc_error_list)

print(mean_RA_perc_error,mean_DE_perc_error)



