"""
Script que utiliza el algoritmo de RK4 predefinido en el modulo scipy.integrate
para resolver el sistema de 3 ecuaciones diferenciales de Lorenz
Ademas grafica en 3D de la solucion (atractor de Lorenz)
"""
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def lorenz(t,y):  #funcion con ecuaciones de Lorenz
    """
    funcion recibe:
    t: arreglo con la discretizacion de tiempo
    y: arreglo que corresponde a los valores de x,y z
    y[0]-> x; y[1]-> y; y[0]-> z
    funcion retorna:
    dy: lista con los valores de las derivdas de x,y z
    """      
    #parámetros         
    sigma=10.0
    beta=8.0/3.0
    rho=28.0
    dy = np.zeros(3)
    #ecuaciones de lorenz
    dy[0] = sigma*(y[1]-y[0]) 
    dy[1] = y[0]*(rho - y[2])-y[1]
    dy[2] = y[0]*y[1]-beta*y[2]
    return dy
#tiempo de integracion
t0 = 0
tfinal = 100.0
dt = 0.01 #paso del tiempo (discretizacion)
#condiciones iniciales
y0 = [1,1,0]      
#listas vacias para colocar valores de x,y,z en el tiempo correspondiente          
Y=[]; T=[]                     
#Se establece la utilizacion de RK4
r = ode(lorenz).set_integrator('dopri5')
#se establecen valores iniciales
r.set_initial_value(y0,t0)
#Se integra y se llenan listas con los valores
while r.successful() and r.t < tfinal:
    r.integrate(r.t+dt)
    Y.append(r.y)        
    T.append(r.t)

Y = array(Y)        
#Grafico  en 3D de las soluciones
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('auto')
x=Y[:,0]
y=Y[:,1]
z=Y[:,2]
ax.plot(x, y, z)
plt.title('Solucion atractor de lorenz')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.draw()
plt.show()

