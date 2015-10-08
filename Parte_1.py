import numpy as np
import matplotlib.pyplot as plt
"""
Parte 1
Script que implementa un algoritmo propio del metodo de runge-kutta de orden 3
para resolver la ecuacion diferencial de segundo grado de 
"El oscilador de van del Pool"

"""
def f( y, u): #ecuacion de oscilador de van der Pool 
    return -y-mu*(y**2-1)*u
    
    
def rk_o3( f, y0, u0, n, h): #metodo runge-kutta de orden 3
    """
    funcion recibe:
        f funcion al lado derecho que se iguala con la segunda derivada de "y"
        y0 condicion inicial para y
        u0 condicion inicial para la derivada de y
        h  paso de integracion
    funcion retorna:
        [y,u] dupla con arreglos de las posiciones y sus derivadas
    """
    #arreglos de ceros donde se reemplazaran por los valores calculados de 
    # posicion "y" y derivada "u".
    y = np.zeros( n )
    u = np.zeros( n )
    #condiciones iniciales
    y[0]=y0
    u[0]=u0
    #implementacion del metodo
    for i in range( n - 1 ):
        #paso de tiempo
        l1=h*u[i]
        k1=h*f( y[i], u[i])
        l2=h*(u[i]+k1/2.0)
        k2=h*f( y[i] + l1/ 2.0, u[i] + k1/ 2.0 )
        l3=h*(u[i]-k1+2*k2)      
        k3=h*f(y[i]-l1+2*l2, u[i]-k1+2*k2)
        #calculo de los siguientes valores para la posicion y la derivada
        y[i+1] = y[i] + (l1 + 4*l2 + l3)/6.0
        u[i+1] = u[i] + (k1 + 4*k2 + k3)/6.0
    return [y,u] 


mu=1.257 #parametro mu*
#condiciones iniciales
y0_1=0.1 #caso 1
y0_2=4 #caso 2
u0=0
#discreizacion de s
n_pasos = 1000
h = 20*np.pi / n_pasos
s=np.linspace(0, 20*np.pi,n_pasos)
#obtencion de las soluciones
[y1,u1]=rk_o3( f, y0_1, u0, n_pasos,h)
[y2,u2]=rk_o3( f, y0_2, u0, n_pasos,h)
#Caso 1:
#grafico de y(s) vs s
plt.figure(1)
plt.plot(s,y1)
plt.title('Grafico y(s) vs s')
plt.xlabel('s')
plt.ylabel('y(s)')
#grafico de la trayectoria dy/ds vs. y(s)
plt.figure(2)
plt.plot(y1,u1)
plt.title('Grafico dy/ds vs y(s)')
plt.xlabel('s')
plt.ylabel('y(s)')
#Caso 2:
#grafico de y(s) vs s
plt.figure(3)
plt.plot(s,y2)
plt.title('Grafico y(s) vs s')
plt.xlabel('s')
plt.ylabel('y(s)')
#grafico de la trayectoria dy/ds vs. y(s)
plt.figure(4)
plt.plot(y2,u2)
plt.title('Grafico dy/ds vs y(s)')
plt.xlabel('s')
plt.ylabel('y(s)')
