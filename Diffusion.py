# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math

def Diff(C,T):
    
    k = 0.0000825; 
    Nc = 6.2e15 * (T ** 1.5)
    Nv = 3.5e15 * (T ** 1.5) 
    Eg = 1.17 - (0.000473 * T ** 2) / (T + 636);
    ni = (Nc * Nv) ** 0.5 * math.exp(-Eg / (2 * k * T)) 
    
    D = 0.37 * math.exp(-3.46 / (k * T)) + (0.72 * C / ni) * math.exp(-3.46 / (k * T));
    
    return D

m = 40
n = 40
# Вводные данные
C0 = 4e20
Cn = 4e17
Timez = 15 * 60
Timer = 15 * 60
Tz = 1150 + 273
Tr = 1150 + 273
Xm = 1
Ym = 3

dx = (Xm * 1e-4)/n
dy = (Ym * 1e-4)/m
dt = 1
#Формирование плоскости XoY
x = np.arange(0, Xm, Xm/n)
y = np.arange(0, Ym, Ym/m)


# Объявление массивов прогоночных коэффициентов
a = np.zeros(n)
b = np.ones(n)
d= np.ones(n)
r = np.zeros(n)
delta = np.zeros(n)
lam = np.zeros(n)
# Объявление и заполнение массива значений концентрации с учетом окна для диффузии
C = np.zeros((n,m))

for i in range(n-1):
    C[0][i] = C0
    

for i in range(0, n-1):
    for j in range(0,m-1):
        if ( y[j] <= Ym/3 or y[j] >= 2*Ym/3 ):
            C[i][j] = 0
            
        


#------------ЗАГОНКА-------------ЗАГОНКА----------ЗАГОНКА------------ЗАГОНКА-------------ЗАГОНКА----------ЗАГОНКА
for t in range(0,Timez,dt):
    for j  in range(1, m-1):
        if (y[j] > Ym / 3) and (y[j] < 2* Ym / 3):
            b[0] = 0
            a[0] = 1
            d[0] = 0
            r[0] = C0
        else:
            b[0] = 0
            a[0] = -1
            d[0] = 1
            r[0] = 0
            
        delta[0] = -d[0] / a[0]
        lam[0] = r[0] / a[0]
   
        for i in range(1, n - 1): 
            a[i] = -(2 + (dx ** 2 / (Diff(C[i, j], Tz) * dt)))
            r[i] = -dx ** 2 / dy ** 2 * (C[i, j + 1] - (2 - dy ** 2 / Diff(C[i, j], Tz) / dt) * C[i, j] + C[i, j - 1])
            delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
            lam[i] = (r[i] - b[i] * lam[i - 1]) / (a[i] + b[i] * delta[i - 1])
       

        for i in range(n - 2, 1, -1):
            C[i, j] = delta[i] * C[i + 1, j] + lam[i]
      
        
    for i in range(1, n-1):
        if (y[j] > Ym / 3) and (y[j] < 2 * Ym / 3):
            b[0] = 0
            a[0] = 1
            d[0] = 0
            r[0] = C0
        else:
            b[0] = 0
            a[0] = -1
            d[0] = 1
            r[0] = 0
        
            
        delta[0] = -d[0] / a[0]
        lam[0] = r[0] / a[0]
        
        for j in range(1, m-1):
            a[j] = -(2 + (dy ** 2 / (Diff(C[i, j], Tz) * dt)))
            r[j] = -dy ** 2 / dx ** 2 * (C[i + 1, j] - (2 - dx ** 2 / Diff(C[i, j], Tz) / dt) * C[i, j] + C[i - 1, j])
            delta[j] = -d[j] / (a[j] + b[j] * delta[j - 1])
            lam[j] = (r[j] - b[j] * lam[j - 1]) / (a[j] + b[j] * delta[j - 1])
       

        for j in range(m - 2,1,-1):
            C[i, j] = delta[j] * C[i, j + 1] + lam[j]

    
#График профиля после загонки
fig = plt.figure(figsize=(9,3), dpi=300)
ax1 = fig.add_subplot(1, 1 , 1, projection='3d')
X, Y = np.meshgrid(x, y)
ax1.view_init(azim=160, elev=25)
ax1.plot_surface(X, Y, C, cmap='inferno', edgecolor='none')
ax1.axes.xaxis.set_ticklabels([])
ax1.axes.yaxis.set_ticklabels([])
ax1.axes.zaxis.set_ticklabels([])             
    







#------------РАЗГОНКА-------------РАЗГОНКА----------РАЗГОНКА------------РАЗГОНКА-------------РАЗГОНКА----------РАЗГОНКА
for t in range(0,Timer,dt):
    for j  in range(1, m-1):
        b[0] = 0
        a[0] = -1
        d[0] = 1
        r[0] = 0
            
        delta[0] = -d[0] / a[0]
        lam[0] = r[0] / a[0]
   
        for i in range(1, n - 1): 
            a[i] = -(2 + (dx ** 2 / (Diff(C[i, j], Tr) * dt)))
            r[i] = -dx ** 2 / dy ** 2 * (C[i, j + 1] - (2 - dy ** 2 / Diff(C[i, j], Tr) / dt) * C[i, j] + C[i, j - 1])
            delta[i] = -d[i] / (a[i] + b[i] * delta[i - 1])
            lam[i] = (r[i] - b[i] * lam[i - 1]) / (a[i] + b[i] * delta[i - 1])
       

        for i in range(n - 2, 1, -1):
            C[i, j] = delta[i] * C[i + 1, j] + lam[i]
      
        
    for i in range(1, n-1):
        b[0] = 0
        a[0] = -1
        d[0] = 1
        r[0] = 0
        
            
        delta[0] = -d[0] / a[0]
        lam[0] = r[0] / a[0]
        
        for j in range(1, m-1):
            a[j] = -(2 + (dy ** 2 / (Diff(C[i, j], Tr) * dt)))
            r[j] = -dy ** 2 / dx ** 2 * (C[i + 1, j] - (2 - dx ** 2 / Diff(C[i, j], Tr) / dt) * C[i, j] + C[i - 1, j])
            delta[j] = -d[j] / (a[j] + b[j] * delta[j - 1])
            lam[j] = (r[j] - b[j] * lam[j - 1]) / (a[j] + b[j] * delta[j - 1])
       

        for j in range(m - 2,1,-1):
            C[i, j] = delta[j] * C[i, j + 1] + lam[j]
   
     
#График профиля после разгонки           
fig = plt.figure(figsize=(9,3), dpi=300)
ax2 = fig.add_subplot(1, 1 , 1, projection='3d')
X, Y = np.meshgrid(x, y)
ax2.view_init(azim=160, elev=25)
ax2.plot_surface(X, Y, C, cmap='inferno', edgecolor='none')
ax2.axes.xaxis.set_ticklabels([])
ax2.axes.yaxis.set_ticklabels([])
ax2.axes.zaxis.set_ticklabels([])  
     


#-----------PN-переход-----------PN-переход-----------PN-переход-----------PN-переход-----------PN-переход--------

junction_y = []
junction_x = []
fig = plt.figure(figsize=(9,3), dpi=300)
for j in range(m):
    for i in range(n-1, 0, -1):
        if C[j][i - 1] > Cn: 
            junction_x.append(x[i])
            junction_y.append(y[j])
            continue

ax_3 = fig.add_subplot()
ax_3.plot(junction_x,junction_y)
plt.grid()










