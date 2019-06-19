#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import math
from collections import defaultdict
import random
import numpy as np

if len(sys.argv) < 5:
    print ('argumentos = secuencia, temp inicial, temp final, intervalo de temp')
    sys.exit(1)
secuencia = sys.argv[1]
t_ini = int(sys.argv[2])
t_fin = int(sys.argv[3])
delta_t = int(sys.argv[4])

secuencia = list(secuencia.strip())

# deltas y deltah para la secuencia de (nucleótidos)
def delta(secuencia, nn_model, extremos):
    ds = 0. # delta S en J/mol*K
    dh = 0. # delta H en J/mol
    
    for v, w in zip(secuencia[:-1], secuencia[1:]):
        h, s = nn_model[v+w]
        dh = dh + h
        ds = ds + s

    for v in [secuencia[0],secuencia[-1]]:
        h, s = extremos[v]
        dh = dh + h
        ds = ds + s
        
    return ds, dh

# probabilidad de la reaccion de ida y vuelta
def calcular_tiempos(sps, k):
    rn1 = random.uniform(0, 1)
    t_ef = fc_1 * fc_1 * (sps[0] * sps[1]) * (k[1]) + fc_1 * sps[2] * (k[2])
    t_ef = t_ef * fc
    t1 = (1/t_ef) * math.log(1/rn1)
    return t1

# calculo de probabilidad de que la óxima reacción sea formación o separación del duplex
def calcular_probabilidad(sps, k):
    # probabilidad de la reaccion de ida y vuelta
    if sps[0] + sps[1] > 0:
        fu = sps[2] / (k[0] * (sps[0] * sps[1]))
        f = fu / fc_1 # conviritiendo los nros de moléculas en concentraciones
        p1 =  1 / (1+f)  # probabilidad de que la reacción sea la de formación de dímero
    elif sps[0] + sps[1] == 0:
        p1 = 0.
    return p1

# actualizo las concentraciones de las distintas especies
def reaccion(i, sps):
    if i == 1:
        sps[0] = sps[0] - 1
        sps[1] = sps[1] - 1
        sps[2] = sps[2] + 1
    elif i == 2:
        sps[0] = sps[0] + 1
        sps[1] = sps[1] + 1
        sps[2] = sps[2] - 1
    return sps

def dat(t, sps):
    return (t, sps[0] * fc_1, sps[1] * fc_1, sps[2] * fc_1)

def gillespie(data, c_prom):

    for k in range(t_ini, t_fin, delta_t):

        t_acum = 0.
        especies = [ s1, s2, s12 ]
    
        temp = float(k)
        dg = dh - temp*ds # delta de energía libre de gibbs en J/mol*K
        fexp = -dg/(r_b*temp)
        kd = math.exp(fexp) # constante de equilibrio de hibridación s1 + s2 -> s12
        k_m = k_h / kd
        kaes=[ kd, k_h, k_m ]
       
        data[temp].append(dat(t_acum, especies))
        prom = [0, 0, 0]

        for t in range(1, pasos, 1):
            
            t1 = calcular_tiempos(especies,kaes)  # tiempo más probable hasta la próxima reacción
            p1 = calcular_probabilidad(especies,kaes)  # probabilidad de que la reacción sea la formación del duplex
            t_acum = t_acum + t1

            # elegir un número al azar y hacer el sorteo
            rn = random.uniform(0, 1)

            # modificar las concentraciones
            if rn <= p1:
                especies = reaccion(1, especies)
            else:
                especies = reaccion(2, especies)
                
            if t > p_prom:
                prom=[x+y for x,y in zip(prom,especies)]
            
            data[temp].append(dat(t_acum, especies))
    
        n_prom = pasos - p_prom
        prom = [float(x)/n_prom for x in prom]
        #total = prom[0]+prom[1]+2*prom[2]
        #prom = [ x/total for x in prom ]
        #prom = [k]+prom
        prom = dat(k, prom)

        c_prom.append(prom)

    return(data, c_prom)

# Constantes

si1 = 0.00001   # concentración inicial de la hebra 1 (M)
si2 = 0.00001   # concentración incial de la hebra 2 (M)
si12 = 0        # concentración inicial de la doble hebra (M) # constante de velocidad de formación del duplex
vol = 1*10e-15  # volumen de reacción (l)
k_h = 3*10e3    # constante cinética de hibridización (M^-1 s^-1)

r_b = 8.3  # constante de gases J/mol*K
fc = 6.02 * 10 ** (23) * vol  # conversión de M a # de moleculas en el volumen de reacción
fc_1 = (1 / fc)  # conversión de moléculas en el volumen de reacción a molar
s1 = int(si1 * fc)  # nro de moleculas hebra 1
s2 = int(si2 * fc)  # nro de moleculas hebra 2
s12 = int(si12 * fc)  # nro de moleculas duplex

pasos = 4*(s1+s2+2*s12)  # numero de pasos totales a hacer. 4 veces el número de moléculas
p_prom = int(pasos-pasos*0.25)  # promediar sobre el último cuarto de los datos

# Tabla modelo de primeros vecinos (Sta Lucía).
# Es una lista con la secuencia de vecinos de un lado y el delta h y delta s del otro

nn_model = {}
nn_model['AA'] = [-33100,-92.9]
nn_model['AT'] = [-30100,-85.4]
nn_model['TA'] = [-30100,-89.1]
nn_model['CA'] = [-35600,-95.0]
nn_model['GT'] = [-35100,-93.7]
nn_model['CT'] = [-32600,-87.9]
nn_model['GA'] = [-34300,-92.9]
nn_model['CG'] = [-44400,-113.8]
nn_model['GC'] = [-41000,-102.1]
nn_model['GG'] = [-33500,-83.3]
nn_model['TT'] = [-33100,-92.9]
nn_model['AT'] = [-30100,-85.4]
nn_model['TA'] = [-30100,-89.1]
nn_model['TG'] = [-35600,-95.0]
nn_model['AC'] = [-35100,-93.7]
nn_model['AG'] = [-32600,-87.9]
nn_model['TC'] = [-34300,-92.9]
nn_model['CG'] = [-44400,-113.8]
nn_model['GC'] = [-41000,-102.1]
nn_model['CC'] = [-33500,-83.3]

extremos = {'A': (9600.0, 17.2), 'T': (9600.0, 17.2), 'C': (400.0, 11.7), 'G': (400.0, 11.7)}
ds, dh = delta(secuencia, nn_model, extremos)

# Aca el programa 

# data va a guardar los valores como un diccionario de temperaturas donde cada temperatura tiene la lista de 
# tiempo s1 s2 s12

data_tiempo = defaultdict(list)
prom_temp = []

data_tiempo, prom_temp = gillespie(data_tiempo, prom_temp)

prom_temp = np.array(prom_temp)

for k in data_tiempo:
    data_tiempo[k] = np.array(data_tiempo[k])
