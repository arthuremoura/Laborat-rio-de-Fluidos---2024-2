import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
from math import pi,e
from uncertainties import ufloat
from uncertainties.umath import sqrt, pow
from scipy.optimize import curve_fit


P_atm = ufloat(694*133.32,0.5*133.32) # Pressão atmosférica
T_atm = ufloat(25,0.5) # Temperatura ambiente
g = ufloat(9.78,0.01) # Aceleração da gravidade




#Dados para a bomba superior. Formato : {Tipo de dado : [Dados]}
medicoes_superior = {
    "Sucção": [ufloat(0.1,0.05) , ufloat(0.1,0.05), ufloat(0.1,0.05), ufloat(0.00,0.05), ufloat(0.00,0.05)],
    "Recalque": [ufloat(1.9,0.05) , ufloat(2.1,0.05), ufloat(2.2,0.05), ufloat(2.5,0.05), ufloat(2.6,0.05)],
    "Vazão": [ufloat(8500,250) , ufloat(8000,250), ufloat(7000,250), ufloat(5000,250), ufloat(4000,250)],
    "Desnivel": ufloat(38/100,0.05/100)
}

#Dados para a bomba inferior. Formato : {Tipo de dado : [Dados]}
medicoes_inferior = {
    "Sucção": [ufloat(38/100,0.05/100) , ufloat(38.5/100,0.05/100), ufloat(39.6/100,0.05/100), ufloat(42.3/100,0.05/100), ufloat(44.6/100,0.05/100), ufloat(46.1/100,0.05/100)],
    "Recalque": [ufloat(1.85,0.025) , ufloat(1.85,0.025), ufloat(1.95,0.025), ufloat(2.2,0.025), ufloat(2.4,0.025), ufloat(2.55,0.025)],
    "Vazão": [ufloat(9000,250) , ufloat(8500,250), ufloat(8000,250), ufloat(7000,250), ufloat(5000,250), ufloat(4000,250)],
    "Desnivel": ufloat(51.5/100,0.05/100)
}

rho_agua = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_atm.n, 'water') # Densidade da água
def bomba_superior(S,R,D,rho,g): #Sucção, Recalque, Desnível
    suc = S*101325/(rho*g)
    rec = R*101325/(rho*g)
    hman = suc + rec + D
    return hman

def bomba_inferior(S,R,D,rho,g): #Sucção, Recalque, Desnível
    rec = R*101325/(rho*g)
    suc = S
    hman = rec - suc + D
    return hman

H_man_superior = []
H_man_inferior = []

for i in range(0,5):
    H_man_superior.append(bomba_superior(medicoes_superior["Sucção"][i],medicoes_superior["Recalque"][i],medicoes_superior["Desnivel"],rho_agua,g))

for i in range(0,6):
    H_man_inferior.append(bomba_inferior(medicoes_inferior["Sucção"][i],medicoes_inferior["Recalque"][i],medicoes_inferior['Desnivel'],rho_agua,g))

# Separar valores nominais para plotagem
vazao_superior = [v.n for v in medicoes_superior["Vazão"]]
vazao_inferior = [v.n for v in medicoes_inferior["Vazão"]]
h_man_superior_nominal = [h.n for h in H_man_superior]
h_man_inferior_nominal = [h.n for h in H_man_inferior]

# Dados do catálogo
H_catalogo = [22, 23, 24, 25, 26, 28]
Vazao_catalogo = [8.6, 8, 7.4, 6.6, 6.0, 4.1]  # m3/h

# Converter vazão do catálogo para L/h
vazao_catalogo_lh = [v * 1000 for v in Vazao_catalogo]

# Plotar gráfico
plt.figure(figsize=(10, 6))
plt.plot(vazao_superior, h_man_superior_nominal, label="H_man Superior", marker='o', linestyle='-')
plt.plot(vazao_inferior, h_man_inferior_nominal, label="H_man Inferior", marker='s', linestyle='--')
plt.plot(vazao_catalogo_lh, H_catalogo, label="H_man Catálogo", marker='^', linestyle='-.')

plt.title("Curvas H_man em Função da Vazão")
plt.xlabel("Vazão (L/h)")
plt.ylabel("H_man (m)")
plt.legend()
plt.grid(True)
plt.show()

print("Q e Hman inferior")
print(medicoes_inferior["Vazão"])
print(H_man_inferior)

print("Q e Hman superior")
print(medicoes_superior["Vazão"])
print(H_man_superior)