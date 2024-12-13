import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
from math import pi,e
from uncertainties import ufloat
from uncertainties.umath import sqrt, pow
from scipy.optimize import curve_fit

## Dados Base ##

#Todos os dados estão no formato ufloat(dado,incerteza)

P_atm = ufloat(694*133.32,0.5*133.32) # Pressão atmosférica
T_atm = ufloat(25.1,0.5) # Temperatura ambiente
g = ufloat(9.78,0.01) # Aceleração da gravidade

d1 = ufloat(7.5/100,0.025/1000) # Diâmetro do tubo
dt = ufloat(4.735/100,0.025/1000) # Diâmetro do orifício

m_esfera = ufloat(19.4/1000,0.1/1000) # Massa da esfera
D_esfera = ufloat(31.05/1000,0.025/1000) # Diâmetro da esfera
D_tubo = ufloat(34.05/1000,0.025/1000) #Diâmetro do tubo da esfera

rho_agua = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_atm.n, 'water') # Densidade da água

#Dicionário para armazenatar os dados coletados. Formato : [Frequência , Delta H, Pressão Estática]
medicoes = [ufloat(17.85,0.01) , ufloat(1.6/1000,0.025/1000) , ufloat(29.75/1000,0.025/1000)]

resultados = [] # Lista para armazenar os resultados]

def calculo_vazao_orificio(d1,dt,rho,Delta_P,mi):
    Re = ufloat(100,0) # Estimativa inicial de Re
    parada = 0.001 # Critério de parada
    erro = 100 # Valor inicial arbitrário para o erro
    betha = dt/d1 # Cálculo do Betha

    L2 = ufloat(0.47,0)
    L1 = ufloat(1,0)
    while erro > parada: #Loop iterativo
        # Equação do C
        A = pow((19000 * betha) / Re,0.8)
        M2 = (2 * L2) / (1 - betha)
        C = 0.5961 + (0.0261*pow(betha,2)) - (0.216*pow(betha,8)) + (0.000521*pow(((pow(10,6)*betha)/Re),0.7)) + ((0.0188 + 0.0063 * A) * pow(betha,3.5) * pow((pow(10,6) / Re),0.3)) + (((0.043) + (0.080 * pow(e,(-10 * L1))) - (0.123 * pow(e,(-7 * L1)))) * (1 - 0.11 * A) * (pow(betha,4) / (1 - pow(betha,4)))) - (0.031 * (M2 - (0.8 * pow(M2,1.1))) * pow(betha,1.3))
        
        #Equação da vazão mássica
        m = (C*(pi*pow(dt,2)/4)*sqrt(2*rho*Delta_P))/sqrt(1-pow(betha,4))

        #Equação do número de Reynolds
        Re_2 = (4*m)/(mi*pi*d1)

        #Cálculo do erro
        erro = abs(Re_2.n - Re.n)
        Re = Re_2
    return [C, m/rho, m, Re]

difrenca_pressao = lambda rho,delta,g: rho*delta*g # Função para cálculo da diferença de pressão
p_estatica  = lambda rho,delta,g,patm: rho*delta*g + patm # Função para cálculo da pressão estática

str_resultados = "\n ---- RESULTADOS ----\n\n -- PLACA DE ORIFÍCIO --\n\n " #String para mostrar os resultados

DeltaP = difrenca_pressao(rho_agua,medicoes[1],g) # Diferença de pressão
P_est = p_estatica(rho_agua,medicoes[2],g,P_atm) # Pressão Estática

rho_ar = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_est.n, 'air') # Densidade do Ar
print(rho_ar)
mi_ar = cp.PropsSI("V", "T", T_atm.n + 273.15, "P", P_est.n, 'air') # Viscosidade do Ar

resultados += calculo_vazao_orificio(d1,dt,rho_ar,DeltaP,mi_ar)

#Formatação da string para apresentação de resultados
str_resultados += "Frequênica = " + str(medicoes[0].n) + "; C = " + str(resultados[0]) + "; Vazão Volumétrica = " + str(resultados[1]) + "; Vazão Mássica = " + str(resultados[2]) + "; Número de Reynolds = " + str(resultados[3]) + "\n"

str_resultados += "\n -- ESFERA --\n\n "

V_esfera = 4*pi*pow(D_esfera/2,3)/3 # Volume da esfera

E = rho_ar*g*V_esfera # Força de empuxo
W = m_esfera*g # Força peso
D = W - E # Força de arrasto

def calculo_vazao_esfera(D_esfera,D_tubo,rho,mi):
    Re = ufloat(100,0) # Estimativa inicial de Re
    parada = 0.001 # Critério de parada
    erro = 100 # Valor inicial arbitrário para o erro

    while erro > parada: #Loop iterativo
        # Equação do Cd
        Cd = (24/Re) + ((2.6*(Re/5))/(1+pow(Re/5,1.52))) + ((0.411*pow(Re/(2.63*pow(10,5)),-7.94))/(1+pow(Re/(2.63*pow(10,5)),-8))) + ((0.25*Re/pow(10,6))/(1+(Re/pow(10,6))))

        #Equação da velocidade
        V = sqrt(D/(0.5*rho*Cd*pi*pow(D_esfera/2,2)))

        #Equação do número de Reynolds
        Re_2 = V*D_esfera*rho/mi

        #Cálculo do erro
        erro = abs(Re_2.n - Re.n)
        Re = Re_2

        #Cálculo da vazão mássica
        m = V*rho*pi*pow(D_tubo/2,2)
    return [Cd,V,Re,m]

resultados += calculo_vazao_esfera(D_esfera,D_tubo,rho_ar,mi_ar)

#Formatação da string para apresentação de resultados
str_resultados += "Cd = " + str(resultados[4]) + "; Velocidade = " + str(resultados[5]) + "; Vazão Mássica = " + str(resultados[7]) + "; Número de Reynolds = " + str(resultados[6]) + " ;Força de Arrasto = " + str(D) + "\n"

print(str_resultados)