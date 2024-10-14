import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
from math import pi
from uncertainties import ufloat
from uncertainties.umath import sqrt, pow
from scipy.optimize import curve_fit

## Dados Base ##

P_atm = ufloat(689*133.32,0.5*133.32)
T_atm = ufloat(25,0.5)
g = ufloat(9.78,0.01)

d1 = ufloat(7.5/100,0.005/1000)
dt = ufloat(4.735/100,0.005/1000)

rho_agua = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_atm.n, 'water')

#Formato : {Frequência : [Delta H, Pressão Estática]}
medicoes = {
    "hz15": [ufloat(0.215/100,0.005/1000) , ufloat(2.1/100,0.005/1000)],
    "hz30": [ufloat(0.77/100,0.005/1000) , ufloat(7.77/100,0.005/1000)],
    "hz45": [ufloat(1.83/100,0.005/1000) , ufloat(17.2/100,0.5/1000)],
    "hz60": [ufloat(3.32/100,0.005/1000) , ufloat(30.1/100,0.5/1000)]
}

resultados = {}

def calculo_vazao(d1,dt,rho,Delta_P,mi):
    Re = ufloat(100,0)
    parada = 0.001
    erro = 100
    betha = dt/d1
    while erro > parada:
        C = 0.5959 + (0.0312*pow(betha,2.1)) - (0.184*pow(betha,8)) + (91.71*pow(betha,2.5)*pow(Re,-0.75)) + (((0.09*pow(betha,4))/(1 - pow(betha,4)))*0.4333) - (0.0337*pow(betha,3)*0.47)
        m = (C*(pi*pow(dt,2)/4)*sqrt(2*rho*Delta_P))/sqrt(1-pow(betha,4))
        Re_2 = (4*m)/(mi*pi*d1)
        erro = abs(Re_2.n - Re.n)
        Re = Re_2
    return {"C":C,"Vazão Volumétrica":m/rho,"Número de Reynolds":Re}

difrenca_pressao = lambda rho,delta,g: rho*delta*g
p_estatica  = lambda rho,delta,g,patm: rho*delta*g + patm

str_resultados = "\n ---- RESULTADOS ----\n\n"
y = []
for frequencia in medicoes:
    DeltaP = difrenca_pressao(rho_agua,medicoes[frequencia][0],g)
    P_est = p_estatica(rho_agua,medicoes[frequencia][1],g,P_atm)
    rho_ar = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_est.n, 'air')
    mi = cp.PropsSI("V", "T", T_atm.n + 273.15, "P", P_est.n, 'air')
    resultados[frequencia] = calculo_vazao(d1,dt,rho_ar,DeltaP,mi)
    y += [resultados[frequencia]["Vazão Volumétrica"].n]
    str_resultados += frequencia + ": C = " + str(resultados[frequencia]["C"]) + "; Vazão Volumétrica = " + str(resultados[frequencia]["Vazão Volumétrica"]) + "; Número de Reynolds = " + str(resultados[frequencia]["Número de Reynolds"]) + "\n"


print(str_resultados)

# Arrays com os dados
x = np.array([15,30,45,60])
y = np.array(y)

# Função para cálculo de R²
def calcular_r2(y, y_pred):
    residuos = y - y_pred
    sq_res = np.sum(residuos**2) #soma dos resíduos ao quadrado
    sq_tot = np.sum((y - np.mean(y))**2) #soma total ao quadrado
    return 1 - (sq_res / sq_tot)

# Função para exibir a equação do polinômio
def equacao_polinomial(p):
    eq = "y = "
    for i, coef in enumerate(p):
        potencia = len(p) - i - 1
        if potencia == 0:
            eq += f"{coef}"
        elif potencia == 1:
            eq += f"{coef}x + "
        else:
            eq += f"{coef}x^{potencia} + "
    return eq

# Ajuste para diferentes graus de polinômio
melhor_R2 = -1
melhor_grau = 0
melhor_fit = None
graus = [1, 2, 3]

plt.scatter(x, y, label='Dados reais', color='blue')

for grau in graus:
    # Obter os coeficientes do ajuste polinomial
    p = np.polyfit(x, y, grau)
    
    # Gera as previsões usando os coeficientes
    y_previsto = np.polyval(p, x)
    
    # Calcula o R²
    r2 = calcular_r2(y, y_previsto)

    # Salva o melhor ajuste
    if r2 > melhor_R2:
        melhor_R2 = r2
        melhor_grau = grau
        melhor_fit = p

# Visualização do melhor ajuste
melhor_y_previsto = np.polyval(melhor_fit, x)

# Mostrar a equação do ajuste final
equacao = equacao_polinomial(melhor_fit)
print(f"Equação do ajuste final (grau {melhor_grau}):\n")
print(equacao + '\n')

plt.plot(x, melhor_y_previsto, label=f'Ajuste Polinomial (grau {melhor_grau})', color='red')
plt.title(f'Melhor Ajuste Polinomial\nGrau: {melhor_grau}, R²: {melhor_R2}\nEquação: {equacao}')
plt.legend()
plt.show()