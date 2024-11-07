import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as cp
from math import pi,e
from uncertainties import ufloat
from uncertainties.umath import sqrt, pow
from scipy.optimize import curve_fit

## Dados Base ##

#Todos os dados estão no formato ufloat(dado,incerteza)

P_atm = ufloat(689*133.32,0.5*133.32) # Pressão atmosférica
T_atm = ufloat(25,0.5) # Temperatura ambiente
g = ufloat(9.78,0.01) # Aceleração da gravidade

d1 = ufloat(7.5/100,0.025/1000) # Diâmetro do tubo
dt = ufloat(4.735/100,0.025/1000) # Diâmetro do orifício

rho_agua = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_atm.n, 'water') # Densidade da água

#Dicionário para armazenatar os dados coletados. Formato : {Frequência : [Delta H, Pressão Estática]}
medicoes = {
    "hz15": [ufloat(0.215/100,0.025/1000) , ufloat(2.1/100,0.025/1000)],
    "hz30": [ufloat(0.77/100,0.025/1000) , ufloat(7.77/100,0.025/1000)],
    "hz45": [ufloat(1.83/100,0.025/1000) , ufloat(17.2/100,0.5/1000)],
    "hz60": [ufloat(3.32/100,0.025/1000) , ufloat(30.1/100,0.5/1000)]
}

resultados = {} # Dicionário para armazenar os resultados

def calculo_vazao(d1,dt,rho,Delta_P,mi):
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
    return {"C":C,"Vazão Volumétrica":m/rho,"Vazão Mássica":m,"Número de Reynolds":Re}

difrenca_pressao = lambda rho,delta,g: rho*delta*g # Função para cálculo da diferença de pressão
p_estatica  = lambda rho,delta,g,patm: rho*delta*g + patm # Função para cálculo da pressão estática

str_resultados = "\n ---- RESULTADOS ----\n\n" #String para mostrar os resultados

y = [] # Lista que será usada no ajuste de curva. Armazena os resultados da vazão volumétrica das frequências na ordem de cálculo (15,30,45,60)

for frequencia in medicoes: #Loop para cálculo de todas as frequências medidas
    DeltaP = difrenca_pressao(rho_agua,medicoes[frequencia][0],g) # Diferença de pressão
    P_est = p_estatica(rho_agua,medicoes[frequencia][1],g,P_atm) # Pressão Estática
    
    rho_ar = cp.PropsSI("D", "T", T_atm.n + 273.15, "P", P_est.n, 'air') # Densidade do Ar
    print(rho_ar)
    mi = cp.PropsSI("V", "T", T_atm.n + 273.15, "P", P_est.n, 'air') # Viscosidade do Ar
    
    resultados[frequencia] = calculo_vazao(d1,dt,rho_ar,DeltaP,mi)
    y += [resultados[frequencia]["Vazão Volumétrica"].n]

    #Formatação da string para apresentação de resultados
    str_resultados += frequencia + ": C = " + str(resultados[frequencia]["C"]) + "; Vazão Volumétrica = " + str(resultados[frequencia]["Vazão Volumétrica"]) + "; Vazão Mássica = " + str(resultados[frequencia]["Vazão Mássica"]) + "; Número de Reynolds = " + str(resultados[frequencia]["Número de Reynolds"]) + "\n"

print(str_resultados)

# Arrays com os dados para o ajuste
x = np.array([15,30,45,60])
y = np.array(y)

# Função para cálculo de R²
def calcular_r2(y, y_pred):
    residuos = y - y_pred
    sq_res = np.sum(residuos**2) # soma dos resíduos ao quadrado
    sq_tot = np.sum((y - np.mean(y))**2) # soma total ao quadrado
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
melhor_ajuste = None
graus = [1, 2, 3] # Possíveis graus de polinômio considereando a quantidade de dados

plt.scatter(x, y, label='Dados reais', color='blue')

for grau in graus: # testagem dos diferentes graus de polinômio
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
        melhor_ajuste = p

# Visualização do melhor ajuste
melhor_y_previsto = np.polyval(melhor_ajuste, x)

# Mostrar a equação do ajuste final
equacao = equacao_polinomial(melhor_ajuste)
print(f"Equação do ajuste final (grau {melhor_grau}):\n")
print(equacao + '\n')

#Plotagem do gráfico com o ajuste
plt.plot(x, melhor_y_previsto, label=f'Ajuste Polinomial (grau {melhor_grau})', color='red')
plt.title(f'Melhor Ajuste Polinomial\nGrau: {melhor_grau}, R²: {melhor_R2}')
plt.xlabel('Frequência (Hz)') 
plt.ylabel('Vazão Volumétrica ($m^3$/s)')  
plt.show()
