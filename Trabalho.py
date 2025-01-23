import CoolProp.CoolProp as CP
import math
from scipy.optimize import fsolve

# Função para calcular o fator de atrito com a equação de Colebrook
def colebrook(Re, D, e=0):
    if Re < 2300:
        return 64 / Re  # Escoamento laminar
    elif Re >= 2300:
        # Define a função implícita de Colebrook
        def colebrook_eq(f):
            return -2.0 * math.log10((e / (3.7 * D)) + (2.51 / (Re * math.sqrt(f[0])))) - 1 / math.sqrt(f[0])
        
        # Chuta um valor inicial para f e resolve numericamente
        f_initial = [0.02]  # Inicializa como uma lista/array
        f = fsolve(colebrook_eq, f_initial)
        return f[0]  # Extrai o valor escalar do array retornado

# Função para calcular o número de Reynolds
def calc_reynolds(densidade, velocidade, diametro, viscosidade):
    return (densidade * velocidade * diametro) / viscosidade

# Função para calcular a perda de carga em um tubo (Darcy-Weisbach)
def perda_carga_tubo(f, L, D, v, g=9.81):
    return f * (L / D) * (v**2) / (2 * g)

# Função para calcular a perda localizada (coeficiente K)
def perda_carga_localizada(K, v, g=9.81):
    return K * (v**2) / (2 * g)

# Propriedades da água a 80 °C
T1 = 20 + 273.15  # Temperatura em Kelvin
P = 101325       # Pressão atmosférica em Pa
rho1 = CP.PropsSI("D", "T", T1, "P", P, "Water")  # Densidade (kg/m³)
mu1 = CP.PropsSI("V", "T", T1, "P", P, "Water")   # Viscosidade dinâmica (Pa.s)

T2 = 80 + 273.15
rho2 = CP.PropsSI("D", "T", T2, "P", P, "Water")  # Densidade (kg/m³)
mu2 = CP.PropsSI("V", "T", T2, "P", P, "Water")   # Viscosidade dinâmica (Pa.s)

# Dados gerais
Q = 1.6667e-5  # Vazão volumétrica em m³/s (60 l/h)
g = 9.81              # Gravidade (m/s²)

# Tubo interno (PVC)
D_int_tubo_interno = 0.015 - 2 * 0.0016  # Diâmetro interno do tubo interno (m)
L_tubo_interno = 0.202  # Comprimento (m)
A_int_tubo_interno = math.pi * (D_int_tubo_interno**2) / 4  # Área da seção transversal (m²)
v_tubo_interno = Q / A_int_tubo_interno  # Velocidade média (m/s)
Re_tubo_interno = calc_reynolds(rho2, v_tubo_interno, D_int_tubo_interno, mu2)
f_tubo_interno = colebrook(Re_tubo_interno, D_int_tubo_interno, e=0)  # PVC: rugosidade ~0 µm
hf_tubo_interno = perda_carga_tubo(f_tubo_interno, L_tubo_interno, D_int_tubo_interno, v_tubo_interno)

# Tubo externo (Cobre)
D_int_tubo_externo = 0.0195 - 2 * 0.00079 - 0.015 # Diâmetro interno do tubo externo (m)
L_tubo_externo = 0.182  # Comprimento (m)
A_int_tubo_externo = math.pi * (D_int_tubo_externo**2) / 4  # Área da seção transversal (m²)
v_tubo_externo = Q / A_int_tubo_externo  # Velocidade média (m/s)
Re_tubo_externo = calc_reynolds(rho2, v_tubo_externo, D_int_tubo_externo, mu2)
f_tubo_externo = colebrook(Re_tubo_externo, D_int_tubo_externo, e=0.0015e-3)  # Cobre: rugosidade ~15 µm
hf_tubo_externo = perda_carga_tubo(f_tubo_externo, L_tubo_externo, D_int_tubo_externo, v_tubo_externo)
print(hf_tubo_externo)

# Tubo externo (Cobre) (Saída)
D_int_tubo_externo = 0.011 # Diâmetro interno do tubo externo (m)
L_tubo_externo = 0.020  # Comprimento (m)
A_int_tubo_externo = math.pi * (D_int_tubo_externo**2) / 4  # Área da seção transversal (m²)
v_tubo_externo = Q / A_int_tubo_externo  # Velocidade média (m/s)
Re_tubo_externo = calc_reynolds(rho2, v_tubo_externo, D_int_tubo_externo, mu2)
f_tubo_externo = colebrook(Re_tubo_externo, D_int_tubo_externo, e=0.0015e-3)  # Cobre: rugosidade ~15 µm
hf_tubo_externo_saida = perda_carga_tubo(f_tubo_externo, L_tubo_externo, D_int_tubo_externo, v_tubo_externo)
print(hf_tubo_externo_saida)

# Mangueira
D_mangueira = 0.012 - 2*0.0005  # Diâmetro interno da mangueira (m)
L_mangueira = 20  # Comprimento total da mangueira (m)
A_mangueira = math.pi * (D_mangueira**2) / 4  # Área da seção transversal (m²)
v_mangueira = Q / A_mangueira  # Velocidade média (m/s)
Re_mangueira = calc_reynolds(rho1, v_mangueira, D_mangueira, mu1)
f_mangueira = colebrook(Re_mangueira, D_mangueira, e=0.01e-3)  # Borracha: rugosidade ~100 µm
hf_mangueira = perda_carga_tubo(f_mangueira, L_mangueira, D_mangueira, v_mangueira)
print(hf_mangueira)

# Perdas localizadas
# 2 saídas tipo borda viva
K_borda_viva_saida = 0.5
hf_borda_viva = 2 * perda_carga_localizada(K_borda_viva_saida, v_mangueira)

# 1 entrada tipo borda viva
K_borda_viva_entrada = 1.05
hf_borda_viva = perda_carga_localizada(K_borda_viva_entrada, v_mangueira)


# 4 válvulas tipo engate rápido
K_valvula_engate = 0.05
hf_valvulas_engate = 4 * perda_carga_localizada(K_valvula_engate, v_mangueira)

# 1 válvula tipo agulha
K_valvula_agulha = 10  # Coeficiente típico para válvulas agulha
hf_valvula_agulha = perda_carga_localizada(K_valvula_agulha, v_mangueira)

# 1 saída de tubo
K_saida_tubo = 1.05
hf_saida_tubo = perda_carga_localizada(K_saida_tubo, v_tubo_interno)

# Perda total
hf_total = (
    hf_tubo_interno +
    hf_tubo_externo +
    hf_mangueira +
    hf_borda_viva +
    hf_valvulas_engate +
    hf_valvula_agulha +
    hf_saida_tubo +
    hf_tubo_externo_saida
)

# Resultado final
print(f"Perda de carga total no sistema: {hf_total:.4f} m (coluna de água)")
