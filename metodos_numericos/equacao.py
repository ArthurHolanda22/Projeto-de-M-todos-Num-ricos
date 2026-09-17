import numpy as np

def f(H, g=9.81, L=4.0, t=2.5, v_alvo=5.0):
    """
    Formulação f(H) = 0 para determinar a altura H.
    v = sqrt(2*g*H) * tanh( (sqrt(2*g*H) / (2*L)) * t )
    """
    if H <= 0:
        return np.nan
        
    termo_raiz = np.sqrt(2 * g * H)
    v_calc = termo_raiz * np.tanh((termo_raiz / (2 * L)) * t)
    return v_calc - v_alvo


def df(H, g=9.81, L=4.0, t=2.5):
    """
    Derivada de f(H) em relação a H para o Método de Newton-Raphson.
    """
    if H <= 0:
        return np.nan
        
    termo_raiz = np.sqrt(2 * g * H)
    arg = (termo_raiz / (2 * L)) * t
    sech2 = 1.0 / (np.cosh(arg) ** 2)
    
    # Aplicação da regra da cadeia
    df_dH = (g / termo_raiz) * np.tanh(arg) + termo_raiz * sech2 * (g * t / (4 * L * termo_raiz))
    return df_dH