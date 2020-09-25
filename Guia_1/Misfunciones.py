# Funciones para la guía 1 de astrometria
def RNG(a, c, M, x1, N):
    """ Devuelve números en [0,1] creados con el método de congruencia lineal
    
    Parameters 
    ----------
    a : .float
        Multiplicador 
    c : .float
        Incremento
    M : .float
        Módulo
    x1 : .float
        Valor inicial
    N : int
        N° de puntos a generar
        
    Returns
    -------
    np.sdarray
        Tupla numpy de números aleatorios en [0,1]
    
    """
    import numpy as np
    # Hago el loop generador
    ij = 0 
    Aux = x1 
    X = np.empty(N)
    while ij<N:
        x0 = Aux             
        Aux = (a*x0 + c) %M
        X[ij] = Aux          
        ij = ij + 1
    # Obtengo la distribución normalizada. i.e) limito al rango [0,1]
    X = X/max(X) 
    return X

