def Modelo(Mags, Phi, Me, alpha):
    """ Modelo para ajustar
    
    Parameters
    ----------
    Mags, ERR : list
        Magnitudes observadas 
    Phi, Me, alpha : .float, .float, .float
        Parámetros del modelo
        
    Returns
    --------
    F : list
        Valores de la función
        
    """
    import numpy as np
    M = Mags # Definición para mejor vizualización
    F = [] # Contendrá valores de la función
    ij = 0
    while ij<len(M):
        # Para que no sea tan larga la def. de "F": parto en factores a la función
        # F = f1*f2*f3
        f1 = 0.4*np.log(10)*Phi
        f2 = 10**(-0.4*(M[ij]-Me)*(alpha+1))
        f3 = np.exp( -10**(-0.4*(M[ij]-Me)) )
        
        F.append( f1*f2*f3 )
        ij = ij + 1
    return F

def Likelihood(Mags, Lum, ERR, Phi, Me, alpha):
    """ Función likelihood para el problema
    
    Parameters
    ----------
    Mags : list
        Magnitudes observadas
    Lum, ERR : list, list
        Luminosidad y sus errores asociados
    Phi, Me, alpha : .float, .float, .float
        Parámetros del modelo
        
    Returns
    --------
    LK : .float
        Valor del likelihood
    """
    import numpy as np
    import scipy.stats as st
    
    Obs = np.array(Lum)
    Calc = np.array( Modelo(Mags=Mags, Phi=Phi, Me=Me, alpha=alpha) )
    
    p = st.norm(loc=Calc, scale=ERR).pdf(Obs)
    LK = p.prod() 
    return LK
    

def PRIOR(Phi, Phimin, Phimax, Me, Memin, Memax, alpha, alphamin, alphamax):
    """Función prior, es un escalón en 3d
    
    Parameters
    ----------
    Phi, Phimin, Phimax : .float, .float, .float
        Valor del parámetro Phi y sus limitens inferiores y superiores para el escalón
    Me, Memin, Memax : .float, .float, .float
        Valor del parámetro Me y sus limitens inferiores y superiores para el escalón
    alpha, alphamin, alphamax : .float, .float, .float
        Valor del parámetro alpha y sus limitens inferiores y superiores para el escalón
        
    Returns
    --------
    Prob_norm: .float
        Mass probability function para el punto definido por Phi, Me y alpha
    
    """
    norm = abs(Phimax - Phimin) * abs(Memax - Memin)* abs(alphamax - alphamin)
    
    rPhi = (Phi < Phimax) * (Phi > Phimin) # Rango para Phi
    rMe = (Me < Memax) * (Me > Memin)
    ralpha = (alpha < alphamax) * (alpha > alphamin)
    
    Prob = 1. * rPhi * rMe * ralpha
    Prob_norm = Prob/norm # Normalizo
    
    return Prob_norm

def POSTERIOR(Mags, Lum, ERR, Phi, Phimin, Phimax, Me, Memin, Memax,
              alpha, alphamin, alphamax):
    """Devuelve el valor de la función posterior en función del likelihood y el prior
    
    Parameters
    ----------
    Tiene los mismos parámetros que las funciones Likelihood() y PRIOR()
    
    Returns
    -------
    post : .float
        Valor de la probalibidad posterior"""
    post = Likelihood(Mags, Lum, ERR, Phi, Me, alpha) * PRIOR(Phi, Phimin, Phimax, 
                                                              Me, Memin, Memax,
                                                              alpha, alphamin, alphamax)
    return post

def Ploteo(C, color, label, fig, ax):
    """Plotea las cadenas
  
    Parameters
    ----------
    C : list
        Lista que se obtiene usando la función CADENAS()
    color : string
        color de las lineas ploteadas
    label : string
        label de las lineas ploteadas
    fig, ax : "no sé"
        Son los objetos que uno define antes de hacer un plot
    Returns
    -------
    Una figura
    
    """
    import numpy as np 
    import matplotlib.pyplot as plt

    # Ploteo mis cadenas
    ax[0].plot(C[0], C[1], color=color, label=label)
    ax[1].plot(C[0], C[2], color=color)
    ax[2].plot(C[0], C[3], color=color)
    
    # Ploteo los resultados de Blanton
    ax[0].hlines(y=0.0146, xmin=C[0][0], xmax=C[0][-1], color='cyan', lw=3)
    ax[1].hlines(y=-20.83, xmin=C[0][0], xmax=C[0][-1], color='cyan', lw=3)
    ax[2].hlines(y=-1.20, xmin=C[0][0], xmax=C[0][-1], color='cyan', lw=3)

    ax[2].set_xlabel('Pasos', fontsize=20)
    ax[0].set_ylabel('Phi', fontsize=20)
    ax[1].set_ylabel('Me', fontsize=20)
    ax[2].set_ylabel('Alpha', fontsize=20)
    ax[0].legend(fontsize=15, loc=2);

def Save_chain(Steps, Phi, Me, alpha, name):
    """ Guarda las cadenas de Markov
    
    Parameters
    ----------
    Steps : list
        Pasos en las cadenas
    Phi, Me, alpha : list, list, list
       Cadenas relacionadas a los parámetros
    name : string
        Nombre del archivo a guardar (incluye extensión). Ej: 'Cadena1.txt'
        
    Returns
    -------
    Crea un archivo en el directorio en el que se guardan los resultados del MCMC
    
    """
    import numpy as np
    save_chain = open( name , "w") 
    # Loop para rellenarlo de datos
    ij = 0
    while ij<len(Steps):
        print(Steps[ij], Phi[ij], Me[ij], alpha[ij],  "\n", file=save_chain, end=' ', sep=' ') 
        ij = ij + 1
    save_chain.close()

