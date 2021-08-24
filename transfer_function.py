import numpy as np 
import trans_vp  
from scipy.special import jv,jn_zeros
from scipy import special
import matplotlib.pyplot as plt  
from multiprocessing import Pool, Pipe, Process 
import json 
from scipy import special
 

def frozen(f, KKs, muus, rho1, phiw1, phii1, n1, layer, r, rmax, th):
    def aimf(ss): 
        def phaseP(term1, term2, term3, term4):
            coef = [term1, -term2, term3, -term4]
            root = 1/np.real(np.sqrt(np.roots(coef))) 
            return  np.abs(root[2])
            
        m_num = 200; a1 = 5/2/100 
        km = jn_zeros(0,m_num)/rmax
        fm=np.zeros(m_num,dtype = float)

        for i in range(m_num): 
            fm[i] = 2*a1*jv(1,km[i]*a1)/(rmax**2*km[i]*jv(1,rmax*km[i])**2) 
        KKi = [3.53*10**9]*n1; Kww = [2.25*10**9]*n1
        muui = [1.80*10**9]*n1; H1 = [th]*n1 
        kappas = 1*10**(-13) 
        kappai = 5*10**(-6)        
        b013 = 0 
        fn = loadcheck(ss)  
        yt, vp = trans_vp.f(np.complex(ss), KKs, KKi, muus, muui, rho1, H1, phiw1, phii1, b013,
                        kappas, kappai, Kww, km, fm, r, fn, m_num, n1, 2,  phaseP,)
        return  yt*10**9, vp

    # the transfer function is independent of the load defined below
    def loadcheck(s): 
        freq = 84*2
        fn1 = -freq*10**3*(np.exp(-np.complex(s)/(freq/10*10**3)))*np.pi/(freq**2*10.0**6 *np.pi**2 + np.complex(s)**2)  
        fn2 = freq*10**3*(np.exp(-np.complex(s)/(freq/10*5/4*10**3)))*np.pi/(freq**2*10.0**6*np.pi**2 + np.complex(s)**2)  
        fn = fn1 + fn2  
        return fn     

    def inv(i): 
        y,vp = aimf(f[i]*1j*2*np.pi)
        load = loadcheck(f[i]*1j*2*np.pi) 
        trans = np.abs(y)/np.abs(load)
        return trans   

    def spawn(f):
        def fun(pipe,x):
            pipe.send(f(x))
            pipe.close()
        return fun

    def parmap(f,X):
        pipe=[Pipe() for x in X]
        proc=[Process(target=spawn(f),args=(c,x)) for x,(p,c) in zip(X,pipe)]
        [p.start() for p in proc]
        [p.join() for p in proc]
        return [p.recv() for (p,c) in pipe] 

    yt = parmap(inv, range(len(f)))

    def velocity_p(i): 
        y, vp = aimf(f[i]*1j*2*np.pi)
        return vp

    vp = velocity_p(1)
    return np.abs(yt), vp


 


