import math
import cmath
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def calcSubCoef(psi, delta, angle):
    ''' Рассчитываем коэффициенты поглощения для фона -- дописать!!
    '''
    n = 0
    k = 0
    n1 = 1
    f6 = math.radians(psi)
    g6 = math.radians(delta)
    h6 = math.radians(angle)
    a = n1*n1*math.sin(h6)**2*(1+(math.tan(h6)**2*(math.cos(2*f6)**2-math.sin(g6)**2*math.sin(2*f6)**2)/(1+math.sin(2*f6)*math.cos(g6))**2))
    b = (n1**2*math.sin(h6)**2*math.tan(h6)**2*math.sin(4*f6)*math.sin(g6))/(1+math.sin(2*f6)*math.cos(g6))**2
    n3 = math.sqrt((a + math.sqrt(a**2+b**2))/2)
    k3 = -b/2/n3
    return(n3, k3)

''' рассчитываем пси и дельту субстрата '''
def calcPsi(psi1, psi2):
    return(270-(psi1+psi2)/2)

def calcDelta(delta1, delta2):
    return(((delta1+delta2)/2 - 360)*2 + 270)

def rms(a, ateor, b, bteor):
    return(math.sqrt((a-ateor)**2 + (b-bteor)**2))


''' рассчитываем пси и дельту по трехслойной модели '''
def calcTeorPsiDelta(ns, ks, n0, k0, n1, k1, n2, k2, n3, k3, angle, wavelen, thick1, thick2, thick3):
    p5 = complex(n0, k0)
    q5 = complex(n1, k1)
    r5 = complex(n2, k2)
    s5 = complex(n3, k3)
    t5 = complex(ns, ks)
    u5 = complex(math.radians(angle),0)
    v5 = cmath.sqrt(1-0j-cmath.sin(u5)**2)
    w5 = (cmath.sqrt(q5**2 - p5**2*cmath.sin(u5)**2))/q5
    x5 = (cmath.sqrt(r5**2 - p5**2*cmath.sin(u5)**2))/r5
    y5 = (cmath.sqrt(s5**2 - p5**2*cmath.sin(u5)**2))/s5
    z5 = (cmath.sqrt(t5**2 - p5**2*cmath.sin(u5)**2))/t5
    aa = (q5*v5-p5*w5)/(q5*v5+p5*w5)
    ab = (p5*v5-q5*w5)/(p5*v5+q5*w5)
    ac = (r5*w5-q5*x5)/(r5*w5+q5*x5)
    ad = (q5*w5-r5*x5)/(q5*w5+r5*x5)
    ae = (s5*x5-r5*y5)/(s5*x5+r5*y5)
    af = (r5*x5-s5*y5)/(r5*x5+s5*y5)
    ag = (t5*y5-s5*z5)/(t5*y5+s5*z5)
    ah = (s5*y5-t5*z5)/(s5*y5+t5*z5)
    ai = 2*math.pi*q5*w5*thick1/wavelen
    al = 2*math.pi*thick2*r5*x5/wavelen
    ao = 2*math.pi*thick3*s5*y5/wavelen
    ap = (ae+ag*cmath.exp((0-2j)*ao))/(1+ae*ag*cmath.exp((0-2j)*ao))
    aq = (af+ah*cmath.exp((0-2j)*ao))/(1+af*ah*cmath.exp((0-2j)*ao))
    am = (ac+ap*cmath.exp((0-2j)*al))/(1+ac*ap*cmath.exp((0-2j)*al))
    an = (ad+aq*cmath.exp((0-2j)*al))/(1+ad*aq*cmath.exp((0-2j)*al))
    aj = (aa+am*cmath.exp(ai*(0-2j)))/(1+aa*am*cmath.exp((0-2j)*ai))
    ak = (ab+an*cmath.exp((0-2j)*ai))/(1+ab*an*cmath.exp((0-2j)*ai))
    #temp = math.fabs(aj).imag/math.fabs(ak).imag

    psirad = math.atan(abs(aj)/(abs(ak)))
    ar = psirad
    if ar < 0:
        at = 180
    else:
        at = 0
    as5 = cmath.log((aj/ak)/cmath.tan(ar)).imag
    if at == 0 and as5 < 0:
        au = 360
    else:
        au = 0
    psiRes = math.degrees(psirad)
    deltaRes = as5/math.pi*180+at+au
    #deltarad = 

    return(psiRes, deltaRes)

        

df = pd.read_csv('steel3.csv')
df['psi'] = calcPsi(df['psi1'],df['psi2']) 
df['delta'] = calcDelta(df['delta1'], df['delta2'])
df['ns'] = 0
df['ks'] = 0
for i in range(1, 5):
    df.loc[(df.nomer == i), 'ns'] = calcSubCoef(df.iloc[(i-1) % 4]['psi'], df.iloc[(i-1) % 4]['delta'], 68.5)[0]
    df.loc[(df.nomer == i), 'ks'] = calcSubCoef(df.iloc[(i-1) % 4]['psi'], df.iloc[(i-1) % 4]['delta'], 68.5)[1]

'''
   ns    ks     oks: ns
1  1.24  -2.86       1.36  -2.7
2  1.3   -2.86       1.31  -2.74
3  1.28  -2.84       1.30  -2.73
4  1.34  -2.87       1.30  -2.73

'''
'''
for i in range(4):
    print(df[df.nomer == i + 1][['psi', 'delta', 'ns']])
    '''

print(df[['nomer', 'psi', 'delta', 'ns', 'ks']])
print(len(df))


def calcThick3(df):
    df['thick3'] = 0
    df['epsilon'] = 0
    for number in range(len(df)): 
        ns = df.iloc[number % 4]['ns']
        ks = df.iloc[number % 4]['ks']
        epsilon = 1000
        result = 0
        for thick3 in range(1000):
            n3 = 1.5
            k3 = -0.05
            psiexp = df.iloc[number, 7]
            deltaexp = df.iloc[number, 8]
            psiteor =(calcTeorPsiDelta(ns,ks, 1,0,1,0,1,0,n3,k3,68.5,5400,0,0,thick3))[0]
            deltateor =(calcTeorPsiDelta(ns,ks, 1,0,1,0,1,0,n3,k3,68.5,5400,0,0,thick3))[1]
            #print(psiexp, psiteor, deltaexp, deltateor, epsilon)
            epsilon1 = rms(psiexp, psiteor, deltaexp, deltateor)
            if epsilon1 < epsilon:
                epsilon = epsilon1
                result = thick3
        df.iloc[number, 11] = result 
        df.iloc[number, 12] = epsilon
            
calcThick3(df)
print(df)
for i in range(4):
    print(df[df.nomer == i + 1][['data', 'nomer', 'thick3', 'epsilon']])
#print(df.iloc[0]['psi'])
# так можно выбрать все точки с одинаковыми номерами

#calcTeorPsiDelta(1.36, -2.7, 1, 0, 1.46, 0, 1.46, 0, 1.36, 0, 68.5, 5400)

