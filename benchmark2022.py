from cmath import tau
import copy
#import imp
import math
import numpy as np

'''
    Pure Python/Numpy implementation of the Nelder-Mead algorithm.
    Reference: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
'''



def ackley_function(x):
  #returns the point value of the given coordinate
  part_1 = -0.2*math.sqrt(0.5*((x[0])*(x[0]) + (x[1])*(x[1])))
  part_2 = 0.5*(math.cos(2*math.pi*(x[0])) + math.cos(2*math.pi*(x[1])))
  value = math.exp(1) + 20 -20*math.exp(part_1) - math.exp(part_2)
  #returning the value
  return value

#----BenchMark 2022 --------------------------------------------------------#  

def benchmark_2022(x):

    if x[0]<0.0032:
        x0=0.0032
    elif x[0]>0.127:
        x0=0.127
    else:
        x0=x[0]
    
    if x[1]<0.0025:
        x1=0.0025
    elif x[1]>0.254:
        x1=0.254
    else:
        x1=x[1]

    if x[2]<0.0025:
        x2=0.0025
    elif x[2]>0.254:
        x2=0.254
    else:
        x2=x[2]

    if x[3]<0.0032:
        x3=0.0032
    elif x[3]>0.127:
        x3=0.127
    else:
        x3=x[3]
    #print(x0)
    #print(x1)
    #print(x2)
    #print(x3)
    L=0.3556
    g=9.81
    F=2722.
    E=200.e9
    G=75.8e9
    c_viga=2935.
    c_solda=67413.
    tau_max=0.09e9
    sigma_max=0.2e9
    delta_max=0.0065
    R=((((x1)**2)/4)+((((x0)+(x2))/2)**2))**0.5
    M=F*(L+(x1/2))
    J=2*((2**0.5)*(x0)*(x1)*((((x1)**2)/12)+((x0+x2)/2)**2))
    tau_l=F/((2**0.5)*x0*x1)
    tau_2l=(M*R)/J

    #Tensão de Corte
    tau_=((tau_l**2)+(2*tau_l*tau_2l*(x1/(2*R)))+tau_2l**2)**0.5

    #Tensão Normal
    sigma_=(6*F*L)/(x3*(x2**2))

    #Deslocamento
    delta_=(4*F*(L**3))/(E*x3*(x2**3))

    #Esforço Resistente ao longo da direção t
    Fc_=((4013*E*((((x2**2)*(x3**6))/36)**0.5))/(L**2))*(1-(((x2)/(2*L))*((E/(4*G))**0.5)))
    

    #-----Restrições-------#
    g1=tau_-tau_max
    g2=sigma_-sigma_max
    g3=F-Fc_
    g4=delta_-delta_max
    g5=x0-x3

    

    #--------------#
    #r1
    if g1>0:
        r1=g1**2
    else:
        r1=0
    
    #r2
    if g2>0:
        r2=g2**2
    else:
        r2=0

    #r3
    if g3>0:
        r3=g3**2
    else:
        r3=0
    
    #r4
    if g4>0:
        r4=g4**2
    else:
        r4=0

    #r5
    if g5>0:
        r5=g5**2
    else:
        r5=0

    

#---------------------------------#

    p=1


    #Custo
    f1=(c_solda*((x0)**2)*(x1))+c_viga*(x2)*(x3)*(L+(x1))

    #Deflexão
    f2=(4*F*(L**3))/(E*(x3)*((x2)**3))

    #Ainda sem restrições
    alfa=1

    return (alfa*f1)+((1-alfa)*f2)+(p*r1)+(p*r2)+(p*r3)+(p*r4)+(p*r5)

#----------------------------------------------------------------------------------#

# Test objetive function with constrains
def fo(x):
    f=((1-x[0])**2)+2*((x[1]-(x[0]**2))**2)
    if x[0]+(3*x[1])-3<=0:
        r=0
    else:
        r=20+10*(x[0]+(3*x[1])-3)
    return f+r

#Init array
in_array=np.array([0.0107,0.0058,0.254,0.0107])

#Optimization Algorithm
def nelder_mead(f, x_start,
                step=0.55, no_improve_thr=10e-6,
                no_improv_break=10, max_iter=0,
                alpha=1., gamma=2., rho=-0.5, sigma=0.5):
    '''
        @param f (function): function to optimize, must return a scalar score
            and operate over a numpy array of the same dimensions as x_start
        @param x_start (numpy array): initial position
        @param step (float): look-around radius in initial step
        @no_improv_thr,  no_improv_break (float, int): break after no_improv_break iterations with
            an improvement lower than no_improv_thr
        @max_iter (int): always break after this number of iterations.
            Set it to 0 to loop indefinitely.
        @alpha, gamma, rho, sigma (floats): parameters of the algorithm
            (see Wikipedia page for reference)

        return: tuple (best parameter array, best score)
    '''

    # init
    dim = len(x_start)
    prev_best = f(x_start)
    no_improv = 0
    res = [[x_start, prev_best]]

    for i in range(dim):
        x = copy.copy(x_start)
        x[i] = x[i] + step
        score = f(x)
        res.append([x, score])

    # simplex iter
    iters = 0
    while 1:
        # order
        res.sort(key=lambda x: x[1])
        best = res[0][1]

        # break after max_iter
        if max_iter and iters >= max_iter:
            return res[0]
        iters += 1

        # break after no_improv_break iterations with no improvement
        print ('Best so far:', best)
        print('Iteration no:', iters)

        if best < prev_best - no_improve_thr:
            no_improv = 0
            prev_best = best
        else:
            no_improv += 1

        if no_improv >= no_improv_break:
            return res[0]

        # centroid
        x0 = [0.] * dim
        for tup in res[:-1]:
            for i, c in enumerate(tup[0]):
                x0[i] += c / (len(res)-1)

        # reflection
        xr = x0 + alpha*(x0 - res[-1][0])
        rscore = f(xr)
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]
            res.append([xr, rscore])
            continue

        # expansion
        if rscore < res[0][1]:
            xe = x0 + gamma*(x0 - res[-1][0])
            escore = f(xe)
            if escore < rscore:
                del res[-1]
                res.append([xe, escore])
                continue
            else:
                del res[-1]
                res.append([xr, rscore])
                continue

        # contraction
        xc = x0 + rho*(x0 - res[-1][0])
        cscore = f(xc)
        if cscore < res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # reduction
        x1 = res[0][0]
        nres = []
        for tup in res:
            redx = x1 + sigma*(tup[0] - x1)
            score = f(redx)
            nres.append([redx, score])
        res = nres


if __name__ == "__main__":

    print (nelder_mead(benchmark_2022, in_array))
    #print(benchmark_2022(np.array([0.0107,0.0058,0.254,0.0107])))
    #print(transform(1,2,0))
    #print(benchmark_2022(np.array([0.0032,0.02797181,0.09526017,0.0032])))