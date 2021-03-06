import copy
import imp
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

def benchmark_2022(x):
    L=0.3556
    g=9.81
    F=2722.*g
    E=200.e9
    c_viga=2935.
    c_solda=67413.

    #Custo
    f1=(c_solda*((x[0])**2)*(x[1]))+c_viga*(x[2])*(x[3])*(L+(x[1]))

    #Deflexão
    f2=(4*F*(L**3))/(E*(x[3])*((x[2])**3))

    #Ainda sem restrições
    peso=0.5

    return peso*f1+(1-peso)*f2



#Objetive function
def f(x):
    return (x[0])**2+(x[1])**2

#Init array
in_array=np.array([0.01,0.1,0.1,0.01])

#Optimization Algorithm
def nelder_mead(f, x_start,
                step=0.1, no_improve_thr=10e-6,
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
