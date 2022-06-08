from cmath import tau
import copy
#import imp
import math
import numpy as np
import matplotlib.pyplot as plt


'''
    Pure Python/Numpy implementation of the Nelder-Mead algorithm.
    Reference: https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
'''



#----f1 - Custo --------------------------------------------------------#  

def f1(x):
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
    R=((((x[1])**2)/4)+((((x[0])+(x[2]))/2)**2))**0.5
    M=F*(L+(x[1]/2))
    J=2*((2**0.5)*(x[0])*(x[1])*((((x[1])**2)/12)+((x[0]+x[2])/2)**2))
    tau_l=F/((2**0.5)*x[0]*x[1])
    tau_2l=(M*R)/J

    #Tensão de Corte
    tau_=((tau_l**2)+(2*tau_l*tau_2l*(x[1]/(2*R)))+tau_2l**2)**0.5

    #Tensão Normal
    sigma_=(6*F*L)/(x[3]*(x[2]**2))

    #Deslocamento
    delta_=(4*F*(L**3))/(E*x[3]*(x[2]**3))

    #Esforço Resistente ao longo da direção t
    Fc_=((4013*E*((((x[2]**2)*(x[3]**6))/36)**0.5))/(L**2))*(1-(((x[2])/(2*L))*((E/(4*G))**0.5)))
    

    #-----Restrições-------#
    g1=(tau_-tau_max)/(tau_max/100)
    g2=(sigma_-sigma_max)/(sigma_max/100)
    g3=(F-Fc_)/(F/100)
    g4=(delta_-delta_max)/(delta_max/100)
    g5=(x[0]-x[3])/(0.0001)

    #Restrições de domínio
    g6=(x[0]-0.0032)/0.0001
    g7=(x[0]-0.127)/0.0001
    g8=(x[3]-0.0032)/0.0001
    g9=(x[3]-0.127)/0.0001
    g10=(x[1]-0.0025)/0.0001
    g11=(x[1]-0.254)/0.0001
    g12=(x[2]-0.0025)/0.0001
    g13=(x[2]-0.254)/0.0001

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

    #r6
    if g6<0:
        r6=g6**2
    else:
        r6=0

    #r7
    if g7>0:
        r7=g7**2
    else:
        r7=0

    #r8
    if g8<0:
        r8=g8**2
    else:
        r8=0
    
    #r9
    if g9>0:
        r9=g9**2
    else:
        r9=0
    
    #r10
    if g10<0:
        r10=g10**2
    else:
        r10=0
    
    #r11
    if g11>0:
        r11=g11**2
    else:
        r11=0

    #r12
    if g12<0:
        r12=g12**2
    else:
        r12=0

    #r13
    if g13>0:
        r13=g13**2
    else:
        r13=0

#---------------------------------#

    p=1000


    #Custo
    f1=(c_solda*((x[0])**2)*(x[1]))+c_viga*(x[2])*(x[3])*(L+(x[1]))

    return (f1)+(p*r1)+(p*r2)+(p*r3)+(p*r4)+(p*r5)+(p*r6)+(p*r7)+(p*r8)+(p*r9)+(p*r10)+(p*r11)+(p*r12)+(p*r13)

#---------------------------------------------------------------------------------------------------------------------------------#

#f2 - Deflexão------------------------------#
def f2(x):
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
    R=((((x[1])**2)/4)+((((x[0])+(x[2]))/2)**2))**0.5
    M=F*(L+(x[1]/2))
    J=2*((2**0.5)*(x[0])*(x[1])*((((x[1])**2)/12)+((x[0]+x[2])/2)**2))
    tau_l=F/((2**0.5)*x[0]*x[1])
    tau_2l=(M*R)/J

    #Tensão de Corte
    tau_=((tau_l**2)+(2*tau_l*tau_2l*(x[1]/(2*R)))+tau_2l**2)**0.5

    #Tensão Normal
    sigma_=(6*F*L)/(x[3]*(x[2]**2))

    #Deslocamento
    delta_=(4*F*(L**3))/(E*x[3]*(x[2]**3))

    #Esforço Resistente ao longo da direção t
    Fc_=((4013*E*((((x[2]**2)*(x[3]**6))/36)**0.5))/(L**2))*(1-(((x[2])/(2*L))*((E/(4*G))**0.5)))
    

    #-----Restrições-------#
    g1=(tau_-tau_max)/(tau_max/100)
    g2=(sigma_-sigma_max)/(sigma_max/100)
    g3=(F-Fc_)/(F/100)
    g4=(delta_-delta_max)/(delta_max/100)
    g5=(x[0]-x[3])/(0.0001)

    #Restrições de domínio
    g6=(x[0]-0.0032)/0.0001
    g7=(x[0]-0.127)/0.0001
    g8=(x[3]-0.0032)/0.0001
    g9=(x[3]-0.127)/0.0001
    g10=(x[1]-0.0025)/0.0001
    g11=(x[1]-0.254)/0.0001
    g12=(x[2]-0.0025)/0.0001
    g13=(x[2]-0.254)/0.0001

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

    #r6
    if g6<0:
        r6=g6**2
    else:
        r6=0

    #r7
    if g7>0:
        r7=g7**2
    else:
        r7=0

    #r8
    if g8<0:
        r8=g8**2
    else:
        r8=0
    
    #r9
    if g9>0:
        r9=g9**2
    else:
        r9=0
    
    #r10
    if g10<0:
        r10=g10**2
    else:
        r10=0
    
    #r11
    if g11>0:
        r11=g11**2
    else:
        r11=0

    #r12
    if g12<0:
        r12=g12**2
    else:
        r12=0

    #r13
    if g13>0:
        r13=g13**2
    else:
        r13=0

#---------------------------------#

    p=1000


    #Deflexão
    f2=(4*F*(L**3))/(E*(x[3])*((x[2])**3))

    return (f2)+(p*r1)+(p*r2)+(p*r3)+(p*r4)+(p*r5)+(p*r6)+(p*r7)+(p*r8)+(p*r9)+(p*r10)+(p*r11)+(p*r12)+(p*r13)


#--------Função Multi-Objetivo------------------------------------------------------------------------------------------------#
def fmo(x,alfa,f1n,f2n):
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
    R=((((x[1])**2)/4)+((((x[0])+(x[2]))/2)**2))**0.5
    M=F*(L+(x[1]/2))
    J=2*((2**0.5)*(x[0])*(x[1])*((((x[1])**2)/12)+((x[0]+x[2])/2)**2))
    tau_l=F/((2**0.5)*x[0]*x[1])
    tau_2l=(M*R)/J

    #Tensão de Corte
    tau_=((tau_l**2)+(2*tau_l*tau_2l*(x[1]/(2*R)))+tau_2l**2)**0.5

    #Tensão Normal
    sigma_=(6*F*L)/(x[3]*(x[2]**2))

    #Deslocamento
    delta_=(4*F*(L**3))/(E*x[3]*(x[2]**3))

    #Esforço Resistente ao longo da direção t
    Fc_=((4013*E*((((x[2]**2)*(x[3]**6))/36)**0.5))/(L**2))*(1-(((x[2])/(2*L))*((E/(4*G))**0.5)))
    

    #-----Restrições-------#
    g1=(tau_-tau_max)/(tau_max/100)
    g2=(sigma_-sigma_max)/(sigma_max/100)
    g3=(F-Fc_)/(F/100)
    g4=(delta_-delta_max)/(delta_max/100)
    g5=(x[0]-x[3])/(0.0001)

    #Restrições de domínio
    g6=(x[0]-0.0032)/0.0001
    g7=(x[0]-0.127)/0.0001
    g8=(x[3]-0.0032)/0.0001
    g9=(x[3]-0.127)/0.0001
    g10=(x[1]-0.0025)/0.0001
    g11=(x[1]-0.254)/0.0001
    g12=(x[2]-0.0025)/0.0001
    g13=(x[2]-0.254)/0.0001

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

    #r6
    if g6<0:
        r6=g6**2
    else:
        r6=0

    #r7
    if g7>0:
        r7=g7**2
    else:
        r7=0

    #r8
    if g8<0:
        r8=g8**2
    else:
        r8=0
    
    #r9
    if g9>0:
        r9=g9**2
    else:
        r9=0
    
    #r10
    if g10<0:
        r10=g10**2
    else:
        r10=0
    
    #r11
    if g11>0:
        r11=g11**2
    else:
        r11=0

    #r12
    if g12<0:
        r12=g12**2
    else:
        r12=0

    #r13
    if g13>0:
        r13=g13**2
    else:
        r13=0

#---------------------------------#

    p=1000


    #Custo
    f1=(c_solda*((x[0])**2)*(x[1]))+c_viga*(x[2])*(x[3])*(L+(x[1]))

    #Deflexão
    f2=(4*F*(L**3))/(E*(x[3])*((x[2])**3))

    return (alfa*(f1/f1n))+((1-alfa)*(f2/f2n))+(p*r1)+(p*r2)+(p*r3)+(p*r4)+(p*r5)+(p*r6)+(p*r7)+(p*r8)+(p*r9)+(p*r10)+(p*r11)+(p*r12)+(p*r13)
    


#Init array
in_array=np.array([0.01,0.005,0.2,0.01])

#Optimization Algorithm
def nelder_mead(f, x_start, alfa, f1n, f2n,
                step=0.05, no_improve_thr=10e-9,
                no_improv_break=30, max_iter=0,
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
    prev_best = f(x_start,alfa,f1n,f2n)
    no_improv = 0
    res = [[x_start, prev_best]]

    for i in range(dim):
        x = copy.copy(x_start)
        x[i] = x[i] + step
        score = f(x,alfa,f1n,f2n)
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
        #print ('Best so far:', best)
        #print('Iteration no:', iters)

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
        rscore = f(xr,alfa,f1n,f2n)
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]
            res.append([xr, rscore])
            continue

        # expansion
        if rscore < res[0][1]:
            xe = x0 + gamma*(x0 - res[-1][0])
            escore = f(xe,alfa,f1n,f2n)
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
        cscore = f(xc,alfa,f1n,f2n)
        if cscore < res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # reduction
        x1 = res[0][0]
        nres = []
        for tup in res:
            redx = x1 + sigma*(tup[0] - x1)
            score = f(redx,alfa,f1n,f2n)
            nres.append([redx, score])
        res = nres




#---------------- Main-----------------------------------------------------------------------------------------------------------#

if __name__ == "__main__":

    #print(nelder_mead(f1, in_array))
    f1_o_sol=nelder_mead(fmo,in_array,1.,1.,1.)[0]
    f2_o_sol=nelder_mead(fmo,in_array,0.,1.,1.)[0]
    f1_o=nelder_mead(fmo,in_array,1.,1.,1.)[1]
    f2_o=nelder_mead(fmo,in_array,0.,1.,1.)[1]
    print('Minimizar Custo... Custo: %f €  /  Deflexão: %f m' %(f1(f1_o_sol),f2(f1_o_sol)))
    print('Minimizar Deflexão... Custo: %f €  /  Deflexão: %f m' %(f1(f2_o_sol),f2(f2_o_sol)))
    custos=[]
    deflexoes=[]
    j=0
    while j<=1:
        alfas_sol=nelder_mead(fmo,in_array,j,f1_o,f2_o)[0]
        cost=f1(alfas_sol)
        deflec=f2(alfas_sol)
        custos.append(cost)
        deflexoes.append(deflec)
        j=j+0.01

    plt.plot(custos,deflexoes,'g^')
    plt.show()
