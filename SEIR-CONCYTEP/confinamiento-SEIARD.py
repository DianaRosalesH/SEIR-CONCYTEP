import numpy as np
import tqdm
import sys

def beta(t):
    global tap
    if t<=tap:
        return beta0
    else:
        return beta0*np.exp(-(t-tap)/t0)

def pars_ini():
    global n
    g=np.random.uniform(0.55,0.2,n)
    s=np.random.uniform(0.2,0.333,n)
    R0=np.random.uniform(2,4.5,n)
    mu=0.02*np.ones(n)
    m=np.random.uniform(0.2,0.8, n)
    b=R0*g
    return b, g, s, mu, m

n=100000
N=6168883
dt=0.01

I0=np.zeros(n)
A0=np.zeros(n)
E0=np.ones(n)
S0=(N-E0)*np.ones(n)
R0=np.zeros(n)
D0=np.zeros(n)
t=0.
tap=15
t0=55

beta0, sigma, gamma, mu, m = pars_ini()

f=float(sys.argv[1])

fg=open('relajado-'+str(f)+'.dat','w')

for i in tqdm.tqdm(range(20000)):
    #k1
    Sk1=-beta(t)*S0*I0/N-f*beta(t)*S0*A0/N
    Ek1=beta(t)*S0*I0/N+f*beta(t)*S0*A0/N-sigma*E0
    Ak1=m*sigma*E0-gamma*A0
    Ik1=(1.-m)*sigma*E0-gamma*I0
    Rk1=gamma*A0+gamma*(1.-mu)*I0
    Dk1=gamma*mu*I0
    #k2
    Sk2=-beta(t+0.5*dt)*(S0+0.5*Sk1*dt)*(I0+0.5*Ik1*dt)/N-f*beta(t+0.5*dt)*(S0+0.5*Sk1*dt)*(A0+0.5*Ak1*dt)/N
    Ek2=beta(t+0.5*dt)*(S0+0.5*Sk1*dt)*(I0+0.5*Ik1*dt)/N+f*beta(t+0.5*dt)*(S0+0.5*Sk1*dt)*(A0+0.5*Ak1*dt)/N-sigma*(E0+0.5*Ek1*dt)
    Ak2=m*sigma*(E0+0.5*Ek1*dt)-gamma*(A0+0.5*Ak1*dt)
    Ik2=(1.-m)*sigma*(E0+0.5*Ek1*dt)-gamma*(I0+0.5*Ik1*dt)
    Rk2=gamma*(A0+0.5*Ak1*dt)+gamma*(1.-mu)*(I0+0.5*Ik1*dt)
    Dk2=gamma*mu*(I0+0.5*Ik1*dt)
    #k3
    Sk3=-beta(t+0.5*dt)*(S0+0.5*Sk2*dt)*(I0+0.5*Ik2*dt)/N-f*beta(t+0.5*dt)*(S0+0.5*Sk2*dt)*(A0+0.5*Ak2*dt)/N
    Ek3=beta(t+0.5*dt)*(S0+0.5*Sk2*dt)*(I0+0.5*Ik2*dt)/N+f*beta(t+0.5*dt)*(S0+0.5*Sk2*dt)*(A0+0.5*Ak2*dt)/N-sigma*(E0+0.5*Ek2*dt)
    Ak3=m*sigma*(E0+0.5*Ek2*dt)-gamma*(A0+0.5*Ak2*dt)
    Ik3=(1.-m)*sigma*(E0+0.5*Ek2*dt)-gamma*(I0+0.5*Ik2*dt)
    Rk3=gamma*(A0+0.5*Ak2*dt)+gamma*(1.-mu)*(I0+0.5*Ik2*dt)
    Dk3=gamma*mu*(I0+0.5*Ik2*dt)
    #k4
    Sk4=-beta(t+dt)*(S0+Sk3*dt)*(I0+Ik3*dt)/N-f*beta(t+dt)*(S0+Sk3*dt)*(A0+Ak3*dt)/N
    Ek4=beta(t+dt)*(S0+Sk3*dt)*(I0+Ik3*dt)/N+f*beta(t+dt)*(S0+Sk3*dt)*(A0+Ak3*dt)/N-sigma*(E0+Ek3*dt)
    Ak4=m*sigma*(E0+Ek3*dt)-gamma*(A0+Ak3*dt)
    Ik4=(1.-m)*sigma*(E0+Ek3*dt)-gamma*(I0+Ik3*dt)
    Rk4=gamma*(A0+Ak3*dt)+gamma*(1.-mu)*(I0+Ik3*dt)
    Dk4=gamma*mu*(I0+Ik3*dt)
    #actualizacion de los valores de S, E, I, R, D
    S0+=dt*(Sk1+2.*Sk2+2.*Sk3+Sk4)/6.
    E0+=dt*(Ek1+2.*Ek2+2.*Ek3+Ek4)/6.
    I0+=dt*(Ik1+2.*Ik2+2.*Ik3+Ik4)/6.
    A0+=dt*(Ak1+2.*Ak2+2.*Ak3+Ak4)/6.
    R0+=dt*(Rk1+2.*Rk2+2.*Rk3+Rk4)/6.
    D0+=dt*(Dk1+2.*Dk2+2.*Dk3+Dk4)/6.
    t+=dt
    print(t, np.mean(S0), np.mean(E0), np.mean(I0), np.mean(A0), np.mean(R0), np.mean(D0), file=fg)