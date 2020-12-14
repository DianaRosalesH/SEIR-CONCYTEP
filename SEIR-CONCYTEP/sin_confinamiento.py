dt=0.01
N=6168883

I0=0
E0=50
S0=N-I0
R0=0
D0=0

R0=4.5
sigma=0.33
gamma=0.2
mu=0.02
beta0=R0*gamma
t=0
f=open('sinconfinamiento2.dat','w')

for i in range(20000):
    #k1
    Sk1=-beta0*S0*I0/N
    Ek1=beta0*S0*I0/N-sigma*E0
    Ik1=sigma*E0-gamma*I0
    Rk1=gamma*(1.-mu)*I0
    Dk1=gamma*mu*I0
    #k2
    Sk2=-beta0*(S0+0.5*Sk1*dt)*(I0+0.5*Ik1*dt)/N
    Ek2=beta0*(S0+0.5*Sk1*dt)*(I0+0.5*Ik1*dt)/N-sigma*(E0+0.5*Ek1*dt)
    Ik2=sigma*(E0+0.5*Ek1*dt)-gamma*(I0+0.5*Ik1*dt)
    Rk2=gamma*(1.-mu)*(I0+0.5*Ik1*dt)
    Dk2=gamma*mu*(I0+0.5*Ik1*dt)
    #k3
    Sk3=-beta0*(S0+0.5*Sk2*dt)*(I0+0.5*Ik2*dt)/N
    Ek3=beta0*(S0+0.5*Sk2*dt)*(I0+0.5*Ik2*dt)/N-sigma*(E0+0.5*Ek2*dt)
    Ik3=sigma*(E0+0.5*Ek2*dt)-gamma*(I0+0.5*Ik2*dt)
    Rk3=gamma*(1.-mu)*(I0+0.5*Ik2*dt)
    Dk3=gamma*mu*(I0+0.5*Ik2*dt)
    #k4
    Sk4=-beta0*(S0+Sk3*dt)*(I0+Ik3*dt)/N
    Ek4=beta0*(S0+Sk3*dt)*(I0+Ik3*dt)/N-sigma*(E0+Ek3*dt)
    Ik4=sigma*(E0+Ek3*dt)-gamma*(I0+Ik3*dt)
    Rk4=gamma*(1.-mu)*(I0+Ik3*dt)
    Dk4=gamma*mu*(I0+Ik3*dt)

    #actualizacion de los valores de S, E, I, R, D
    S0+=dt*(Sk1+2.*Sk2+2.*Sk3+Sk4)/6.
    E0+=dt*(Ek1+2.*Ek2+2.*Ek3+Ek4)/6.
    I0+=dt*(Ik1+2.*Ik2+2.*Ik3+Ik4)/6.
    R0+=dt*(Rk1+2.*Rk2+2.*Rk3+Rk4)/6.
    D0+=dt*(Dk1+2.*Dk2+2.*Dk3+Dk4)/6.
    t+=dt
    print(t, S0, E0, I0, R0, D0, file=f)