import numpy as np
import matplotlib.pyplot as plt
from scipy import signal

import matplotlib.cm as cm

## adaptado do curso do Schuster

def model1(migi,nx,nz,ntime,dt,app,rick,dx,dz,c):
    data=np.zeros([nx,ntime]);
    nl=len(rick);
    data1=np.zeros([nx,nl+ntime-1]);
    for ixtrace in range(0,nx):
        istart=int( (ixtrace-1)-app);
        iend=int(1+ (ixtrace-1)+app);
        if istart<0:
            istart=0;
        if iend>nx:
            iend=nx;
        for ixs in range(istart-1,iend):
            for izs in range(0,nz):
                r = np.sqrt((ixtrace*dx-ixs*dx)**2+(izs*dz)**2);
                time = int(np.round( r/c/dt ));
                data[ixtrace,time] = migi[ixs,izs]/r + data[ixtrace,time];
                
        data1[ixtrace,:]=signal.convolve(data[ixtrace,:],rick)
        data[:,0:ntime]=data1[:,nl-1:nl+ntime]
    return data


def mexhat (t_a, temp_aq, cfreq):
    #mexh = wavelet Mexican Hat (Ricker)
    #
    #t_a = taxa de amostragem
    #temp_aq = tempo
    #cfreq = frequência central
    s_t = []
    t = []
    a = temp_aq/t_a
    for i in np.arange(0,temp_aq,a):
        t.append(i)
        alpha = (i-1.1/cfreq)**2*cfreq**2*np.pi**2
        s_tpp = (1-alpha)*np.exp(-alpha)
        s_t.append(s_tpp)
    return s_t, t


def migrate(cdp1,nx,nz,ntime,dt,app,rick,dx,dz,c):
    nl=len(rick);
    cdp2=np.zeros([nx,2*ntime-1]);
    migi=np.zeros([nx,nz]);
    for q in range(0,nx):
        a=np.correlate(rick,cdp1[q],'full')
        cdp2[q,:]=np.block([np.zeros(ntime-nl),np.flip(a)])
    cdp3=cdp2*0;
    for i in range (0,ntime-1):
        cdp3[:,i]=cdp2[:,ntime+i-1]
    for ixtrace in range(0,ntrace):
        istart=1 + (ixtrace-1)-app;
        iend=1+ (ixtrace-1)+app;
        if istart<0:
            istart=0
        if iend>nx:
            iend=nx
        for ixs in range(istart,iend):
            for izs in range(0,nz):
                r = np.sqrt((ixtrace*dx-ixs*dx)**2+(izs*dx)**2);
                time = int(np.round( 1 +  r/c/dt ));
                if r!=0:
                    migi[ixs,izs] = migi[ixs,izs] + cdp3[ixtrace,time]/r;
                else:
                    migi[ixs,izs] = migi[ixs,izs] + cdp3[ixtrace,time];
    
    return migi

###Criando o modelo

dx = 25; dz = 25; #dx = intervalo espacial de amostragem (assume dx=dz)
c = 1500 #c = velocidade (cm/ns)*10**2
 
#nit = 30
nx = 100; nz = 100; #(nx,nz) = dimensões do grid do modelo 
ntrace = nx #(um traço por ponto em superfície)
dt = 0.0025
time_table = 1
ntime = int(np.round(1.2*np.sqrt(nx**2+nz**2)*dx/c/dt))
app=np.round(nx/2) #app = abertura da migração

MIG=np.zeros([nx,nz]); 

for i in range(0,int(nx/3)):
    MIG[i,int(np.round(nz/3))]=1
for i in range(round(nx*2/3),nx):
    MIG[i,int(round(nz/3))]=1
for i in range(int(np.round(nx/3))+1,int(np.round(nx*2/3-1))):
    iz=15*np.sin((i-nx/3+1)*1.9*np.pi*3/nx/2);
    MIG[i,int(np.round(iz+nz/3))]=1

rick,trick=mexhat(nz,3000,dt)

######################################

#Gerando a matriz 

cdp1= model1(MIG,nx,nz,ntime,dt,app,rick,dx,dz,c)

rick,t_rick=mexhat(100,100*0.025,3)

refl1=np.block([np.diff(cdp1),np.zeros([nx,1])])

app=int(np.round(nx/2))

dado_migr=migrate(refl1,nx,nz,ntime,dt,app,rick,dx,dz,c)