from readgssi import readgssi
import matplotlib.pyplot as plt
import numpy as np

def readgprdata(input_file,zero_gpr_value):

    #hdr = header;
    #arr = array(matriz);
    #gps = gps_info;
    #zero = time-zero
    #data = matriz com os valores de amplitude obtidos a cada amostragem

    hdr, arr, gps = readgssi.readgssi(infile=input_file, zero=[zero_gpr_value])
    data = arr[0]

    #Ajustando as variaveis

    #nst = número de amostras por traço
    #nt = número traços
    #dts =  intervalo de tempo entre amostras
    #jt = tempo (s)
    #nt = numero de traços
    #d = distância (m)
    #d_trace = intervalo de espaço entre os traços (m)
    #pa = profundidade de aquisição- GPR (m)
    
    data_y = np.transpose(data)
    nst = len(data_y[0])
    nt = len(data_y)
    dts = hdr['ns_per_zsample']
    jt = nst*dts
    d = nt/(hdr['rhf_spm'])
    d_trace = d/nt
    pa = hdr['rhf_depth']

    return data_y, nst, nt, dts, jt, d, d_trace, pa

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


def migrate (data,nx,nz,dt,app,wavelet,dx,dz,v):
    
    #data = matriz de entrada
    #nx = traços
    #nz = samples
    #dt = tempo por sample
    #app = abertura de migração
    #wavelet = wavelet para correlação cruzada
    #dx = espaço entre traços
    #dz = espaço entre samples
    #v = velocidade
    
    #matriz_t1 e matriz_t2 = matrizes temporárias
    
    matriz_t1 = np.zeros([nx,2*nz-1])
    matriz_t2 = np.zeros([nx,nz])
    
    #nl = samples na wavelet
    
    nl = len(wavelet)
    
    #Realizando a correlação cruzada e a inversão do traço gerado
    
    for i in range(0,nx):
        trace_corr = np.correlate(wavelet,data[i],'full')
        matriz_t1[i] = np.block([np.zeros(nz-nl),np.flip(trace_corr)])
    
    matriz_t3 = matriz_t1*0
    
    #Corrigindo o efeito da wavelet no traço gerado
    
    for j in range (0,nz-1):
        matriz_t3[:,j] = matriz_t1[:,nz+j-1]
    
    #loop para realização da migração com correção do valor de abertura (app) 
    
    for trace in range (0,nx):
        istart = 1 + (trace - 1) - app
        iend = 1 + (trace -1) + app
        if istart < 0:
            istart = 0
        if iend > nx:
            iend = nx
        for x_val_point in range (int(istart),int(iend)):
            for point in range(0,nz):
                r = np.sqrt((x_val_point*dx - trace*dx)**2+(point*dz)**2)
                time = int(np.round(1 + r/v/dt))
                if r!=0:
                    matriz_t2[x_val_point,point] = matriz_t2[x_val_point,point] + matriz_t3[trace, time]/r
                else:
                    matriz_t2[x_val_point,point] = matriz_t2[x_val_point,point] + matriz_t3[trace, time]
    
    return matriz_t2


matriz, samples,traces, dt_samples, time, dist, ds_trace, prof  = readgprdata('Dados/FILE____346.DZT', 65)

wv, twv = mexhat(samples, time,2.6*10**9)

matriz_diff = np.block([np.diff(matriz),np.zeros([traces,1])])

matriz_mig = migrate(matriz_diff,traces,samples,dt_samples,(traces/2),wv,ds_trace,(prof/samples),(10*10**-2)/(1*10**-9))