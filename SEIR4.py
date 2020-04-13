import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime

##########SEIR#################
#時間刻み幅
dt = 0.01
#シミュレーションの日数
Day = 200
T=int(Day/dt)

#β、I0を変化して、最もエラーが小さくなる組合せを採用する
betalist=pd.Series(range(0,50))/100*dt
I0list=range(0,50)
errorlist=pd.DataFrame(np.zeros([len(betalist),len(I0list)]))
for i, beta0 in enumerate(betalist):
    for j, I0 in enumerate(I0list):
        print(i,j)
        noise1 = np.random.randn(T)
        noise2 = np.random.randn(T)
        #N=9723000
        N=126500000
        drate=0.02
        
        tau =500
        #beta0 = 0.26*dt
        beta0 = 0.29*dt
        beta1 = beta0*0.2
        gamma = 0.2*dt
        delta = 0.1*dt
        #I0=40
        I0=23
        
        #R0list=[1,1.5,2,2.5,3,4,5]
        S = np.array(np.zeros(T+tau),dtype=float)
        E = np.array(np.zeros(T+tau),dtype=float)
        I = np.array(np.zeros(T+tau),dtype=float)
        Iq = np.array(np.zeros(T+tau),dtype=float)
        beta = np.array(np.zeros(T+tau),dtype=float)
        for t in range(0,int((pd.datetime(2020,4,12)-pd.datetime(2020,2,6)).days/dt)):
            beta[t]=beta0
        for t in range(int((pd.datetime(2020,4,12)-pd.datetime(2020,2,6)).days/dt),len(beta)):
            beta[t]=beta1
        
        for t in range(0,tau):
            I[t]=I0
            S[t]=N
        
        for t in range(tau-1,T-1):
            S[t+1]=S[t]-beta[t]*(S[t]*I[t])/N
            E[t+1]=E[t]+beta[t]*(S[t]*I[t])/N-gamma*E[t]
            I[t+1]=I[t]+gamma*E[t]-delta*I[t]
            Iq[t+1]=Iq[t]+delta*I[t]
        
        #Use alternative source of patient number(data has been updated)
        Idata=pd.read_csv("https://raw.githubusercontent.com/kaz-ogiwara/covid19/master/data/summary.csv")
        Idata["date"] = 0
        for row in range(0,Idata.shape[0]):
            Idata.ix[row,"date"] = (pd.datetime(Idata.iloc[row,0],Idata.iloc[row,1],Idata.iloc[row,2])-pd.datetime(Idata.iloc[0,0],Idata.iloc[0,1],Idata.iloc[0,2])).days
        res = Idata[["date","PCR検査陽性者"]]
        
        res2 = pd.DataFrame(np.zeros((res.iloc[-1,0]+1,1)))
        for t in range(0,len(res2)):
            if sum(res.iloc[:,0]==t)==1:
                res2.iloc[t,0]=int(res.ix[res.iloc[:,0]==t,1])
        
        Iday = np.zeros(Day+int(tau*dt))
        Allday = np.zeros(Day+int(tau*dt))
        
        error=0
        for d in range(0,Day-10):
            Iday[d] = Iq[int(d/dt+tau)]
            Allday[d]=Iq[int(d/dt+tau)]+E[int(d/dt+tau)]+I[int(d/dt+tau)]
        
        for d in range(0,len(res2)):
            error += (Iday[d]-float(res2.iloc[d]))**2
        
        errorlist.iloc[i,j] = error

errorlist.to_csv("C:/Users/Hideto Kamei/Documents/10_Environment/covid19/errorlist.csv")

plt.plot(np.log(res2))
plt.plot(np.log(Iday))
plt.plot(np.log(Allday))

plt.plot(res2)
plt.plot(Iday)
plt.plot(Allday)

plt.plot(np.array(res2)[1:]-np.array(res2)[:-1])
plt.plot(np.array(Iday)[1:-30]-np.array(Iday)[:-31])
plt.plot(np.array(Iday2)[1:-30]-np.array(Iday2)[:-31])
plt.plot(np.array(Allday)[1:-30]-np.array(Allday)[:-31])


##################################################

