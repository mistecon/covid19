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
        #beta0 = 0.29*dt
        beta0 = 0.7*dt
        beta1 = beta0*0.07
        gamma = 0.2*dt
        delta = 0.2*dt
        #I0=40
        I0=6
        
        td = int(4/dt)
        td2 = int(14/dt)
        
        #R0list=[1,1.5,2,2.5,3,4,5]
        S = np.array(np.zeros(T+tau+td2),dtype=float)
        E = np.array(np.zeros(T+tau+td2),dtype=float)
        I = np.array(np.zeros(T+tau+td2),dtype=float)
        Iq = np.array(np.zeros(T+tau+td2),dtype=float)
        beta = np.array(np.zeros(T+tau+td2),dtype=float)
        dbegin = pd.datetime(2020,3,1)
        dld = pd.datetime(2020,3,7)
        ctime = 8
        tld = int((dld-dbegin).days/dt)+tau+td2
        for t in range(0,tld):
            beta[t]=beta0
        for t in range(tld,tld+int(ctime/dt)):
            beta[t]=beta0+(beta1-beta0)*((t-tld)/(ctime/dt))
        for t in range(tld+int(ctime/dt),len(beta)):
            beta[t]=beta1
        
        for t in range(0,tau):
            I[t]=I0
            S[t]=N
            #To get the actual number of Q we should add the value of Q at time 0.
        
        for t in range(tau-1,T-1):
            #SEIR model
            #S[t+1]=S[t]-beta[t]*(S[t]*I[t])/N
            #E[t+1]=E[t]+beta[t]*(S[t]*I[t])/N-gamma*E[t]
            #I[t+1]=I[t]+gamma*E[t]-delta*I[t]
            #Iq[t+1]=Iq[t]+delta*I[t]
            
            #delayned model
            S[t+1]=S[t]-beta[t]*(S[t]*I[t])/N
            E[t+1]=E[t]+beta[t]*(S[t]*I[t])/N-beta[t-td]*(S[t-td]*I[t-td])/N
            I[t+1]=I[t]+beta[t-td]*(S[t-td]*I[t-td])/N-beta[t-td2]*(S[t-td2]*I[t-td2])/N
            Iq[t+1]=Iq[t]+beta[t-td2]*(S[t-td2]*I[t-td2])/N


        
        #Use alternative source of patient number(data has been updated)
        #Idata=pd.read_csv("https://raw.githubusercontent.com/kaz-ogiwara/covid19/master/data/summary.csv")
        Idata=pd.read_csv("C:/Users/Hideto Kamei/Documents/10_Environment/covid19/summary.csv")
        Idata["date"] = 0
        for row in range(0,Idata.shape[0]):
            Idata.ix[row,"date"] = (pd.datetime(Idata.iloc[row,0],Idata.iloc[row,1],Idata.iloc[row,2])-pd.datetime(Idata.iloc[0,0],Idata.iloc[0,1],Idata.iloc[0,2])).days
        res = Idata[["date","PCR検査陽性者"]]

       
        res2 = pd.DataFrame(np.zeros((res.iloc[-1,0]+1,1)))
        for t in range(0,len(res2)):
            if sum(res.iloc[:,0]==t)==1:
                res2.iloc[t,0]=int(res.ix[res.iloc[:,0]==t,1])

        #using the world data instead of domestic data
        ctr = "Austria"
        positive = pd.read_csv("C:/Users/Hideto Kamei/Documents/10_Environment/covid19/world/covid-confirmed-cases-since-100th-case.csv")
        res2 = positive.ix[positive.ix[:,0]==ctr,3]
        res2 = res2[res2>10]
        res2 = res2.reset_index().iloc[:,1]

        
        Iday = np.zeros(Day+int(tau*dt)+int(td2*dt))
        Allday = np.zeros(Day+int(tau*dt)+int(td2*dt))
        
        error=0
        for d in range(0,Day-30):
            Iday[d] = Iq[int(d/dt+tau+td2)]+res2[0]
            Allday[d]=Iq[int(d/dt+tau+td2)]+E[int(d/dt+tau)]+I[int(d/dt+tau)]
        
        for d in range(0,len(res2)):
            error += (Iday[d]-float(res2.iloc[d]))**2
        
        errorlist.iloc[i,j] = error
        
        plt.plot(np.log(res2))
        plt.plot(np.log(Iday))  
        plt.plot(np.log(Allday))
        
        plt.plot(res2)
        plt.plot(Iday)
        plt.plot(Allday)
        
        plt.plot(np.array(res2)[1:]-np.array(res2)[:-1])
        plt.plot(np.array(Iday)[1:-60]-np.array(Iday)[:-61])
        #plt.plot(np.array(Iday)[1:-30]-np.array(Iday)[:-31])
        plt.plot(np.array(Iday2)[1:-30]-np.array(Iday2)[:-31])
        plt.plot(np.array(Iday3)[1:-30]-np.array(Iday2)[:-31])
        plt.plot(np.array(Allday)[1:-30]-np.array(Allday)[:-31])
        
        plt.plot(np.log(np.array(res2)[1:]-np.array(res2)[:-1]))
        plt.plot(np.log(np.array(Iday)[1:-30]-np.array(Iday)[:-31]))
        
errorlist.to_csv("C:/Users/Hideto Kamei/Documents/10_Environment/covid19/errorlist.csv")

##################################################

