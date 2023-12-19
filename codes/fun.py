import numpy as np
from numpy.linalg import inv
from scipy.stats import norm

def VARls(y,p,inc):
    Traw=np.size(y,axis=0)
    K=np.size(y,axis=1)
    T=Traw-p
    Y=y[p:,:].T
    a=np.empty((1,K,))
    a[:]=np.NaN
    rhs=[]
    r=y
    for i in range(p):
        r=np.vstack((a,r))
        h=np.delete(r,np.s_[Traw:],0)
        if i>0:
            rhs=np.hstack((rhs,h))
        else:
            rhs=h
    Z=rhs[p:,:].T
    if inc==1:
        Z=np.vstack((np.ones(T),Z))
    Bhat = np.dot(np.dot(Y,Z.T),inv(np.dot(Z,Z.T)))
    Uhat = (Y-np.dot(Bhat,Z))
    Sigmahat = np.dot(Uhat,Uhat.T)/(T-K*p)
    corr = (np.array([Sigmahat[np.triu_indices(K, 1)[::-1]]] ).T)
    corr_mean = np.mean(np.abs(corr))
    return Bhat,Sigmahat,Uhat,Traw,K,Z,corr,corr_mean


def inf_criteria(Yraw, pmax, inc):
    FPEraw=np.zeros(pmax+1)
    AICraw=np.zeros(pmax+1)
    HQraw=np.zeros(pmax+1)
    BICraw=np.zeros(pmax+1)
    [TpPmax , K]=np.size(Yraw,axis=0), np.size(Yraw,axis=1)
    T = TpPmax - pmax
    for m in range(1,pmax+1):
        Y=Yraw[(1+(pmax-m)):,:] # same estimation length for all 
        [_,_,Uhat,_,_,_,_,_]=VARls(Y,m,inc) #function does not adjust for degrees of freedom in the VAR estimators
        Sigmahat = np.dot(Uhat,Uhat.T)/T
        FPEraw[m]=np.power(((T+K*m+1)/(T-K*m-1)),K) * np.linalg.det(Sigmahat) #Criteria Values
        AICraw[m]=np.log(np.linalg.det(Sigmahat))+(2*m*K**2)/T
        HQraw[m]=np.log(np.linalg.det(Sigmahat)) + (2 * np.log(np.log(T)))/T * m * K**2
        BICraw[m]=np.log(np.linalg.det(Sigmahat))+ (np.log(T))/T * m* K**2

    FPE=np.argmin(FPEraw[1:]) #Index of the Minimum 
    FPE=FPE+1
    AIC=np.argmin(AICraw[1:])
    AIC=AIC+1
    HQ=np.argmin(HQraw[1:])
    HQ=HQ+1
    BIC=np.argmin(BICraw[1:])
    BIC=BIC+1
    return FPE, AIC, HQ, BIC


def HVARls1(y_d, p_q, Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12):
    Traw=np.size(y_d,axis=0)
    p_d = 66
    if Q2 == 1:
        p_d=132
    if Q3 == 1:
        p_d=198
    if Q4 == 1:
        p_d=264
    if Q5 == 1:
        p_d=330
    if Q6 == 1:
        p_d=396
    if Q7 == 1:
        p_d=66*7
    if Q8 == 1:
        p_d=66*8
    if Q9 == 1:
        p_d=66*9
    if Q10 == 1:
        p_d=66*10
    if Q11 == 1:
        p_d=66*11
    if Q12 == 1:
        p_d=66*12
    K=np.size(y_d,axis=1)
    a=np.empty((1,K,))
    a[:]=np.NaN
    rhs=[]
    r=y_d
    for i in range(p_d):
        r=np.vstack((a,r))
        h=np.delete(r,np.s_[Traw:],0)
        if i>0:
            rhs=np.hstack((rhs,h))
        else:
            rhs=h
    Z=rhs[p_d:,:].T
    Y=y_d[p_d:,:].T
    C = np.zeros((K*3,Z.shape[1]))
    if Q2 == 1:
        C = np.zeros((K*4,Z.shape[1]))
    if Q3 == 1:
        C = np.zeros((K*5,Z.shape[1]))
    if Q4 == 1:
        C = np.zeros((K*6,Z.shape[1]))
    if Q5 == 1:
        C = np.zeros((K*7,Z.shape[1]))
    if Q6 == 1:
        C = np.zeros((K*8,Z.shape[1]))
    if Q7 == 1:
        C = np.zeros((K*9,Z.shape[1]))
    if Q8 == 1:
        C = np.zeros((K*10,Z.shape[1]))
    if Q9 == 1:
        C = np.zeros((K*11,Z.shape[1]))
    if Q10 == 1:
        C = np.zeros((K*12,Z.shape[1]))
    if Q11 == 1:
        C = np.zeros((K*13,Z.shape[1]))
    if Q12 == 1:
        C = np.zeros((K*114,Z.shape[1]))
    for i in range(Z.shape[1]):
        for j in range(K):
            #C[j,i]= (Z[j,i] + Z[j+4,i] + Z[j+8,i] + Z[j+12,i] + Z[j+16,i])
            
            for m in range(4):
                if m == 0:
                    M = Z[j+1*K,i]
                else:
                    M = np.hstack((M,Z[j+(m+1)*K,i]))
            C[j,i] = np.sum(M)
            
            #C[j+4,i]= (Z[j+20,i] + Z[j+24,i] + Z[j+28,i] + Z[j+32,i] + Z[j+36,i] + Z[j+40,i] + Z[j+44,i]+ Z[j+48,i] + Z[j+52,i] + Z[j+56,i] + Z[j+60,i] + Z[j+64,i] + Z[j+68,i]+Z[j+72,i] + Z[j+76,i] + Z[j+80,i] + Z[j+84,i])           
            
            for m in range(17):
                if m == 0:
                    M = Z[j+5*K,i]
                else:
                    M = np.hstack((M,Z[j+(m+5)*K,i]))
            C[j+K,i] = np.sum(M)
            
            #C[j+8,i]= (Z[j+88,i] + Z[j+92,i] + Z[j+96,i] + Z[j+100,i]+Z[j+104,i] + Z[j+108,i] + Z[j+112,i] + Z[j+116,i] + Z[j+120,i] + Z[j+124,i] + Z[j+128,i] + Z[j+132,i]+Z[j+136,i] + Z[j+140,i] + Z[j+144,i]+ Z[j+148,i] + Z[j+152,i] + Z[j+156,i] + Z[j+160,i] + Z[j+164,i] + Z[j+168,i]+Z[j+172,i] + Z[j+176,i] + Z[j+180,i] + Z[j+184,i] + Z[j+188,i] + Z[j+192,i] + Z[j+196,i] + Z[j+200,i]+Z[j+204,i] + Z[j+208,i] + Z[j+212,i] + Z[j+216,i] + Z[j+220,i] + Z[j+224,i] + Z[j+228,i] + Z[j+232,i] + Z[j+236,i] + Z[j+240,i] + Z[j+244,i]+ Z[j+248,i] + Z[j+252,i] + Z[j+256,i] + Z[j+260,i])                                      
            
            for m in range(44):
                if m == 0:
                    M = Z[j+22*K,i]
                else:
                    M = np.hstack((M,Z[j+(m+22)*K,i]))
            C[j+K*2,i] = np.sum(M)
            
            if Q2 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+66*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+66)*K,i]))
                C[j+K*3,i] = np.sum(M)
            
            if Q3 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+132*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+132)*K,i]))
                C[j+K*4,i] = np.sum(M)
            
            if Q4 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+198*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+198)*K,i]))
                C[j+K*5,i] = np.sum(M)
            
            if Q5 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+264*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+264)*K,i]))
                C[j+K*6,i] = np.sum(M)
            
            if Q6 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+330*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+330)*K,i]))
                C[j+K*7,i] = np.sum(M)
            
            if Q7 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+396*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+396)*K,i]))
                C[j+K*8,i] = np.sum(M)
            
            if Q8 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+(66*7)*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+66*7)*K,i]))
                C[j+K*9,i] = np.sum(M)
            
            if Q9 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+(66*8)*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+66*8)*K,i]))
                C[j+K*10,i] = np.sum(M)
                
            if Q10 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+(66*9)*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+66*9)*K,i]))
                C[j+K*11,i] = np.sum(M)
                
            if Q11 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+(66*10)*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+66*10)*K,i]))
                C[j+K*12,i] = np.sum(M)
                
            if Q12 == 1:
                for m in range(66):
                    if m == 0:
                        M = Z[j+(66*11)*K,i]
                    else:
                        M = np.hstack((M,Z[j+(m+66*11)*K,i]))
                C[j+K*13,i] = np.sum(M)
            
            
    Z_het = np.vstack((Z[0:K,:],C))        
    
    if p_q == 0:
        Z_het = Z_het[:K*3,:]
    
    
    Bhat = np.dot(np.dot(Y,Z_het.T),inv(np.dot(Z_het,Z_het.T)))
    Uhat = (Y-np.dot(Bhat,Z_het))
    #T=Y.shape[0]
    Sigmahat = np.dot(Uhat,Uhat.T)/(np.size(y_d,axis=0)-p_d)
    corr = (np.array([Sigmahat[np.triu_indices(K, 1)[::-1]]] ).T)
    corr_mean = np.mean(np.abs(corr))
    Bhat1 = Bhat
    Bhat2 = []
    Bhat3 = []
    Bhat4 = []
    Bhat5 = []
    Bhat6 = []
    Bhat7 = []
    Bhat8 = []
    Bhat9 = []
    Bhat10 = []
    Bhat11 = []
    Bhat12 = []

    if Q2 == 1:
        Bhat1 = Bhat[:,:K*4]
        Bhat2 = Bhat[:,K*4:K*4+K]
    if Q3 == 1:
        Bhat3 = Bhat[:,K*4+K:K*4+2*K]
    if Q4 == 1:
        Bhat4 = Bhat[:,K*4+2*K:K*4+3*K]
    if Q5 == 1:
        Bhat5 = Bhat[:,K*4+3*K:K*4+4*K]
    if Q6 == 1:
        Bhat6 = Bhat[:,K*4+4*K:K*4+5*K]
    if Q7 == 1:
        Bhat7 = Bhat[:,K*4+5*K:K*4+6*K]
    if Q8 == 1:
        Bhat8 = Bhat[:,K*4+6*K:K*4+7*K]
    if Q9 == 1:
        Bhat9 = Bhat[:,K*4+7*K:K*4+8*K]
    if Q10 == 1:
        Bhat10 = Bhat[:,K*4+8*K:K*4+9*K]
    if Q11 == 1:
        Bhat11 = Bhat[:,K*4+9*K:K*4+10*K]
    if Q12 == 1:
        Bhat12 = Bhat[:,K*4+10*K:]
    
    return Bhat1,Bhat2,Bhat3,Bhat4,Bhat5,Bhat6,Bhat7,Bhat8,Bhat9,Bhat10,Bhat11,Bhat12,Sigmahat,Uhat,Traw,K,Z_het,corr,corr_mean


def r_to_z(corr,T,conf):
    crit_val = np.zeros((2,corr.shape[0]))
    for i in range(corr.shape[0]):
        z = 0.5*np.log((1+corr[i])/(1-corr[i]))
        up = z + norm.ppf(conf)*(1/np.sqrt(T-3))
        low = z - norm.ppf(conf)*(1/np.sqrt(T-3)) 
        up = (np.exp(2*up) -1)/(np.exp(2*up) +1)
        low = (np.exp(2*low) -1)/(np.exp(2*low) +1)
        crit_val[:,i] = [up, low]
    return crit_val


def diff(x):
    return x[1:]-x[:-1]

def diff_q(x):
    return x[66:]-x[:-66]
