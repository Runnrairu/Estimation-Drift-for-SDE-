import numpy as np
import matplotlib.pyplot as plt

def drift(t):#ドリフトの決定
    return 0.5*t #ここではu_t=0.5tとした

def h(n,t):#カメロンマルティン空間の基底（詳しくは参考文献[1][4]を参照）
    insin=(n-0.5)*np.pi*t/T
    return np.power(2*T,0.5)*np.sin(insin)/(sigma*np.pi*(n-0.5))

def h_dot(n,t):#基底の微分を表記する（詳しくは参考文献[1][4]を参照）
    incos=(n-0.5)*np.pi*t/T
    return np.power(2/T,0.5)*np.cos(incos)/sigma
    
def lamda(n):#基底の固有値（詳しくは参考文献[1][4]を参照）
    return sigma*T/(np.pi*(n-0.5))

def X_h(k):#独自の工夫ポイント。カメロンマルティンの元の微分を定義関数で近似したうえでXの線形性に頼る
    sum=0
    for i in range(m-1):
        t=i*T/m
        sum += h_dot(k,t)*(X[i+1]-X[i])
    return sum

def Pi_n_X(t,n,X):#ここの役割については文献[1]を参照。この実装の核である
    sum=0
    for j in range(n):
        k=j+1
        lam_k=lamda(k)
        sum += X_h(k)*sigma*h(k,t)/lam_k
    return sum

def Pi_n_L2norm(n,X):
    norm=0
    for i in range(m):
        t=i*T/m
        norm += np.power(Pi_n_X(t,n,X),2)*delta_t
    return norm

T=1.0#終端時刻
sigma=1#σの値
m=1000#時間の分割数
t=0#初期時刻
X=[0]*(m+1) #確率微分方程式のサンプルパスを格納
uhat=[0]*(m+1) #推定量を格納
loss=[0]*(m) #誤差の記録用配列
delta_t = T/m  #Δt
sigma_t = np.power(delta_t,0.5)#ブラウン運動実装のため
delta_W = np.random.normal(0,sigma_t,m)#ブラウン運動実装用の乱数
for i in range(m):#サンプルパスの実装
    t=t+delta_t
    X[i+1]=X[i]+drift(t)*delta_t+sigma*delta_W[i]
    loss[i]=np.absolute(X[i]-drift(t))
plt.plot(loss)
for l in range(3):#近似の次元をあげていく
    n=l+10
    Pi_n_norm=Pi_n_L2norm(n,X)
    for i in range(m):#推定量の構成と誤差の計測
        t=i*T/m
        uhat[i]=X[i]-(n-2)*Pi_n_X(t,n,X)/Pi_n_norm
        loss[i]=np.absolute(uhat[i]-drift(t))
    plt.plot(loss)#誤差の表示
    n += 1
plt.show()
