# Estimation-of-the-Drift-for-SDE-
## ドリフト推定
　次のような[0,T]区間上の確率微分方程式を考える。  
<img src="https://latex.codecogs.com/gif.latex?dX_t=\dot{u}_tdt+\sigma&space;dW_t" />

　σは既知であるとして、このような確率微分方程式に対してuのシュタイン推定量をノンパラメトリックで構成したい。
<img src="https://latex.codecogs.com/gif.latex?E[X_t]=u_t" />
が明らかに成り立つため、これは平均値の推定と等しい。  

　d次元ガウス型確率変数(d>2)の値Xを観測したうえで、平均値のシュタイン推定量はこう書ける。  
 
<img src="https://latex.codecogs.com/gif.latex?\hat{\mu}=X+\frac{2-d}{||X||}X" />

　それに対して、連続時間ガウス過程(すなわち実数濃度の個数のガウス型確率変数)に対しては、平均値uのシュタイン推定量はこうなる。

<img src="https://latex.codecogs.com/gif.latex?\hat{u}_t=X_t+D_tlogF" />

　この理論を用いて、実際にコンピューター上でドリフトuを推定してみる。
  　ただしD_tはマリアヴァン微分である。Fの具体的な構成や証明に関しては[1]、マリアヴァン微分については[2][3]、数値計算におけるブラウン運動の構成については[4]を参照した。

## 参考文献
[1] NICOLAS PRIVAULT,ANTHONY RÉVEILLAC,STEIN ESTIMATION FOR THE DRIFT OF GAUSSIAN PROCESSES USING THE MALLIAVIN CALCULUS(2008)  
[2]重川一郎,確率解析(2008)  
[3]Giulia Di Nunno,Bernt Øksendal,Frank Proske,Malliavin Calculus for Lévy Processes with Applications to Finance(2009)  
[4]谷口説男,松本裕行,確率解析(2013)  
