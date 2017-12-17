# Estimation-of-the-Drift-for-SDE-
## ドリフト推定
次のような確率微分方程式を考える。  
<img src="https://latex.codecogs.com/gif.latex?dX_t=\dot{u}_tdt+\sigma&space;dW_t" />

σは既知であるとして、このような確率微分方程式に対してuのシュタイン推定量をノンパラメトリックで構成したい。  
非可算無限次元のガウス型確率変数に対してのシュタイン推定は、マリアヴァン解析を用いて定式化される。
　  
d次元ガウス型確率変数(d>2)の値Xを観測したうえで、平均値のシュタイン推定量はこう書ける。
<img src="https://latex.codecogs.com/gif.latex?\hat{\mu}=X+\frac{2-d}{||X||}X" />


##参考文献
[1] NICOLAS PRIVAULT,ANTHONY RÉVEILLAC,STEIN ESTIMATION FOR THE DRIFT OF GAUSSIAN PROCESSES USING THE MALLIAVIN CALCULUS(2008)
[2]重川一郎,確率解析(2008)
[3]Giulia Di Nunno,Bernt Øksendal,Frank Proske,Malliavin Calculus for Lévy Processes with Applications to Finance(2009)
[4]谷口説男,松本裕行,確率解析(2013)
