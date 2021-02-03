# KTOptical

ラゲールガウスとエルミートガウスのパターンを2D-arrayで返す。  
```  
using Plots
using KTOptical  
  
um = 10^-6  
  
#setParamでω、ｚ、λを入れる。  
#Exportしてないので、モジュール名は頭につける。
  
KTOptical.setParam(5um, 0.0, 1.06um)  
  
x = range(-10um,10um,step = 0.1um)  
y = range(-10um,10um,step = 0.1um)  
#引数はHG_I(l,m,x,y)
#ブロードキャストを使って使用する。

I = HG_I.(1,1,x,y')  
  
gr()  
  
plot(x,y,I,fill = true)  
```