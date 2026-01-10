
实现代码如下(MMA):
```mathematica
Clear["`*"]
(*定义三个测试函数: 洛伦兹, 指数, 高斯*)
f[x_] := (x - 1) Exp[-1/2 (x - 1)^2];
fp[x_] := -E^(-(1/2) (-1 + x)^2) (-2 + x) x;
fpp[x_] := E^(-(1/2) (-1 + x)^2) (2 - 3 x^2 + x^3);

g[x_] := 1/(1 + (x - 5)^2);
gp[x_] := -((2 (-5 + x))/(1 + (-5 + x)^2)^2);
gpp[x_] := (2 (74 - 30 x + 3 x^2))/(26 - 10 x + x^2)^3;

h[x_] := Exp[x];
hp[x_] := Exp[x];
hpp[x_] := Exp[x];

x = Range[0., 10., 0.05];
dx = x[[2]] - x[[1]];
n = Length@x;

fnvals = f[x];
fnders = fp[x];
fndders = fpp[x];

gnvals = g[x];
gnders = gp[x];
gndders = gpp[x];

hnvals = h[x];
hnders = hp[x];
hndders = hpp[x];
```

#linebreak()

```mathematica
(*构建一二阶微分算符矩阵, 边界处使用非对称有偏差分公式. 
端点处微分不求, 实际情况需要的是边界条件, 因而微分算符矩阵的行数比函数空间维度少2*)
D1 = {};
AppendTo[D1,Join[{-3, -10, 18, -6, 1}/12/dx, Table[0., n - 5]]];
AppendTo[D1,Join[{3, -30, -20, 60, -15, 2}/60/dx, Table[0., n - 6]]];
Do[
  AppendTo[D1,ArrayPad[{-1, 9, -45, 0., 45, -9, 1}/60/dx, {i, n - i - 7}]];
  , {i, 0, n - 7}];
AppendTo[D1,Join[Table[0, n - 6], {-2, 15, -60, 20, 30, -3}/60/dx]];
AppendTo[D1,Join[Table[0, n - 5], {-1, 6, -18, 10, 3}/12/dx]];
D1 // Dimensions

D2 = {};
AppendTo[D2, Join[{11, -20, 6, 4, -1}/12/dx^2, Table[0., n - 5]]];
AppendTo[D2, Join[{-1, 16, -30, 16, -1}/12/dx^2, Table[0., n - 5]]];
Do[
  AppendTo[D2,ArrayPad[{2,-27,270,-490,270,-27,2}/180/dx^2, {i,n-i-7}]];
  , {i, 0, n - 7}];
AppendTo[D2, Join[Table[0, n - 5], {-1, 16, -30, 16, -1}/12/dx^2]];
AppendTo[D2, Join[Table[0, n - 5], {-1, 4, 6, -20, 11}/12/dx^2]];
D2 // Dimensions
```

#linebreak()

```mathematica
(*误差分析*)
fddfder = D1 . fnvals;
fddgder = D1 . gnvals;
fddhder = D1 . hnvals;
Print["一阶数值微分绝对误差: ", Min@Abs[fddfder - fnders[[2 ;; -2]]], "\t", 
 Max@Abs[fddfder - fnders[[2 ;; -2]]]]
Print["一阶数值微分绝对误差: ", Min@Abs[fddgder - gnders[[2 ;; -2]]], "\t", 
 Max@Abs[fddgder - gnders[[2 ;; -2]]]]
Print["一阶数值微分绝对误差: ", Min@Abs[fddhder - hnders[[2 ;; -2]]], "\t", 
 Max@Abs[fddhder - hnders[[2 ;; -2]]]]

fddfdder = D2 . fnvals;
fddgdder = D2 . gnvals;
fddhdder = D2 . hnvals;
Print["二阶数值微分绝对误差: ", Min@Abs[fddfdder - fndders[[2 ;; -2]]], "\t", 
  Max@Abs[fddfdder - fndders[[2 ;; -2]]]];
Print["二阶数值微分绝对误差: ", Min@Abs[fddgdder - gndders[[2 ;; -2]]], "\t", 
  Max@Abs[fddgdder - gndders[[2 ;; -2]]]];
Print["二阶数值微分绝对误差: ", Min@Abs[fddhdder - hndders[[2 ;; -2]]], "\t", 
  Max@Abs[fddhdder - hndders[[2 ;; -2]]]];

Out[]: 一阶数值微分绝对误差: 2.872963293548163*10^-20	2.573846246703426*10^-6
       一阶数值微分绝对误差: 3.552713678800501*10^-15	4.662078710304662*10^-7
       一阶数值微分绝对误差: 1.297304486058692*10^-10	0.006282381797063863
       二阶数值微分绝对误差: 5.599256121459719*10^-20	0.00008775644002922967
       二阶数值微分绝对误差: 7.202571872255703*10^-15	1.086634426883393*10^-6
       二阶数值微分绝对误差: 3.267741632839716*10^-11	0.2115191904267704
```


#linebreak()
#linebreak()
Python实现:
```py
# 测试函数及其一二阶导函数
f   = lambda x: (x - 1) * np.exp(-0.5 * (x - 1)**2)
fp  = lambda x: -np.exp(-0.5 * (x - 1)**2) * (x - 2) * x
fpp = lambda x: np.exp(-0.5 * (x - 1)**2) * (2 - 3 * x**2 + x**3)

# 离散化
x = np.arange(0., 10.+0.05, 0.05)
dx = x[1] - x[0]
n = len(x)

# 解析解
fnvals = f(x)
fnders = fp(x)
fndders = fpp(x)

# 构建微分算符矩阵
D1 = np.zeros((n-2, n))
D2 = np.zeros((n-2, n))

D1[0, 0:5] = np.array([-3, -10, 18, -6, 1]) / 12 / dx
D1[1, 0:6] = np.array([3, -30, -20, 60, -15, 2]) / 60 / dx
D1[n-4, -6:] = np.array([-2, 15, -60, 20, 30, -3]) / 60 / dx
D1[n-3, -5:] = np.array([-1, 6, -18, 10, 3]) / 12 / dx
kernel_7d1 = np.array([-1, 9, -45, 0., 45, -9, 1]) / 60 / dx
for i in range(2, n-5+1):
    D1[i, i-2:i+5] = kernel_7d1

D2[0, 0:5] = np.array([11, -20, 6, 4, -1]) / 12 / dx / dx
D2[1, 0:5] = np.array([-1, 16, -30, 16, -1]) / 12 / dx / dx
D2[n-4, -5:] = np.array([-1, 16, -30, 16, -1]) / 12 / dx / dx
D2[n-3, -5:] = np.array([-1, 4, 6, -20, 11]) / 12 / dx / dx
kernel_7d2 = np.array([2, -27, 270, -490, 270, -27, 2]) / 180 / dx / dx
for i in range(2, n-5+1):
    D2[i, i-2:i+5] = kernel_7d2


# 有限差分数值微分. 注意两个端点处的微分不予计算(掐头去尾), 因为他们是边界条件
fddder = D1 @ fnvals
fddder2 = D2 @ fnvals

print(f"一阶微分有限差分绝对误差: {np.max(abs(fddder - fnders[1:-1]))}")
print(f'二阶微分有限差分绝对误差: {np.max(abs(fddder2 - fndders[1:-1]))}')
```



#linebreak()
#line(length: 100%, stroke: (dash: "dashed"))
#linebreak()
