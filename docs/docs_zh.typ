// #show heading: set block(above: 3em, below: 1em)
// #show heading: set text(weight: "bold")
// #let int = math.integral
// #set par(spacing: 1.15em, first-line-indent: 2.0em, justify: true, linebreaks: "simple")
// #set page(numbering: "1")
// #set heading(numbering: "1.")b
// #set text(font: "Times New Roman", size: 10pt, lang: "zh")
// #set math.equation(numbering: "(1)")

// #line(length: 100%, stroke: (dash: "dashed"))

#import "@preview/rubber-article:0.3.1": *
#import "@preview/rich-counters:0.2.1": *

// #show: great-theorems-init
#show link: text.with(fill: rgb("#cb1e5a"))

#show: article.with(
  show-header: true,
  //header-titel: "The Title of the Paper", // 会出现在页眉
  eq-numbering: "(1.1)",
  eq-chapterwise: true,
)

#set par(spacing: 1.15em, first-line-indent: 2.0em, justify: true, linebreaks: "simple")
#set math.equation(block: true)
#set heading(numbering: "1.")

#set text(lang: "zh")
#set text(size: 10pt)
#set text(font: "Times New Roman")
#show heading: set block(above: 2em, below: 1em)
#show heading: set text(weight: "bold")
#show raw.where(block: true): block.with(
  fill: luma(240),
  inset: 10pt,
  radius: 4pt,
  width: 100%,
)



#let hbar = math.planck
#let int = math.integral
#let bvec(arg) = {$bold(arrow(#arg))$}
#let dirac(arg) = {$chevron.l #arg chevron.r$}
#let diracc(bra, ket) = {$chevron.l #bra mid(|) ket #chevron.r$}
#let diraccc(bra, H, ket) = {$chevron.l #bra mid(|) #H mid(|) #ket chevron.r$}
#let indent=h(2em)
#let mathcounter = rich-counter(
  identifier: "mathblocks",
  inherited_levels: 1
)



#set par(first-line-indent: 0em)
#line(length: 50%, stroke: white)
lyuhang428\@gmail.com
#line(length: 50%, stroke: white)


#heading(numbering: none, level: 1)[*符号约定*]
#table(
    columns: 2,
    stroke: none,
    align: left,
    [eDFT], [电子密度泛函理论], 
    [cDFT], [经典密度泛函理论], 
    [$psi(bvec(r))$, $f(bvec(r))$], [波函数],
    [$phi.alt(bvec(r))$], [基函数], 
    [$U(bvec(r))$], [库伦势能], 
    [$u(r)$], [库伦势球谐展开系数], 
    [$rho(bvec(r))$], [电子密度(eDFT)，粒子数密度(cDFT)]

)





#outline()


#set par(first-line-indent: 2em)
/////////////////
/// main text ///
/// /////////////
= *数学基础*
== 高斯型基函数
=== 高斯乘积定理
$
& g_1(arrow(r)) = e^(-alpha |arrow(r)-bold(R_A)|^2) \ 
& g_2(arrow(r)) = e^(-beta |arrow(r)-bold(R_B)|^2) \ 
& g_1 dot.c g_2 = e^(-(alpha beta)/(alpha+beta) |bold(R_A)-bold(R_B)|^2) e^(-(alpha+beta) space |arrow(r) - (alpha bold(R_A)+beta bold(R_B))/(alpha+beta)|^2) \ 
$

定义一些参数: 
$
cases(
  K = e^(-(alpha beta)/(alpha+beta) |bold(R_A)-bold(R_B)|^2), 

  p = alpha+beta, 

  bold(R_C) = (alpha bold(R_A) + beta bold(R_B))/(alpha+beta)
)
$


当高斯函数的中心不再是原点时, 被积函数的形式需要稍作调整: 
$
    int_R^infinity e^(-alpha |bvec(r) - bvec(R)|^2) (r-R)^2 sin theta d r d theta d phi <=> int_0^infinity e^(-alpha bvec(r)^2) r^2 d r d theta d phi \ 
$

如果将该积分在直角坐标下计算, 高斯函数的中心在哪里无所谓, 因为所有的积分变量的积分范围为负无穷到正无穷, 积分区间对于任何中心都对称. 而在球坐标下, 积分变量$r$的积分范围为零到正无穷, 该区间对函数中心不再是对称的, 最终会给出高斯误差函数的错误结果. 

函数中心改变不过是平移变换, 平移变换不改变积分结果（曲线下面积）. 

但是如果是两个高斯函数乘在一起, 需要额外考虑由于高斯乘积定理产生出的额外的常数项:
$
& g_1 dot.c g_2 = e^(-(alpha beta)/(alpha+beta) |bold(R_A)-bold(R_B)|^2) e^(-(alpha+beta) space |bold(r) - bold(R_C)|^2) \ 

& int g_1 dot.c g_2 space \d^3 bold(r) = K_(A B)^(alpha beta) int_(R_C)^infinity  e^(-p_(alpha beta) space |bold(r) - bold(R_C)|^2) (r-R_C)^2 \dr int \dOmega \ 

& <=> K_(A B)^(alpha beta) int_0^infinity e^(-p r^2) r^2 \dr int \dOmega \ 
$

也就是说, 在使用高斯求积计算球坐标下的积分时, 无需考虑由高斯乘积定理带来的新的函数中心, 因为两条高斯函数相乘后的新的函数中心永远可以视为新的原点. 唯一需要显示计算的是前置系数$K_(A B)^(alpha beta)$, 该系数又依赖于两条高斯函数各自的指数和中心.



=== 纯基函数与笛卡尔基函数
#linebreak()

Primitive function (Cartesian form): $(x-A_x)^(l_x) (y - A_y)^(l_y) (z - A_z)^(l_z) e^(-alpha |bvec(r) - bvec(R)|^2)$

Contracted function (Cartesian form): $sum_i^n c_i N_i (x-A_x)^(l_(x,i)) (y - A_y)^(l_(y,i)) (z - A_z)^(l_(z,i)) e^(-alpha_i |bvec(r) - bvec(R)|^2)$

通常来说同一组收缩中的所有组分使用同一套角动量：$l_(alpha, i) == l_(alpha, j) ; space alpha=x,y,x$. 当显式线性组合成球谐基函数时不同组分会使用不同的角动量. 例如 $d_(x^2-y^2)$: 
$
    d_(x^2-y^2) = (sqrt(3))/(2)(d_(x\x) - d_(y\y))
$

即，球谐形式的 $d_(x^2-y^2)$ 组分中有两组角动量. 但是实操中不会在事先就进行线性组合，而是在笛卡尔基函数上进行操作（例如计算基函数在网格上的函数值），最后通过转换矩阵得到球谐基函数. 


基函数(原子轨道), 是角度部分与径向部分的乘积: $phi.alt(bold(r)) = Y_l^m (theta,phi.alt) R(r)$. 其中角度部分是球谐函数, 而径向部分是指数函数. 对于角度部分, 球谐函数是拉普拉斯方程在球坐标下经过变量分离后角度部分的解, 通常是复的:
$
& nabla^2 f = 1/r^2 (partial)/(partial r) (r^2 (partial f)/(partial r)) + 1/(r^2 sin theta) (partial)/(partial theta) (sin theta (partial f)/(partial theta)) + 1/(r^2 sin^2 theta) (partial^2 f)/(partial phi.alt^2) = 0 \ 

& f(bvec(r)) = R(r) Y_l^m (theta, phi.alt) \ 
& cases(
1/R d/(d\r) (r^2 (d R)/(d\r)) = lambda ,
1/Y 1/(sin theta) partial/(partial theta) (sin theta (partial Y)/(partial theta)) + 1/Y 1/(sin^2 theta) (partial^2 Y)/(partial phi.alt^2) = -lambda
) \ 
$

上式中的变量分离提示我们对于一个实空间函数 $rho(bvec(r))$，我们可以将其拆解成径向部分和角度部分，这实际上就是球谐展开，此时球谐函数作为基函数，而径向部分则是展开系数. 即使径向变量和角度变量耦合，这种展开也可以进行. 

角动量算符：
$
& cases(
bold(hat(L))^2 = -hbar^2 (1/(sin theta) partial/(partial theta) (sin theta (partial)/(partial theta)) + 1/(sin^2 theta) (partial^2)/(partial phi.alt^2)),

bold(hat(L))^2 Y_l^m = hbar^2 l(l+1) Y_l^m ,
) \ 

& cases(
    bold(hat(L)_z) = -i hbar partial/(partial phi.alt) ,
    bold(hat(L)_z) Y_l^m = m hbar Y_l^m\, space m = -l\, -l+1\, dots.c\, l ,
) \ 

& nabla^2 = 1/r^2 partial/(partial r) (r^2 partial/(partial r)) - 1/(hbar^2 r^2) bold(hat(L)^2) \ 
$

球谐函数解析表达式：
$
& Y_l^m (theta, phi.alt) = sqrt((2l+1)/(4pi) ((l-m)!)/((l+m)!)) P_l^m (cos theta) e^(i m phi.alt) \ 

& int_0^(2pi) int_0^pi (Y_(l^prime)^(m^prime))^ast Y_l^m sin theta d theta d phi.alt = delta_(l l^prime) delta_(m m^prime) \ 
$

但是无论是笛卡尔基函数还是球谐基函数，都不会直接使用球谐函数（包括线性组合后的实球谐函数），二是使用笛卡尔坐标下的多项式等价表示（包括球谐基函数）. 对应关系见 #link("https://en.wikipedia.org/wiki/Table_of_spherical_harmonics"). 


形式上: 
$
& cases(
phi.alt(r,theta,phi) = Y_l^m (theta,phi) |bold(r)|^l sum_i c_i n_i e^(-alpha_i |bold(r)|^2) ,

phi.alt(r,theta,phi\;bold(R)) = Y_l^m (theta,phi) |bold(r)-bold(R)|^l sum_i c_i n_i e^(-alpha_i |bold(r)-bold(R)|^2) , 
) \ 

& cases(
phi.alt(x,y,z) = x^(l_\x) y^(l_y) z^(l_z) sum_i c_i n_i e^(-alpha_i (x^2+y^2+z^2)) ,

phi.alt(x,y,z\;a_x,\a_y,\a_z) = (x-\a_x)^(\l_x) (y-\a_y)^(\l_y) (z-\a_z)^(\l_z) sum_i c_i n_i e^(-alpha_i [(x-\a_x)^2+(y-\a_y)^2+(z-\a_z)^2])
) \ 
$

注意球谐函数形式在实际计算中不会使用. 

在纯基函数(pure basis function)表达式中取的是向量的模$|bvec(r)-bvec(R)|$；而在笛卡尔基函数表达式*不能*取绝对值 $times.o |x-\a_x|$. 

纯基函数对任意角动量 $l$ 不存在冗余, 例如对于$l=2$共有5条d轨道, 对于$l=3$共有7条f轨道, 以此类推. 因为在构造时纯基函数的角度部分直接取的是球谐函数, 自然不存在冗余. 

而对于笛卡尔基函数, 情况有所不同. 此时的角度部分是多项式$x^(l_\x) y^(\l_y) z^(\l_z)$. 当角动量$l>=2$时, 笛卡尔基函数相较于纯基函数存在冗余, 且角动量越高冗余越多: 

#table(
    columns: 4,
    column-gutter: 1em,
    row-gutter: 2em,
    stroke: none,
    [$l$], [Cartesian], [Pure], [],
    [2], [$\xx, \xy, \xz, \yy, \yz, \zz$], [$\xy, \yz, z^2, \xz, x^2-y^2$], [$6d->5d$],
    [3], [$x x x, x x y, x x z, x y y, x y z, \ x z z, y y y, y y z, y z z, z z z$], [$y(3x^2-y^2), x y z, y z^2, z^3, x z^2, \ z(x^2-y^2), x(x^2-3y^2)$], [$10f->7f$],
    [$dots.v$]
)

笛卡尔基函数角动量的排序通常采用字典顺序（lexicographic, 纯基函数的顺序常用的有两种，见下文）: ${x\x, x\y, x\z, y\y, y\z, z\z}$, ${x\x\x, x\x\y, x\x\z, x\y\y, x\y\z, x\z\z, y\y\y, y\y\z, y\z\z, z\z\z}$. 
```py
for lx in range(lmax, 0-1, -1):
    for ly in range(lmax-lx, 0-1, -1):
        lz = lmax - lx - ly
        sx = lx * 'x'
        sy = ly * 'y'
        sz = lz * 'z'
        print(f"({lx}, {ly}, {lz})    ({sx}{sy}{sz})")
```



=== 归一化问题
笛卡尔基函数的归一化系数: @obarasaika
$
& phi.alt(x,y,z\; alpha, l_x, l_y, l_z, bold(R)) = N (x-R_x)^(l_x) (y-R_y)^(l_y) (z-R_z)^(l_z) "Exp"{-alpha |bold(r) - bold(R)|^2} \

& N = N(alpha, l_x, l_y, l_z) = ((2 alpha)/pi)^(3/4) (4 alpha)^((l_x+l_y+l_z)/2) 1/sqrt((2l_x-1)!! space (2l_y-1)!! space (2l_z-1)!!) \
& (-1)!! = 1 \ 
$<os>

纯基函数的归一化系数: 
$
& phi.alt(r, theta, phi) = N Y_l^m r^l e^(-alpha r^2) \ 
& N = sqrt((2 (2 alpha)^(l + 3/2)) / (Gamma(l + 3/2)))
$

*基函数既不正交，也不归一！* 

位于不同中心的基函数不正交，这也是为什么会存在重叠矩阵 $S_(i\j)$ 的原因. 然而更反直觉但数学上允许的是，基函数本身不需要归一，只要这些函数满足线性无关的条件，就是数学上允许的基函数. 不同的计算软件有不同的归一化规定，笛卡尔基函数的每一组分(primitive function)除了基组文件规定的系数外，还存在：
1. 每一组分都乘上各自的归一化系数(@os). 即所有组分都是归一化的（但是缩并后的基函数不保证归一化）. 
2. 以轴对齐组分（axis-aligned, i.e. xx, yy, xxx, ...）为基准，所有组分乘上相同的常数 $N = ((2 alpha)/(pi)) (4 alpha)^((l_"tot")/(2)) sqrt(1/(2 l_"tot" - 1)!!)$. 即只有轴对齐的组分是归一化的. 
3. 计算自重叠积分后，强制每一条基函数归一化. 

*重叠矩阵的对角元素，不保证全部为1!*

最简单的情况，是根本不管归一化，直接进行计算. 但是提前对基函数/组分进行归一化操作有其好处，其一是确保数值稳定性，防止出现极小/极大值（`1e-49`, `1e+17`, etc.. 重叠矩阵条件数更好）；其二是归一化的基函数可解释性更好. 但是无论采用何种归一化约定（甚至是根本不归一化），*确保在整个计算过程中使用相同的约定*. 特别是需要基函数在空间网格点上的函数值的情形（例如数值计算交换相关势）. 

归一化约定不影响物理可观测量(observable. 或期望值 expectation value)，但是会影响算符矩阵的构建. 假设 $diracc(phi.alt_mu, phi.alt_mu) equiv 1$, 那么未归一化的基函数是归一化基函数乘上缩放系数：
$
& phi.alt_mu^prime = a_mu phi.alt_mu \ 
& S_(mu nu)^prime = diracc(phi.alt_mu^prime, phi.alt_mu^prime) = a_mu a_nu S_(mu nu) \ 
& "let " C_(mu i)^prime = C_(mu i) / a_mu \ 
& => f_i = sum_mu C_(mu i) phi.alt_mu = sum_mu C_(mu,i)/a_mu  a_mu phi.alt_mu = sum_mu C_(mu,i)^prime phi.alt_mu^prime \ 
$

即波函数不受影响. 

但是算符的矩阵形式发生改变，因为算符的矩阵表示取决于基向量的选择：
$
& h_(mu nu)^prime = a_mu a_nu h_(mu nu) \ 
& (mu nu|lambda sigma)^prime = a_mu a_nu a_lambda a_sigma (mu nu|lambda sigma) \ 
$

$
& S_(mu nu)^prime = a_mu a_nu S_(mu nu) \ 
& S^prime = int mat(
              a_1,    ,          ,     ; 
                 , a_2,          ,     ; 
                 ,    , dots.down,     ; 
                 ,    ,          , a_N ;
             ) vec(phi.alt_1, phi.alt_2, dots.v, phi.alt_N) mat(phi.alt_1, phi.alt_2, dots.c, phi.alt_N) 
mat(
  a_1,    ,          ,     ; 
     , a_2,          ,     ; 
     ,    , dots.down,     ; 
     ,    ,          , a_N ;
) d^3 bvec(r) \ 
& = A S A \ 
$

$
& C_(mu i)^prime = C_(mu i)/a_mu \ 
& C^prime = mat(
  1/a_1,    ,          ,       ; 
     , 1/a_2,          ,       ; 
     ,      , dots.down,       ; 
     ,      ,          , 1/a_N ;
) 
mat(
  C_11  , C_12  , dots.c   , C_(1N) ; 
  C_21  , C_22  , dots.c   , C_(2N) ; 
  dots.v, dots.v, dots.down, dots.v ; 
  C_(N 1)  , C_(N 2)  , dots.c   , C_(N N) ; 
) = A^(-1) C \ 
$

系数矩阵随之变化. 

Kohn-Sham方程形式不变，且特征值不变：
$
& F^prime C^prime = A F A dot.c A^(-1) C = A F C equiv A S C epsilon \ 
& S^prime C^prime epsilon = A S C epsilon \ 
& => F^prime C^prime = S^prime C^prime epsilon \ 
$

密度矩阵变换后为($n_i$为占据数)：
$
& P_(mu nu) = sum_i^"occ" n_i C_(mu i) C_(nu i) = a_mu a_nu sum_i n_i C_(mu i)^prime C_(nu i)^prime = a_mu a_nu P_(mu nu)^prime \ 
& P_(mu nu)^prime = 1/a_(mu) P_(mu nu) 1/a_(nu) \ 
& => P^prime = A^(-1) P A^(-1)
$

期望值保持不变（trace invariant）: 
$
& dirac(hat(O)) = "Tr"{P^prime O^prime} = "Tr"{A^(-1) P A^(-1) A O A} <=> "Tr"{O A A^(-1) P} equiv "Tr"{P O}
$

波函数归一性（拉格朗日乘子法的约束条件）保持不变：
$
& diracc(f_i, f_i) = sum_mu sum_nu diracc(C_(mu i) phi.alt_mu, C_(nu i) phi.alt_nu) \ 
& = sum_(mu nu) C_(mu i) C_(nu i) S_(mu nu) equiv 1 \ 
& <=> C^T S C equiv I \ 
& (C^prime)^T S^prime C^prime = (A^(-1) C)^T A S A A^(-1) C = C^T S C equiv I
$



#linebreak()
== 高斯求积
注：本节内容不是严格的数学定义. 

数值积分是将函数 $f(x)$ 在区间 $[a,b]$ 上的离散点处的函数值做加权平均作为积分近似值（基本思想为积分中值定理）: 
$
& a <= x_0 <= x_1 <= dots.c <= x_n <= b \ 
& int_a^b f(x) d\x approx sum_(i=0)^n A_i f(x_i) \  
$
其中 ${A_i}$ 为系数（权重），${x_i}$ 为求积节点. 

最简单的梯形法使用的是等分节点，各点处的权重相同，需要非常多的离散点才能达到满意的数值精度. 为了用更少的节点获得更高的精度，节点的选择和权重的确定很重要. 

如果存在节点 $x_i in [a,b]$ 以及系数 ${A_i}$ 使得 $int_a^b f(x) d\x approx sum_(i=0)^n A_i f_i$ 具有 $2n+1$ 次代数精度，节点 ${x_i}$ 称为*高斯点*，$A_i$ 称为*高斯系数*, 求积公式称为*高斯型求积公式*. 

通常高斯点是各种正交多项式在各自正交区间上的零点，系数（权重）有专门的公式进行计算，且进行数值求积时需要考虑各种正交多项式的权重函数： 
$
& int_a^b f(x) omega(x) d\x approx sum_(i=1)^n f(x_i) w_i \ 
$




=== 正交多项式
如果一族多项式中的任意两个不同的多项式的带权内积为零，那么这族多项式是正交多项式: 
$
& chevron.l p_i dot.c p_j chevron.r = int_a^b omega(x) p_i (x) p_j (x) d\x = N_(i j) delta_(i j) \ 
$

三角函数 ${1, cos\x, sin\x, cos\2x, sin\2x, dots.c}$ 在区间 $[-pi, pi]$ 上构成正交函数族（不是多项式）. 

- 如果 ${phi.alt_k}_(k=0)^infinity$ 在区间 $[a,b]$ 上带权 $omega$ 正交，$H_n$ 是所有次数不超过 $n$ 的多项式组成的线性空间，那么 ${phi.alt_i | i=0, 1, dots.c, n}$ 构成 $H_n$ 的一组基. 即对所有 $p(x) in H_n$, $p(x) = sum_(j=0)^n c_j phi.alt_j$. 

- 正交多项式满足三项递推关系：
$
& phi.alt_(-1) = 0 \ 
& phi.alt_0 = 1 \ 
& phi.alt_(n+1) = (x - alpha_n) phi.alt_n - beta_n phi.alt_(n-1), space n = 0, 1, 2, dots.c \ 
& alpha_n = ((x phi.alt_n, phi.alt_n))/((phi.alt_n, phi.alt_n)) \ 
& beta_n = ((phi.alt_n, phi.alt_n))/((phi.alt_(n-1), phi.alt_(n-1)))
$

三项递推关系是高效计算正交多项式的手段. 

- 如果 ${phi.alt_n}_0^infinity$ 在区间 $[a,b]$ 带权函数 $omega$ 正交, 那么 $phi.alt_n$ 在区间 $[a,b]$ 内有 $n$ 个不同零点. 




=== 常用正交多项式
- 勒让德求积
积分区间 $[-1,1]$, 权重函数 $omega(x) = 1$, $P_0=1$, $P_1 = x$, $P_n = 1/(2^n n!) (d^n)/(d x^n) (x^2-1)^n$. 
1. $
& int_(-1)^1 P_n P_m d\x = cases(0\, space m != n, 2/(2n+1)\, space m = n) \ 
$

2. $
& P_n (-x) = (-1)^n P_n (x) \ 
$

3. $
& (n+1) P_(n+1) = (2n+1) x P_n - n P_(n-1), space n=1,2,dots.c \ 
$

4. $P_n (x)$ 在区间 $[-1,1]$ 内有 $n$ 个不同的零点. 

求积公式：$int_(-1)^1 f(x) d\x approx sum_(i=1)^n f(x_i) w_i$




- 拉瓜尔求积
积分区间 $[0, infinity)$, 权重函数 $omega(x) = e^(-x)$, $L_n = e^x (d^n)/(d x^n) (x^n e^(-x))$. 

1. $
& int_0^infinity e^(-x) L_n L_m d\x = cases(0\, space m!=n , (n!)^2\, space m=n) \ 
$

2. $
& L_0 = 1 \ 
& L_1 = 1 - x \ 
& L_(n+1) = (1+2n-x) L_n - n^2 L_(n-1), space n=1,2,dots.c 
$

求积公式：$int_(0)^infinity f(x) e^(-x) d\x approx sum_(i=1)^n f(x_i) w_i$




- 厄米求积
区间 $(-infinity, infinity)$, 权重函数 $omega(x) = e^(-x^2)$, $H_n = (-1)^n e^(x^2) (d^n)/(d x^n) (e^(-x^2))$. 

1. $
& int_(-infinity)^infinity e^(-x^2) H_n H_m = cases(0\, space m!=n , 2^n n! sqrt(pi)\, space m = n) \ 
$

2. $
& H_0 = 1 \ 
& H_1 = 2 x \ 
& H_(n+1) = 2x H_n - 2n H_(n-1), space n=1,2,dots.c \ 
$

求积公式：$int_(-infinity)^infinity f(x) e^(-x^2) d\x approx sum_(i=1)^n f(x_i) w_i$



- 切比雪夫求积(第一类)
区间 $[-1,1]$, 权重函数 $omega(x) = 1/sqrt(1-x^2)$

1. $
& T_n (x) = cos(n arccos\x) space |x| <= 1 \ 
& T_0 = 1 \ 
& T_1 = x \ 
& T_(n+1) = 2x T_n - T_(n-1), space n=1,2,dots.c \ 
$

2. $
& int_(-1)^1 (T_m T_n)/(sqrt(1-x^2)) d\x = cases(0\, space n!=m , pi/2\, space n=m!=0 , pi\, n=m=0)
$

3. $
& T_n (-x) = (-1)^n T_n (x) \ 
$

4. $T_n$ 的 $n$ 个零点
$
& x_k = cos((2k-1)/(2n) pi) space k=1,2,dots.c, n \ 
$

求积公式：$int_(-1)^1 f(x) / sqrt(1-x^2) d\x approx sum_(i=1)^n f(x_i) w_i$





- 切比雪夫求积(第二类)
区间 $[-1,1]$, 权重函数 $omega(x) = sqrt(1-x^2)$

1. $
& U_n = (sin[(n+1) arccos\x])/(sqrt(1-x^2)) \ 
& U_0=1 \ 
& U_1 = 2x \ 
& U_(n+1) = 2x U_n - U_(n-1), space n=1,2,dots.c \ 
$

2. $
& int_(-1)^1 U_n U_m sqrt(1-x^2) = cases(0\, space m!=n , pi/2\, space m=n) \ 
$

求积公式：$int_(-1)^1 f(x)  sqrt(1-x^2) d\x approx sum_(i=1)^n f(x_i) w_i$


*总结*

#align(center)[
#table(
  columns: 3,
  stroke: black,
  align: left, 
  [类别], [区间], [权重函数], 
  [Legendre], [$[-1,1]$],        [$omega(x)=1$], 
  [Laguerre], [$[0, infinity)$], [$omega(x) = e^(-x)$], 
  [Hermite], [$(-infinity, infinity)$], [$omega(x) = e^(-x^2)$], 
  [Chebyshev I], [$(-1,1)$], [$omega(x) = 1/sqrt(1-x^2)$], 
  [Chebyshev II], [$[-1,1]$], [$omega(x) = sqrt(1-x^2)$]
)
]

*例子*

- 计算 $int_(-1)^1 1+x+x^2+x^3+x^4+x^5 \dx$


```py
f = lambda x: 1+x+x**2+x**3+x**4+x**5

# Gauss-Legendre quadraure
import numpy.polynomial.legendre as L
# order 5 needs only 3-order roots, as 2x3-1=5, but here we use 5 roots
res = L.leggauss(5) # tuple, res[0] is roots, res[1] is weights
quad = np.sum(f(res[0]) * res[1]) # 3.0666666666666673, exact value 46/15

# Gauss-ChebyshevII quadrature
n = 10
z = np.arange(1, n+1)
x = np.cos(z / (n+1) * np.pi)
w = np.pi / (n+1) * np.sin(z/(n+1) * np.pi)**2
quad = np.sum(f(x) / np.sqrt(1-x**2) * w) # big error!
```

使用高斯-切比雪夫求积时需要注意，如果目标积分为 $int_(-1)^1 f(x) sqrt(1-x^2) \dx$, 那么对于 $2n-1$ 次函数 $f(x)$, 仅需 $n$ 个节点即可准确计算积分. 但是对于 $int f(x) \dx$ 而言，使用切比雪夫求积需要先除以权重函数：$f(x) -> f(x) / sqrt(1-x^2)$, 而此时的目标函数 $(f(x))/(sqrt(1-x^2))$ 不是有限次数多项式！仅靠几个节点得到的积分结果误差很大. 



=== 变量变换
数学物理中的很多函数是定义在无界/半无界区间上的平方可积函数 $f in LL^2$. 然而计算这些函数的积分时通常使用的不是拉瓜尔或厄米积分. 

以厄米求积为例，由于其积分区间是整个实数域，其零点的分布为 $-infinity < x_i < infinity$. 如果某函数恰好可以写成一个多项式函数乘以高斯函数那么使用厄米求积是非常自然的选择：$F(x) = sum_j c_j x^j e^(-x^2)$. 但是如果不能写成这种形式，求积时需要除以高斯权重函数：$int f(x) d\x = int f(x)/e^(-x^2) e^(-x^2) d\x$. 由于厄米多项式的零点分布在整个实数域，分母 $e^(-x^2)$ 可能相当小并导致 $f(x)\/e^(-x^2)$ 数值不稳定. 另外，当我们需要通过除以高斯权重函数“凑出”要求的求积公式时，积分核 $f(x)\/e^(-x^2)$毫无疑问已经不是多项式了(无穷级数)，有限项求和只能得到近似结果，厄米求积已经失去优越性. 增加项数可以提高精度，但是又会遇到数值稳定性问题. 

但是，仍然可以使用(而且效果相当好). 例如：$f(x) = 1/(x^2+1) e^(-x^2) + (x-1)^3 e^(-5 (x+1)^2)$

```py
import numpy.polynomial.hermite as H

# exact -43 sqrt(pi/5) / 5 + e pi Erfc[1] = -5.4736295302356037725
f = lambda x: (1/(x**2+1)) * np.exp(-x**2) + (x-1)**3 * np.exp(-5 * (x+1)**2)
res = H.hermgauss(128)
quad = np.sum(f(res[0]) / np.exp(-res[0]**2) * res[1]) # -5.47362953023580089962, AbsErr = 1e-13
# N.B. np.exp(-res[0]**2) min=1e-102
```

另一种数值计算反常积分的方法是使用勒让德求积. 勒让德求积的优势之一是其具有最简单的权重函数（$omega(x)=1$），另外由于其积分区间为 $[-1,1]$, 其零点也被限制在该范围内，数值上不会出现病态小或病态大值. 同样由于其有效区间有界，需要通过变量变换将有界区间映射到（伪）无界区间. 

对于半无界区间平方可积函数的反常积分，数值求积无法处理无穷边界，但是考虑到 $f|_(x=infinity)=0$, 反常积分的计算可以近似成上限很大的定积分的计算：
```
N[Integrate[Exp[-x^2], {x, 0, Infinity}] - Integrate[Exp[-x^2], {x, 0, 10^4}]]; 
Out[]= 0.
```

为了将 $(-1,1)$ 映射到 $(0, "Big")$, 定义（一种可能的）映射函数：
$
& x in (-1,1) mapsto T(x) in (0, "Big") \ 
& "left" = 1e-3 \ 
& "right" = 1e+4 \ 
& alpha = ln("right"/"left") \ 
& T(x) = "left" (e^(alpha (x+1)/2) - 1) \ 
& T^prime (x) = alpha/2 "left" e^(alpha (x+1)/2) \ 
$

那么积分 $int_0^infinity F(x) d\x$ 的计算可以近似为：
$
& int_0^infinity F(x) d\x approx int_0^"Big" F[T(x)] d T(x) <=> int_(-1)^1 F[T(x)] T^prime (x) d\x \ 
& approx sum_(i=1)^n F[T(x_i)] underbrace((T^prime (x_i)), "Jacobian") w_i
$

其中 $T^prime (x)$ 是由于变量变换产生的雅可比行列式（1D scalar）. 

*例子* 

- 通过勒让德求积+变量映射计算 $int_0^infinity 1/(x^2+1) e^(-x^2) + (x-1)^3 e^(-5 (x+1)^2)$

```py
left = 1e-3
right = 1e+4
alpha = np.log(right/left)
f = lambda x: (1/(x**2+1)) * np.exp(-x**2) + (x-1)**3 * np.exp(-5 * (x+1)**2) # integrand
T = lambda x: left * (np.exp(alpha * (x+1)/2) - 1.) # mapping function
jacobian = lambda x: alpha/2. * left * np.exp((x+1)/2 * alpha) # Jacobian det

res = L.leggauss(128)
x_mapped = T(res[0])
fnvals = f(x_mapped)
quad = np.sum(fnvals * res[1] * jacobian(res[0])) # 0.67116241937370002546, AbsErr = 1e-12
# exact 1/50 (36/e^2 + 25 e pi Erfc[1] - 43 sqrt{5 pi} Erfc[sqrt{5}]) = 0.67116241937195611168
```


同理，相同的手段可以套用在切比雪夫求积(第二类). 与勒让德求积相比，切比雪夫求积需要考虑权重函数，但是其优势在于节点和权重的计算很简单，*且节点的分布可以从等差数列得到*: 
$
& z = 1,2, dots.c, N space "等距网格" \ 
& x_i = cos(z_i/(N+1) pi) in (1,-1) space "节点" \ 
& w_i = pi/(N+1) sin^2(z_i/(N+1) pi) space "权重" \ 
$

使用Becke建议的映射函数：$x mapsto r=r(x) = r_m (1+x)/(1-x) in ("Big", 0) $. 其中 $r_m$ 是元素的 Bragg-Slater 半径@becke1988@slaterradii. 一方面使得切比雪夫求积可以应用于半无界区间的积分，另一方面该映射函数在原子核处密集，远离原子核处稀疏，符合电子密度的分布特点.

#figure(
  image("pics/tf.svg", width: 50%), 
  caption: "Becke mapping function.", 
  supplement: "Figure"
)

$z in [1,N]$, $x in (1,-1)$, $r in (infinity, 0)$. 最终变量${r_i|i=1,2,dots.c,N}$是最初变量${z_i}$的函数，所有以$r$为变量的函数经过链式法则都是变量$z$的函数. 这带来一个好处：如果需要通过有限差分计算$f(r)$对$r$的导数，可以通过链式法则将问题转变为计算$f(r)$对$z$的导数，而$Delta z=1$. 如果使用勒让德求积，原始变量$z$不是均匀网格，有限差分在非均匀网格上不好操作. 

切比雪夫求积需要除以权重函数. 再考虑上变量变换带来的雅可比行列式，最终可以把这些因素全部包含进权重中：
$
& (d\r)/(d\x)|_(x_i) = 2 r_m 1/(1-x_i)^2 \ 
& int_0^infinity f(r) r^2 d\r approx sum_(i=1)^N (f[r(x_i)] r_i^2) / sqrt(1-x_i^2) pi/(N+1) sin^2(z_i/(N+1) pi) 2 r_m 1/(1-x_i)^2 \ 
& = sum_(i=1)^N f(r_i) r_i^2 pi/(N+1) sin^2(z_i/(N+1) pi) (2 r_m) / (sqrt(1-x_i^2) (1-x_i)^2) \ 
& => w_i = pi/(N+1) sin^2(z_i/(N+1) pi) (2 r_m) / (sqrt(1-x_i^2) (1-x_i)^2) \ 
$

注意节点是倒序排列：${x_i} in (1,-1)$. 求积时节点倒序，映射 ${r_i} in ("Big", 0)$ 也是倒序，权重 ${w_i}$ 需要与之对应. 

*例子*
- $int_0^infinity 1/(r^2+1) e^(-r^2) + (r-1)^3 e^(-5 (r+1)^2) r^2 d\r$. 对一个三维函数的径向部分求积（即1D问题）
```py
def gaussCheby2(norder:int=32, rm:float=1.) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z = np.arange(1, norder+1)
    x = np.cos(z / (norder + 1) * np.pi) # z ∈ [1, norder] ↦ x ∈ (1, -1)
    r = rm * (1 + x) / (1 - x)           # x ∈ (1, -1) ↦ r ∈ (∞, 0)
    w = np.pi / (norder + 1) * np.sin(z / (norder + 1) * np.pi)**2 / np.sqrt(1 - x**2) * 2 * rm / (1 - x)**2
    return x, r, w

f = lambda x: (1/(x**2+1)) * np.exp(-x**2) + (x-1)**3 * np.exp(-5 * (x+1)**2)

x, r, w = gaussCheby2(128)
quad = np.sum(f(r) * r**2 * w) # 0.21457600129729723; AbsErr = 1e-14
# exact 1/250(277/e^5 - 125e pi Erfc[1] + sqrt{pi} (125 - 301 Erfc[sqrt{5}])) = 0.2145760012973228
```



== 球谐函数与数值求积
=== 实球谐函数<rsh> 

球谐函数是角动量算符($hat(L)^2$)及其 $z$-分量 ($hat(L)_z$)的特征函数: 
$
& hat(L)^2 Y_l^m (theta, phi) = l(l+1) hbar^2 Y_l^m (theta, phi) \
& hat(L)_z Y_l^m (theta, phi) = m hbar Y_l^m (theta, phi) \ 
$

球谐函数通常是复的：
$
& Y_l^m (theta, phi.alt) = (-1)^m N_(l m) P_l^m (cos theta) e^(i m phi.alt) \ 
& N_(l m) = sqrt((2l+1)/(4pi) ((l-m)!)/((l+m)!)) \ 
$

$(-1)^m$ is Condon-Shortley phase and $P_l^m$ is associated-Legendre polynomial. 

且在单位球面上正交归一：
$
& int_0^pi int_0^(2pi) Y_l^m space Y_(l^prime)^(m^prime) sin theta d theta d phi = delta_(l l^prime) delta_(m m^prime) = chevron.l Y_l^(m)mid(|)Y_(l^prime)^(m^prime) chevron.r_Omega \ 
$

为了实际计算的便利性, 且微分方程解的线性叠加仍是原微分方程的解, 在实际使用中通常对球谐函数进行线性组合以消除虚部, 得到的就是所谓的实球谐函数(real spherical harmonics). 无论是复球谐函数还是实球谐函数, 他们通常都是角度$(theta,phi)$的函数: $Y_l^m (Omega) = Y_l^m (theta,phi)$. 角度变量可以通过反三角函数转化为直角坐标, 因此$Y_l^m (theta,phi) equiv Y_l^m (x,y,z)$.

实球谐函数的构造不唯一，但是必须满足正交归一性：$& int_0^(2pi) int_0^pi (Y_(l^prime)^(m^prime))^ast Y_l^m sin theta d theta d phi.alt = delta_(l l^prime) delta_(m m^prime)$. 

下面是几种可能的构造方法：
1. $
& Y_l^m = cases(
  (-1)^m/sqrt(2) {Y_l^m + (Y_l^m)^dagger} space space m > 0 ,
  ,
  Y_l^0 space space m = 0, 
  ,
  (-1)^m/(i sqrt(2)) {Y_l^(|m|) - (Y_l^(|m|))^dagger} space space m < 0
) \ 
$

2. $
& Y_l^m = cases(
  sqrt(2) (-1)^m Im[Y_l^(|m|)] space space m < 0, 
  , 
  Y_l^0 space space m = 0, 
  , 
  sqrt(2) (-1)^m Re[Y_l^m] space space m > 0, 
) \ 
$

可能的代码实现：
1. ```py
# Copied from PyDFT: [https://github.com/ifilot/pydft]. 
# no Condon-Shortley phase
def rsh(l, m, theta, phi):
    if m < 0:
        val = np.sqrt(2) * np.imag(sph_harm_y(l, np.abs(m), phi, theta)) # ifilot swaps theta and phi
    elif m > 0:
        val = np.sqrt(2) * np.real(sph_harm_y(l, m, phi, theta))
    else:
        val = np.real(sph_harm_y(l, m, phi, theta))
    
    return val
```

2. ```py
# https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/ 
# https://en.wikipedia.org/wiki/Spherical_harmonics#Real_form
def rsh(l:int, m:int, theta:float, phi:float) -> float:
    if abs(m) > l:
        raise ValueError('magnetic quantum number m should not be larger than angular qunautm number l')
    
    phase = -1 if abs(m) % 2 else 1

    if m > 0:
        return phase * np.sqrt(2.) * sph_harm_y(l, m, theta, phi).real
    elif m < 0:
        return phase * np.sqrt(2.) * sph_harm_y(l, -m, theta, phi).imag
    elif m == 0:
        return sph_harm_y(l, 0, theta, phi).real # a+0j, forcefully remove null imaginary part
```

3. ```py
# https://docs.abinit.org/theory/spherical_harmonics/
def rsh2(l:int, m:int, theta:float, phi:float) -> float:
    if abs(m) > l:
        raise ValueError('magnetic quantum number m should not be larger than angular qunautm number l')
    
    phase = 1 if (abs(m) % 2 == 0) else -1 # (-1)^m

    if m > 0:
        return phase / np.sqrt(2.) * (sph_harm_y(l, m, theta, phi) + sph_harm_y(l, m, theta, phi).conjugate())
    elif m < 0:
        return phase / np.sqrt(2.) / 1j * (sph_harm_y(l, -m, theta, phi) - sph_harm_y(l, -m, theta, phi).conjugate())
    elif m == 0:
        return sph_harm_y(l, 0, theta, phi)
```


4. C++标准库中定义了勒让德多项式，可以拿来直接用：

  ```cpp
double rsh(const int l, const int m, const double theta, const double phi) noexcept
{
    assert((l >= 0) && (std::abs(m) <= l));
    const    int        mm = (m < 0 ? -m : m);
    const double     phase = (mm & 1) ? -1. : 1.;
    const double         P = std::sph_legendre(l, mm, theta);
    constexpr double root2 = 1.4142135623730951454746218587388284504414;

    if (m == 0) return P;
    if (m > 0)  return phase * root2 * P * std::cos(mm * phi);
    return phase * root2 * P * std::sin(mm * phi);
}
```



=== 列别杰夫求积
相较于径向求积有多种正交多项式可供选择，角度（球面）求积最常用的是列别杰夫求积@lebedev1975values. 

- *例子*@beentjes2015quadrature

计算在单位球面上的积分: $int_Omega f(x,y,z) d Omega = int_Omega 1+x+y^2+x^2 y + x^4 + y^5 + x^2 y^2 z^2$

```py
import numpy as np
from scipy.integrate import lebedev_rule

f = lambda x,y,z: 1 + x + y**2 + x**2 * y + x**4 + y**5 + x**2 * y**2 * z**2

res = lebedev_rule(17) # 110 angular pts
quad = np.sum(f(*res[0]) * res[1]) # 19.388114662154152
# exact 216pi/35 = 19.388114662154152
```

通过 Mathematica 验证一下（单位球面 $r=1$）: 
```txt
f[x_,y_,z_]:=1+x+y^2+x^2 y+x^4+y^5+x^2 y^2 z^2;
Integrate[(f[x,y,z]/.{x->Sin[\[Theta]] Cos[\[Phi]],y->Sin[\[Theta]] Sin[\[Phi]],z->Cos[\[Theta]]}) Sin[\[Theta]],{\[Theta],0,Pi},{\[Phi],0,2 Pi}]
Out[]=216Pi/35
```





== 泊松方程数值解
=== 泛函一阶导数

类比于函数的导数 $f^prime (x) = (d f(x))/(d\x)$, $d\f = f^prime d\x$. 泛函的一阶导数为：
$
  & F[f(x)] \
  & F'[f(x)] = (delta F[f(x)]) / (delta f(x)) \
  & delta F = F[f(x) + delta f(x)] - F[f(x)] = underbrace(integral ((delta F) / (delta f)) delta f \dx, "contributions from all x") + cal(O)(delta^2 f) \
$

与函数导数不同，泛函导数不存在公式. 导数结果取决于泛函的具体形式. 提取积分中 $delta f$ 的系数即为泛函导数. 

*例子*

- 计算泛函 $F[n(x)] = integral_0^1 n^2 \dx$ 的一阶导数: 
$
  & delta F = integral_0^1 (n + delta n)^2 - n^2 \dx \
  & delta F = integral_0^1 2n delta n \dx \
  & => (delta F) / (delta n) = 2n \
$

*例子*

- 积分型泛函的一阶导数：

一般认为 $delta n$ 在端点（或边界）处固定: $delta n = 0$，因此 $(partial f) / (partial n') delta n mid(|)_"边界" = 0$ 
$
  & F[n] = integral f(n, n', x) \dx \
  & F[n+delta n] = integral f(n+delta n, d / (\dx) (n + delta n), x) \dx \
  & \
  & "泰勒展开到一次项" \
  
  & => F[n+delta n] = int f(n, n^prime, x) + {(partial f)/(partial n) delta n + (partial f)/(partial n^prime) (delta n)^prime} space \dx \
  
  & F[n+delta n] - F[n] = delta F = int (partial f)/(partial n) delta n + (partial f)/(partial n^prime) (delta n)^prime space \dx \
  & \

  & "对第二项进行分部积分" \
  & => delta F = int (partial f)/(partial n) delta n space \dx + int (partial f)/(partial n^prime) d (delta n) = int (partial f)/(partial n) delta n space \dx + (partial f)/(partial n^prime) delta n |_("边界") - int delta n space d ((partial f)/(partial n^prime)) \

  & => delta F = int (partial f)/(partial n) delta n space \dx - int delta n space d/(\dx) ((partial f)/(partial n^prime)) space \dx \

  & delta F = int (delta F[n])/(delta n) space delta n space \dx = int {(partial f)/(partial n) - d/(\dx) ((partial f)/(partial n^prime))} space delta n space \dx \

  & therefore (delta F[n])/(delta n) = (partial f)/(partial n) - d/(\dx) ((partial f)/(partial n^prime))
$

这是欧拉-拉格朗日方程的形式. 



库伦势能 $E[n(bvec(r))]_"Coul"$ 对电子密度的泛函导数得到库伦势 $v[n(bvec(r))]_"Coul" = (delta E[n])/(delta n)$: 

$
& E_"Ha" [n] = 1/2 int int (n(bvec(r)) n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r) \d^3 bvec(r)^prime \
& delta E_"Ha" [n] = E_"Ha" [n+delta n] - E_"Ha" [n] \

& E_"Ha" [n+delta n] = 1/2 int int ({n(bvec(r)) + delta n(bvec(r))} {n(bvec(r)^prime) + delta n(bvec(r)^prime)})/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r) \d^3 bvec(r)^prime \

& approx 1/2 int int (n(bvec(r)) n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) + underbrace({(n(bvec(r)) delta n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) + (n(bvec(r)^prime) delta n(bvec(r)))/(|bvec(r)-bvec(r)^prime|)}, bvec(r) <--> bvec(r)^prime) \d^3 bvec(r) \d^3 bvec(r)^prime \

& = E_"Ha" [n] + int int (n(bvec(r)^prime) delta n(bvec(r)))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r) \d^3 bvec(r)^prime \

& => delta E_"Ha" [n] = int int (n(bvec(r)^prime) delta n(bvec(r)))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r) \d^3 bvec(r)^prime \

& <=> int {int (n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r)^prime} delta n(bvec(r)) \d^3 bvec(r) \

& therefore delta E_"Ha" [n(bvec(r))] = int (delta E_"Ha" [n(bvec(r))])/(delta n(bvec(r))) delta n(bvec(r)) \d^3 bvec(r) \

& => v_"Ha" (bvec(r)) = (delta E_"Ha" [n(bvec(r))])/(delta n(bvec(r))) = int (n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r)^prime \
$<hartreepotential>


=== 泊松方程

重要等式：
$ 
& nabla^2_bvec(r) 1/(|bvec(r)-bvec(r)^prime|) = -4 pi delta(bvec(r)-bvec(r)^prime) \
$<green>

其中 $1/(|bvec(r) - bvec(r)^prime|)$ 是格林函数(Green's function). 角标 $nabla^2_bvec(r)$ 指的是对变量 $bvec(r)$ 施加拉普拉斯算符. 

@hartreepotential 本质上是库伦核函数对电子密度的卷积：

$
& (f * g)(t) = int f(tau) g(-tau+t) d tau \
& f = n(bvec(r)) \
& g = 1/(|bvec(r)|) \ 
& => (f * g)(bvec(r)) = int n(bvec(r^prime)) 1/(|-bvec(r^prime) + bvec(r)|) \d^3 bvec(r^prime) = int (n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r)^prime \
$

将@green 带入@hartreepotential 可以得到：
$
& nabla_bvec(r)^2 v_"Ha" (bvec(r)) = nabla^2 int (n(bvec(r^prime)))/(|bvec(r-r^prime)|) d^3 bvec(r^prime) = int n(bvec(r^prime)) nabla^2 1/(|bvec(r)-bvec(r)^prime|) d^3 bvec(r^prime) \ 

& = int n(bvec(r^prime)) dot.c -4 pi delta(bvec(r)-bvec(r)^prime) d^3 bvec(r^prime) = -4 pi n(bvec(r)) \
$<hartreepotential2>

由此我们得到*库伦势的积分与微分表达式*，二者等价：
$
& v_"Ha" (bvec(r)) = int (n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r)^prime \

& nabla^2 v_"Ha" (bvec(r)) = -4 pi n(bvec(r)) \
$<poisson>

*其中库伦势的微分形式即为“泊松方程”. 求解泊松方程就是在给定电子密度时计算库伦势. *






=== 球谐展开

由于（实）球谐函数在单位球面上正交归一 ($theta in (0, pi)$, $phi.alt in (0, 2 pi)$, $r equiv 1$)，球谐函数在单位球面上构成完备集，单位球面上的平方可积函数都可以表为球谐函数的线性组合: 
$
f(theta, phi.alt) = sum_(l=0)^infinity sum_(m=-l)^l a_(l m) Y_l^m (theta, phi), space space l in bb(N) \
$

球坐标下积函数是球谐函数和高斯函数线性组合的乘积 (N.B. primitive gaussian normalized, contracted gaussian normalization not guaranteed!)：
$
phi.alt(r,theta,phi) = Y_l^m (theta,phi) |bvec(r)-bvec(R)|^l sum_i c_i n_i e^(-alpha_i |bvec(r)-bvec(R)|^2)
$ 

数学上与笛卡尔基函数等价：
$
phi.alt(x,y,z) = (x-R_x)^(l_x) (y-R_y)^(l_y) (z-R_z)^(l_z) sum_i c_i n_i e^(-alpha_i |bvec(r)-bvec(R)|^2) \ 
$

基函数径向部分和角度部分耦合在一起无法显式分离, 但是可以利用球谐展开将笛卡尔基函数的径向部分剥离出来：
$
& |phi.alt(r,theta,phi)|^2 equiv |phi.alt(x,y,z)|^2 = sum_(l m) rho_(l m) (r) Y_l^m (theta,phi) \ 
& rho_(l m) (r) = int Y_l^m (theta,phi)^ast |phi.alt(r,theta,phi)|^2 d Omega \ 
& rho_(l m) (r_i) approx sum_j^M |phi.alt(r_i, Omega_j)|^2 Y_l^m (Omega_j)^ast w_j
$<sphericalexpansion>

@sphericalexpansion 中的$Omega_j$是单位球面上的采样点 (通过列别杰夫求积构建), $w_j$ 是求积权重. 


// *不需要根据原子中心位置进行平移或缩放*，原子网格在构建时已经平移缩放过了，这里出现的$Y_l^m$只是为了通过投影计算展开系数. 

@sphericalexpansion 需要展开到无穷项才使等式成立；除非被展开函数本身就是球谐函数或球谐函数的幂，此时展开项有限，且展开系数可以通过Clebsch-Gordan系数解析确定. 对于一般情况, Becke建议 $l$ 取列别杰夫求积阶数的一半. 例如 `nleb=29 (302 angular pts)`, `lmax=14`, `nlm=225` @becke19882. 

// 这里以氧原子2pz原子轨道为例，角动量最高取到$l=4$，通过球谐展开计算其重叠积分. 

// 首先构建以原子为中心的网格, 径向网格数`nrad`, 在每一径向格点$r_i$处有角度网格数`nang`. 确定角度网格格点数后可以预先计算$l_"max"=4$各实球谐函数在单位球面上的函数值（注意，既不缩放，也不平移）：

// #figure(
//   image(
//     "./pics/demo1.svg", width: 40%
//   )
//   , caption: ""
// )

// 球谐函数在各$(l,m)$处的函数值记作矩阵$Y_(j,k)$, 每行代表一个单位球面上的格点在不同角动量情况下的函数值. 

// `ao_vals`则是被展开函数($phi.alt(x,y,z)$)的矩阵. $"val"_"i,j"$代表在径向格点$r_i$，角度格点$Omega_j$处的被展开函数函数值. 这里的角度格点以原子为中心经过了缩放和平移. 

// 球谐展开利用了球谐函数的正交归一性：
// $
// & rho_(i k) = int Y_(j k) |phi.alt(x,y,z)|^2 d Omega \ 
// & "rho"_"i,k" = "val"_"ij"^2 "@" ("Y"_"jk" "wleb")
// $

// $"rho"_"ik"$是关于径向网格${r_i}$的函数，每一行是某一径向格点在不同角动量处球谐函数的展开系数. 

// 展开系数$rho_(i,k)$确定后，原子轨道重叠积分可以通过下式计算:
// $
// & angle.l phi.alt|phi.alt angle.r = int sum_(l m) rho_(l m) (r) Y_l^m (theta,phi) r^2 \dr d Omega \ 
// & approx sum_i^"nrad" r_i^2 w_i sum_(l m) rho_(i,l m) sum_j^"nang" Y_(j, l m) w_j \ 
// & = sum_i sum_(l m) sum_j (r_i^2 w_i) rho_(i,l m) (Y_(j,l m) w_j) \ 
// $

// 当`nrad`=32, `nang`=110，绝对误差为1e-5. 

// 总结一下上述两种方法：
// 1. 物理量在网格上的函数值已知，在每一径向网格处对该径向格点处的角度网格尽心求积，遍历所有径向网格. 

// 2. 球谐展开剥离出函数的径向部分，然后对展开式进行求积. 

// 当然，这两种方法不是仅用来计算基函数的自重叠积分的. 这里只是用其作为例子. 


=== 通过球谐展开和有限差分数值求解泊松方程

对电子库伦势 $U(bvec(r))$ 和电子密度 $rho(bvec(r))$ 进行球谐展开:
$
& U(bvec(r)) = U(r,theta,phi) = sum_(\lm) (u_(\lm)(r))/r Y_l^m (theta, phi) \
& rho(bvec(r)) = rho(r,theta,phi) = sum_(\lm) rho_(\lm)(r) Y_l^m (theta,phi) \
$

带入@poisson

$
& nabla^2 = 1/r^2 partial/(partial r) (r^2 partial/(partial r)) - L^2/r^2 \ 
& nabla^2 sum_(\lm) u_(\lm)/r Y_l^m = -4pi sum_(\lm) rho_(\lm) Y_l^m \ 
& "lhs" = sum_(\lm) {1/r^2 partial/(partial r) (r^2 partial/(partial r)) - L^2/r^2} u/r Y \

& = sum_(\lm) 1/r^2 partial/(partial r) (r^2 partial/(partial r) u/r) Y - 1/r^2 u/r L^2 Y \

& = sum_(\lm) 1/r^2 partial/(partial r) (r^2 (u^prime r - u)/r^2) Y - u/r^3 l(l+1) Y \
& = sum_(\lm) 1/r^2 (u^(prime prime) r + u^prime - u^prime) Y - u/r^3 l(l+1) Y \
& => sum_(\lm) u_(\lm)^(prime prime)/r Y_l^m - u_(\lm)/r^3 l(l+1) Y_l^m = sum_(\lm) -4pi rho_(\lm)(r) Y_l^m \

& => d^2/(\dr^2) u_(\lm) - (l(l+1))/(r^2) u_(\lm) = -4pi r rho_(\lm) \ 
$

每确定一组 $(l,m)$ 得到一个二阶常微分方程：
$
& d^2/(d r^2) u_(\lm) (r) - (l(l+1))/r^2 u_(\lm) (r) = -4pi r rho_(\lm) (r) \
$

问题就从求解一个二阶偏微分方程变为求解一组二阶常微分方程. 

均匀网格 ${z_i|i=1,2,dots.c,N}$ 通过余弦变换转变为切比雪夫求积(第二类)的节点: ${x_i = cos (z_i)/(N+1) pi | i=1,2,dots.c,N} in (1,-1)$, 再通过变换映射到半无穷区间上: ${r_i = r_m (1+x_i)/(1-x_i) | i=1,2,dots.c, N} in ("Big",0)$. 

二阶微分方程需要至少两个边界条件，由 $(z_0, x_0, r_0)$ 以及 $(z_(N+1), x_(N+1), r_(N+1))$ 确定 (Dirichlet condition). 

根据链式法则，库伦势的径向部分是关于 $r$ 的标量函数: $u_(\lm) = u_(\lm) (r)$, 因此也是关于均匀网格 ${z_i}$ 的标量函数: $u_(\lm) = u_(\lm)[r(z)]$:

$
& (d^2 u(r))/(d r^2) = d/(\dr) (d u)/(d r) = d/(\dr) ((\du)/(\dz) (\dz)/(\dr)) \
& = d/(\dr) ((\du)/(\dz)) (\dz)/(\dr) + (\du)/(\dz) d/(\dr) ((\dz)/(\dr)) \ 
& = d/(\dz) (\dz)/(\dr) ((\du)/(\dz)) (\dz)/(\dr) + (\du)/(\dz) (d^2z)/(\dr^2) \
& = d/(\dz) ((\du)/(\dz)) ((\dz)/(\dr))^2 + (\du)/(\dz) (d^2z)/(\dr^2) \

& => (d^2 u)/(\dr^2) = (d^2u)/(\dz^2) ((\dz)/(\dr))^2 + (\du)/(\dz) (d^2z)/(\dr^2) \

& (d^2u_(\lm))/(\dz^2) ((\dz)/(\dr))^2 + (\du_(\lm))/(\dz) (d^2z)/(\dr^2) - (l(l+1))/(r^2) u_(\lm) = -4pi r rho_(\lm) \ 

& => {((\dz)/(\dr))^2 (d^2)/(\dz^2) + (d^2z)/(\dr^2) d/(\dz) - (l(l+1))/r^2} space u_(\lm) = -4pi r rho_(\lm) \ 
$

*上式等价于求解线性方程组: $A arrow(x) = arrow(b)$. *

$((\dz)/(\dr))^2$ 以及 $(d^2z)/(\dr^2)$ 是已知量: 
$
& ((\dz)/(\dr))^2 = (N+1)^2 (r_m)/(pi^2 r (r_m+r)^2) \ 
& ((d^2z)/(\dr^2)) = (N+1) (r_m^2 (3r+r_m))/(2 pi (r dot.c r_m)^(3\/2) (r+r_m)^2) \ 
$

${z}$ 是等距网格, 可以很自然地使用有限差分求解 $(d^2u_(\lm))/(\dz^2)$ 和 $(\du_(\lm))/(\dz)$. 

根据 ref@becke19882 的建议，使用7阶有限差分：

 - *一阶导数有限差分*
#set math.mat(delim: "{")
#align(left)[
  #grid(
    columns: 3,
    column-gutter: 3em,
    row-gutter: 1.5em,
    [采样点] , [微分公式], [矩阵表示] , 

    [-3,-2,-1,0,1,2,3] , 
    [$(d f)/(d x) approx (-f_(i-3) + 9f_(i-2) - 45f_(i-1) + 45f_(i+1) - 9f_(i+2) + f_(i+3)) / (60 \dx)$] , 
    [$mat(-1, 9, -45, 0, 45, -9, 1)$]
  )
]

- *二阶导数有限差分*
#align(left)[
  #grid(
    columns: 3,
    column-gutter: 3em,
    row-gutter: 1.5em,
    [采样点] , [微分公式], [矩阵表示] , 
    [-3,-2,-1,0,1,2,3] , 
    [$(d^2 f)/(d x^2) \ approx (2f_(i-3) - 27f_(i-2) + 270f_(i-1) - 490f_i + 270f_(i+1) - 27f_(i+2) + 2f_(i+3)) / (180 \dx^2)$] , 
    [$mat(2, -27, 270, , -490, 270, -27, 2)$]
  )
]

对于离散函数空间 ${f_i|i=0,1,2,3,4,dots.c,}$, 七阶中心差分公式只能计算到 $f^prime_3$. 在边界处使用非对称的有限差分公式计算边界点的微分. 注意通过有限差分构建微分算符矩阵不是数值求解带*初始条件*的微分方程的好方法, 因为构建的微分算符很可能是奇异的, 导致微分方程对应的线性方程组解不出来. 但是通过微分算符矩阵求解带*边界条件*的微分方程问题很合适, 只需要在微分算符矩阵首尾两行需要提供两个边界条件即可. 

- *边界非对称有限差分*
对于函数 $arrow(f) = {f_1, f_2, f_3, f_4, dots.c, f_(N-3), f_(N-2), f_(N-1), f_(N)}^T$, $f_4 space dots.c space f_(N-3)$ 处的一二阶导数可由七阶中心差分计算, 对于 $(f_2, f_3)$ 与 $(f_(N-2), f_(N-1))$ 处的微分可由下式计算:
$
& D^1 = cases(
  (\df)/(\dx)|_(x_2) = (-3f_1 - 10f_2 + 18f_3 - 6f_4 + f_5)/(12 \dx) ,
  (\df)/(\dx)|_(x_3) = (3f_1 - 30f_2 - 20f_3 + 60f_4 - 15f_5 + 2f_6)/(60 \dx) ,
  ,
  (\df)/(\dx)|_(x_(N-2)) = (-2f_(N-5) +15f_(N-4) - 60f_(N-3) + 20f_(N-2) + 30f_(N-1) - 3f_N)/(60 \dx) ,
  (\df)/(\dx)|_(x_(N-1)) = (-f_(N-4) + 6f_(N-3) - 18f_(N-2) + 10f_(N-1) + 3f_N)/(12 \dx) 
) \ 

& D^2 = cases(
  (\d^2f)/(\dx^2)|_(x_2) = (11f_1 - 20f_2 + 6f_3 + 4f_4 - f_5)/(12 \dx^2) ,
  (\d^2f)/(\dx^2)|_(x_3) = (-f_1 + 16f_2 - 30f_3 + 16f_4 - f_5)/(12 \dx^2) ,
  , 
  (\d^2f)/(\dx^2)|_(x_(N-2)) = (-f_(N-4) + 16f_(N-3) - 30f_(N-2) + 16f_(N-1) - f_N)/(12\dx^2) ,
  (\d^2f)/(\dx^2)|_(x_(N-1)) = (-f_(N-4) + 4f_(N-3) + 6f_(N-2) - 20f_(N-1) + 11f_N)/(12\dx^2) ,
) \ 
$

对于端点(边界) $(f_1, f_N)$ 处使用边界条件: $cases(f(x_1) = y_1, f(x_N) = y_N)$, 而非计算其微分. 

*例子*
#align(center)[
  #table(
    columns: 3,
    [Exponential decay], [Lorentian], [Exponential], 
    [#image("pics/fdd-example1.svg", )], 
    [#image("pics/fdd-example2.svg", )], 
    [#image("pics/fdd-example3.svg", )]
  )
]

// #include "fdd-code.inc"



= *分子DFT计算*

Kohn-Sham方程：
$
& F S = S C epsilon \ 
& F = underbrace(T + V_"ext", "analytic") + underbrace(U, "analytic\nor\nnumeric") + underbrace(K+C, "numeric")  
$

构建网格的目的在于数值计算交换相关泛函. 库伦势的计算同样可以基于网格数值计算. 

一些符号约定：
#table(
  columns: 2, 
  [`natom`], [number of atoms], 
  [`nlm`], [number of spherical pair], 
  [`nrad`], [number of radial shells/grids], 
  [`nang`], [number of angular grids], 
  [`natgrid`], [number of atomic grids, `natgrid = nrad x nang`], 
  [`ngrid`], [number of total grids, `ngrid = natom x nrad x nang = natom x natgrid`], 
  [`wi, wj`], [weight of radial point/grid i, angular point/grid j], 
  [`wa`], [Becke weight], 
  [`i, j, k`], [indices for radial grid, angular grid, spherical pair `((0,0), (1,-1), (1,0), ..., (l,l))`]

)

== 原子网格

利用之前关于高斯求积，球谐展开，球面求积的技术，*以原子为中心构建 Gauss-Chebyshev-Lebedev 网格*. 

角度网格通过列别杰夫球面采样点构建. 通过第二类切比雪夫求积构建径向网格，在每一层径向网格处，以该层半径对角度网格进行缩放，然后以原子位置为原点进行平移. 这样就得到了一层原子网格，最终共有`nrad`层原子网格，越接近原子核网格越密集，总网格数为`natgrid=nrad x nang`. 每个网格有自己的笛卡尔坐标，原子网格可以通过矩阵形式表示：

#set math.mat(delim: "(")
$
"nrad"{  overbrace(mat(
  mat(delim: "{", x_(0,0), y_(0,0), z_(0,0)), mat(delim: "{", x_(0,1), y_(0,1), z_(0,1)), dots.c, mat(delim: "{", x_(0,"nang-1"), y_(0,"nang-1"), z_(0,"nang-1")) ; 
  mat(delim: "{", x_(1,0), y_(1,0), z_(1,0)), mat(delim: "{", x_(1,1), y_(1,1), z_(1,1)), dots.c, mat(delim: "{", x_(1,"nang-1"), y_(1,"nang-1"), z_(1,"nang-1")) ; 
  dots.v, dots.v, dots.down, dots.v ; 
  mat(delim: "{", x_("nrad-1",0), y_("nrad-1",0), z_("nrad-1",0)), 
  mat(delim: "{", x_("nrad-1",1), y_("nrad-1",1), z_("nrad-1",1)), dots.c, 
  mat(delim: "{", x_("nrad-1","nang-1"), y_("nrad-1","nang-1"), z_("nrad-1","nang-1")) ; 
), "nang")
$

三维函数在网格上的函数值则是 `(nrad, nang)` 维度的矩阵：$f({bvec(r)_(i j)}) = f({x_(i j),y_(i j),z_(i j)}) = A_(i j)$. 

将径向雅可比行列式并入径向权重，通过与角度网格的权重做向量外积将总权重也组织成矩阵形式：
$
w_"atgrid" = vec(r_0^2 w^"rad"_0 , r_1^2 w^"rad"_1 , dots.v, r_"nrad-1"^2 w^"rad"_"nrad-1") dot.c mat(w^"ang"_0, w^"ang"_1, dots.c, w^"ang"_"nang-1") \ 
$

遍历径向网格，在每一层径向网格通过列别杰夫求积计算角度积分，然后乘以径向权重 $w_i$ 和径向雅可比行列式 $r_i^2$ 并求和，就可以数值计算三维函数对整个实空间的积分: 
$
& int_0^infinity r^2 d\r integral.double_Omega f(r, theta, phi) sin theta d Omega approx sum_i^"nrad" w_i r_i^2 sum_j^"nang" A_(i j) w_j  \ 
& int arrow.double.l w_i r_i^2 A_(i j) w_j \ 
$

其中最后一行是爱因斯坦求和约定. 

更有效的方法是将函数值矩阵 $A_(i j)$ 和权重矩阵 $w_"atgrid"$ 逐元素相乘然后求和. 

*例子*

- 计算 sto-3g基组下碳原子 `2px` 轨道的自重叠积分：$phi.alt(x,y,z\; a\x, a\y, a\z) = sum_i^3 c_i n_i (x-a\x) e^(-alpha_i [(x-a\x)^2 + (y-a\y)^2 + (z-a\z)^2])$. 

事实上原子的位置不影响计算结果，因为对于该原子的原子网格而言，原子的位置就是坐标原点，所以无论原子实际位于哪里，等价于原子位于坐标原点. 

计算代码如下：
```py
center = np.array([5, 1., 2.])
expns = np.array([2.9412494, 0.6834831, 0.2222899])
coefs = np.array([0.15591627, 0.60768372, 0.39195739])
nrfcs = (2. * expns / np.pi)**0.75 * (4 * expns)**0.5
f = lambda x, y, z: (x - center[0]) * coefs[0] * nrfcs[0] * np.exp(-expns[0] * ((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)) + \
                    (x - center[0]) * coefs[1] * nrfcs[1] * np.exp(-expns[1] * ((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)) + \
                    (x - center[0]) * coefs[2] * nrfcs[2] * np.exp(-expns[2] * ((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2))

xrad, rrad, wrad = gaussCheby2(75, 0.35 * 1.88972)
xleb, wleb = lebedev_rule(29) # 302 pts

xcoords = np.outer(rrad, xleb[0]) + center[0] # scaling + translation
ycoords = np.outer(rrad, xleb[1]) + center[1]
zcoords = np.outer(rrad, xleb[2]) + center[2]
A = f(xcoords, ycoords, zcoords)**2 # (nrad, nang) fnvals mat

atweight = np.outer(wrad * rrad**2, wleb) # (nrad, nang) atweight mat

quad = np.sum(A * atweight) # <2px, 2px>
print(abs(quad - 1.)) # 2.7e-8
```

利用相同的方法可以*在原子网格上计算算符矩阵元. *

如果某算符（微分算符除外）在原子网格上的函数值（通过某些手段）已知, 其可表示成矩阵形式：$hat(O) = O(bvec(r)) approx O({bvec(r)_(i j)}) = O_(i j)$. 

基函数在原子网格上的函数值也可以表示成矩阵形式：${phi.alt_mu (bvec(r)) | mu=1,2,dots.c} approx {phi.alt^mu_(i j) | mu=1,2,dots.c}$. 

矩阵元不过是一个积分：
$
& O_(mu nu) = diraccc(phi.alt_mu, hat(O), phi.alt_nu) = int underbrace(phi.alt_mu (bvec(r)) O(bvec(r)) phi.alt_nu (bvec(r)), f(bvec(r))) d^3 bvec(r) \ 
& approx sum_i^"nrad" r_i^2 w_i sum_j^"nang" underbrace(phi.alt^mu_(i j) phi.alt^nu_(i j) O_(i j), A_(i j)) w_j \ 
& O_(mu nu) arrow.double.l w_i r_i^2 {phi.alt^mu_(i j) phi.alt^nu_(i j) O_(i j)} w_j \ 
$


基函数 $phi.alt_mu$ 的中心在原子C, $phi.alt_nu$ 的中心在原子B, 要计算原子A的网格上的矩阵元只需要计算 $phi.alt_mu$, $phi.alt_nu$, $hat(O)$ 在原子A网格上的函数值. 基函数的中心不需要与网格中心一致. 




== Becke分子网格@becke1988

对于多原子分子, 空间中任意点同时受到所有原子的影响，区别在于有的原子影响很大，有的原子的影响可忽略. 即，原子在空间各点处有权重. 

引入原子权重函数（不同于 $w_i$, $w_j$ !）: $w_a (bvec(r))$. 且 $sum_a^"natom" w_a (bvec(r)) equiv 1$. 即在空间中任意一点各原子的权重相加为一. 

在 `natom` 系统中对某函数进行积分：
$
& int f(bvec(r)) d^3 bvec(r) = int {sum_a^"natom" w_a (bvec(r)) equiv 1} f(bvec(r)) d^3 bvec(r) \ 
& = sum_a^"natom" int underbrace(w_a (bvec(r)) f(bvec(r)), f_a (bvec(r))) d^3 bvec(r) \ 
$

这样就把对整个分子的积分转化为对组成分子的各个原子的积分: 
$
& int f(bvec(r)) d^3 bvec(r) = sum_a^"natom" underbrace({sum_i^"nrad" r_i^2 w_i sum_j^"nang" (w_(i j)^a f_(i j)) w_j}, "atgrid") \ 
$

权重函数 ${w_a (bvec(r)) | a=1,2,dots.c,"natom"}$ 的引入相当于在空间中按原子划分边界，要求：
1. 不能是刚性边界，边界处原子权重应光滑过渡 (not step function)
2. 边界不能过于光滑，否则无法区分原子的贡献 (switching function)

是一种边界“模糊”的沃罗诺伊原胞 (fuzzy Voronoi cells). 

定义椭圆坐标：$mu_(i j) = (r_i - r_j)/(R_(i j)) in [-1,1]$. 其中 $r_i$ 是空间点到原子i的距离，$r_j$ 是空间点到原子j的距离，$R_(i j)$ 是核间距. 

#figure(
  image("pics/fig2.svg", width: 50%), 
  caption: "Confocal elliptical coordinate system.", 
  supplement: "Figure"
)

定义 cutoff profile: $s(mu_(i j))$. 要求：
1. $s(-1) = 1$
2. $s(1) = 0$
3. $(d\s)/(d mu)|_(-1) = (d\s)/(d mu)|_1 = 0$

即要求其类似于阶梯函数，但是没有阶跃，且在原子核处连续(no cusp). 满足条件的 $s(mu)$ 的定义不唯一. 

在定义一个函数 $f(mu)$, 满足:
1. $f(-1) = -1$
2. $f(1) = 1$
3. $(d\f)/(d mu)|_(-1) = (d\f)/(d mu)|_1 = 0$

$s(mu)$ 通过 $f(mu)$ 得到：$s(mu) = 1/2 [1-f(mu)]$. 

一种满足上述要求的函数：$p(mu) = f(mu) = 3/2 mu - 1/2 mu^3$. 但是该函数过于光滑，不能很好的区分边界. 对 $p(mu)$ 进行迭代可以系统地提高区分度：
$
& f_1(mu) = p(mu) \ 
& f_2(mu) = p[p(mu)] \ 
& f_3(mu) = p{p[p(mu)]} \ 
& dots.v \ 
& s_k (mu) = 1/2 [1-f_k (mu)] \ 
$

当迭代次数非常高时 $s_k$ 的极限即为阶梯函数. 

#table(
  columns: 2,
  align: center,
  [#image("pics/fig3.svg")], 
  [#image("pics/fig4.svg")]
)

最后权重函数定义为：
$
& w_n (bvec(r)) = P_n (bvec(r)) \/ sum_m^"natom" P_m (bvec(r)) \ 
& P_i (bvec(r)) = product_(j != i) s(mu_(i j))
$



== 库伦势的构建
=== 通过电子排斥积分解析构建
=== 通过数值求解泊松方程构建
== 交换相关泛函数值求积
== 笛卡尔基函数矩阵变换到纯基函数
#set math.equation(numbering: none)
电子密度在空间任意点处的函数值已知：
$
& rho(arrow(r)) equiv rho(x,y,z) \ 
$

构建以原子为中心的分子网格：`grid[natom,nrad,(nang,3)]`. 或者分子中某一原子的网格:`grid[nrad,nang,3]`. 即径向`nrad`个网格，每一层`nang`网格，每个网格点坐标`(x,y,z)`



#set math.mat(delim: "[")
*密度矩阵*

Fock矩阵的特征向量组成LCAO分子轨道的系数: 
$
& F C = S C epsilon \
& C_(N times N) = mat(arrow(C_1) , arrow(C_2) , dots.c , arrow(C_N)) \ 
$


$
& psi_i = sum_mu C_(mu, i) phi.alt_mu \ 
& |psi_i|^2 = (sum_mu C_(mu, i) phi.alt_mu) (sum_nu C_(nu, i) phi.alt_nu)^ast \ 

& rho = sum_i 2 |psi_i|^2 = sum_mu sum_nu (sum_i 2 C_(mu,i) C^ast_(nu,i)) phi.alt_mu phi.alt_nu => sum_(mu,nu) P_(mu,nu) phi.alt_mu phi.alt_nu \ 
$

$P_(mu,nu) = sum_i^(N\/2) 2 C_(mu,i) C_(nu,i)$ 即为密度矩阵. 

通过矩阵乘法计算
```py
Puv = 2 * C[:,:ne//2] @ C[:,:ne//2].T
```



#linebreak()
*原子网格上电子密度的计算*

网格的维度为 `(natom, nrad, nang, 3+natom)`. 最后一维的前3列是网格点的笛卡尔坐标, 后几列是各原子1,2,3...i,j,...在原子i网格上的权重. 

基函数在网格上的函数值组织成维度为 `(natom, nao, nrad, nang)` 的矩阵形式. 电子密度在网格上组织成维度为 `(natom, nrad, nang)` 的矩阵. 电子密度的计算为: 
$
& rho(bold(r)) = sum_i^(N\/2) 2 |psi_i|^2 \ 
& = sum_i^(N\/2) 2 sum_mu C_(mu,i) phi.alt_mu sum_nu C_(nu,i) phi.alt_nu \ 
& = sum_(mu,nu) (sum_i^(N\/2) 2 C_(mu,i) C_(nu,i)) phi.alt_mu phi.alt_nu \ 
& = sum_(mu,nu) P_(mu,nu) phi.alt_mu phi.alt_nu \ 
$

矩阵形式: 
$
& sum_(mu,nu) P_(mu,nu) phi.alt_(A,mu,i,j) phi.alt_(A,nu,i,j)
$

其中 `A` 是原子指标, $mu, nu$ 是基函数指标, $i,j$ 分别是径向网格和角度网格指标. 

矩阵形式的好处在于可以通过矩阵乘法和求和规约避免显示嵌套循环: 
```py
rho = np.einsum('uv,auij,avij->aij', Puv, aos_vals, aos_vals)    
```

然后通过电子密度在整个实空间上的积分可以重构出电子数:
$
int rho(bold(r)) d^3 bold(r) = "ne"  \ 
$

实空间上网格是以各个原子为中心构建的, 所以上述积分是对各原子积分后的总和:
$
& int rho(bold(r)) d^3 bold(r) = sum_A w^A int rho_A (bold(r)) d^3 bold(r) \ 
& approx sum_A w^A sum_i w^A_i (r^A_i)^2 sum_j w^A_j rho^A_(i j)
$

其中 $w^A$ 是原子A在自己网格点上的权重, 原子权重 (Becke权重) 也组织为矩阵形式, 维度为 `(natom, nrad, nang, natom)`. 所以有: 
$
& int rho(bold(r)) d^3 bold(r) approx sum_A w^A sum_i w^A_i (r^A_i)^2 sum_j w^A_j rho^A_(i j) \ 
& = sum_(A i j) w^"atom"_(A,i,j,A) (w^"rad"_(A,i) (r^"rad"_(A,i))^2) rho_(A,i,j) w^"ang"_j \ 
$

可以通过求和规约避免显示循环:
```py
ne = np.einsum('aija,ai,aij,j->', grid[:,:,:,3:], becke.chebs[:,:,0]**2 * becke.chebs[:,:,1], rho, becke.wleb)
```











#set math.equation(numbering: "(1)")
= 其他

- *密度矩阵*
#set math.mat(delim: "(")
$
& "ne" = sum_i^"occ" 2 int |psi_i (bold(r))|^2 d^3 bold(r) \ 
& = sum_i^"occ" int 2 (sum_mu |C_(mu,i) phi.alt_mu|) (sum_nu C_(nu,i)  phi.alt_nu) d^3 bold(r) \ 
& = sum_mu sum_nu 2 sum_i^"occ" C_(mu,i) C_(nu,i) int phi.alt_mu^dagger phi.alt_nu d^3 bold(r) \ 
& = sum_(mu nu) P_(mu,nu) chevron.l phi.alt_mu | phi.alt_nu chevron.r \ 
& <=> sum_(mu nu) P_(mu,nu) S_(mu,nu) \ 
& equiv "Tr"{P S^T} \ 
$

直接計算"密度矩陣" $P_(mu,nu)$ 的跡不會給出電子數, 因爲*基函數不正交*. 

備註: 
- 是否 $P S^T$ 才是真正的"密度矩陣"?
- 使用該密度矩陣可否進行動力學演化? - NO


- *拉普拉斯算符*
球坐标下的拉普拉斯算符:
$
& nabla^2 = 1/r^2 partial/(partial r) (r^2 partial/(partial r)) - L^2/r^2 \ 
& L^2 = -1/(sin theta) partial/(partial theta) (sin theta partial/(partial theta)) - 1/(sin theta) partial^2/(partial phi^2)
$

- *广义对角化*
1. Symmetric orthogonalization
$
& F C = S C epsilon \ 
& C = S^(-1/2) C^prime \ 
& F S^(-1/2) C^prime = S S^(-1/2) C^prime epsilon \ 
& F S^(-1/2) C^prime = S^(1/2) C^prime epsilon \ 
& (S^(-1/2) F S^(-1/2)) C^prime = C^prime epsilon \ 
& => F^prime C^prime = C^prime epsilon \ 
$

$
& P^(-1) S P = s \ 
& S^(-1/2) = P s^(-1/2) P^(-1) \ 
because & S^(-1/2) S^(1/2) = P s^(-1/2) P^(-1) P s^(1/2) P^(-1) \ 
& = P s^(-1/2) s^(1/2) P^(-1) = I \
$

2. Canonical orthogonalization
TO BE CONTINUED...


- *Slater 交换泛函* 
$
& E_x [n(bold(r))] = int -3/4 (3/pi)^(1/3) n^(4/3) (bold(r)) d^3 bold(r) \ 
& V_x [n(bold(r))] = (delta E_x [n])/(delta n) = -(3/pi)^(1/3) n^(1/3) \ 
& E_x [n] = int n epsilon_x d^3 bold(r) \ 
& => epsilon_x [n(bold(r))] = -3/4 (3/pi)^(1/3) n^(1/3) (bold(r)) = -3/4 V_x \ 
$

交换泛函的能量: 
$
& E_"exchange" = int V_x n d^3 bold(r) = int V_x sum_i^("occ") sum_(mu,nu) 2 C_(mu,i) C_(nu,i) phi.alt_(mu) phi.alt_(nu) \ 
& = sum_(mu,nu) P_(mu,nu) int phi.alt_mu V_x phi.alt_nu d^3 bold(r) \ 
& = sum_(mu,nu) P_(mu,nu) V^x_(mu,nu) \ 
& = "Tr{P @ V.T}" \ 
$

最后乘上 $-3/4$ 得到交换能量. 

在 `libxc` 中编号为1. 


- *谱方法, 范德蒙德矩阵*



#set text(lang: "en")
#bibliography("./refs.bib")
