
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

#set text(lang: "en")
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


#heading(numbering: none, level: 1)[*Notations*]
#table(
    columns: 2,
    stroke: none,
    align: left,
    [eDFT], [electronic DFT], 
    [cDFT], [classical DFT], 
    [$psi(bvec(r))$, $f(bvec(r))$], [wavefunction],
    [$phi.alt(bvec(r))$], [basis function], 
    [$U(bvec(r))$], [Coulomb potential], 
    [$u(r)$], [Coulomb potential spherical expansion coefficient], 
    [$rho(bvec(r))$], [electron density (eDFT)，particle number density (cDFT)]

)





#outline()


#set par(first-line-indent: 2em)
/////////////////
/// main text ///
/// /////////////
= *Mathematical background*
== Gaussian type basis function
=== Gaussian product theorem
$
& g_1(arrow(r)) = e^(-alpha |arrow(r)-bold(R_A)|^2) \ 
& g_2(arrow(r)) = e^(-beta |arrow(r)-bold(R_B)|^2) \ 
& g_1 dot.c g_2 = e^(-(alpha beta)/(alpha+beta) |bold(R_A)-bold(R_B)|^2) e^(-(alpha+beta) space |arrow(r) - (alpha bold(R_A)+beta bold(R_B))/(alpha+beta)|^2) \ 
$

Define some parameters: 
$
cases(
  K = e^(-(alpha beta)/(alpha+beta) |bold(R_A)-bold(R_B)|^2), 

  p = alpha+beta, 

  bold(R_C) = (alpha bold(R_A) + beta bold(R_B))/(alpha+beta)
)
$

When Gaussian function does not center at the origin, the form of integrand needs to modifies a little bit: 
$
    int_R^infinity e^(-alpha |bvec(r) - bvec(R)|^2) (r-R)^2 sin theta d r d theta d phi <=> int_0^infinity e^(-alpha bvec(r)^2) r^2 d r d theta d phi \ 
$

The center of Gaussian function does not matter if this integral is computed in Cartesian coordinate, because all variables go from $-infinity$ to $infinity$, which is symmetric range. But in the spherical coordinate, the variable $r$ goes from 0 to $infinity$, which is not symmetric to the function center and gives solution of Gaussian error function, which is wrong. 

The change of function center is translational, which does not change the integral (are under the curve). 

But the product of two Gaussian functions give an extra term:
$
& g_1 dot.c g_2 = e^(-(alpha beta)/(alpha+beta) |bold(R_A)-bold(R_B)|^2) e^(-(alpha+beta) space |bold(r) - bold(R_C)|^2) \ 

& int g_1 dot.c g_2 space \d^3 bold(r) = K_(A B)^(alpha beta) int_(R_C)^infinity  e^(-p_(alpha beta) space |bold(r) - bold(R_C)|^2) (r-R_C)^2 \dr int \dOmega \ 

& <=> K_(A B)^(alpha beta) int_0^infinity e^(-p r^2) r^2 \dr int \dOmega \ 
$

Namely, when using Gaussian quadrature to compute the integral in the spherical coordinate, we don't need to take care of the new function center created by Gaussian product theorem, as the new function center can always be seen as the new coordinate origin, and only the pre-factor needs to be computed explicitly $K_(A B)^(alpha beta)$, which depends on the exponents and centers of each Gaussian functions involved in the product.



=== Pure basis function and Cartesian basis function<Cartbf>
#linebreak()

Primitive function (Cartesian form): $(x-A_x)^(l_x) (y - A_y)^(l_y) (z - A_z)^(l_z) e^(-alpha |bvec(r) - bvec(R)|^2)$

Contracted function (Cartesian form): $sum_i^n c_i N_i (x-A_x)^(l_(x,i)) (y - A_y)^(l_(y,i)) (z - A_z)^(l_(z,i)) e^(-alpha_i |bvec(r) - bvec(R)|^2)$

In general all primitive functions in a contraction share one set of angular momentum: $l_(alpha, i) == l_(alpha, j) ; space alpha=x,y,x$. But different primitive functions may use different angular momentums in a spherical contraction that is created via explicit linear combination. For example $d_(x^2-y^2)$: 
$
    d_(x^2-y^2) = (sqrt(3))/(2)(d_(x\x) - d_(y\y))
$

where two sets of angular momentums exist in the spherical form of $d_(x^2-y^2)$ function. But in practical calculations linear combinations are not performed a prior. Calculations are done in the Cartesian basis functions (such as compute the function values on the grid), and converted to spherical basis function values via transformation matrix after. 

Basis function (atomic orbital), is the product of angular part and radial part:  $phi.alt(bold(r)) = Y_l^m (theta,phi.alt) R(r)$, where the angular part is spherical harmonics, and the radial part is exponential. 

Spherical harmonics are the angular part solutions to the Laplace equation in the spherical coordinate via separation of variables, which is in general complex: 
$
& nabla^2 f = 1/r^2 (partial)/(partial r) (r^2 (partial f)/(partial r)) + 1/(r^2 sin theta) (partial)/(partial theta) (sin theta (partial f)/(partial theta)) + 1/(r^2 sin^2 theta) (partial^2 f)/(partial phi.alt^2) = 0 \ 

& f(bvec(r)) = R(r) Y_l^m (theta, phi.alt) \ 
& cases(
1/R d/(d\r) (r^2 (d R)/(d\r)) = lambda ,
1/Y 1/(sin theta) partial/(partial theta) (sin theta (partial Y)/(partial theta)) + 1/Y 1/(sin^2 theta) (partial^2 Y)/(partial phi.alt^2) = -lambda
) \ 
$

The separation of variables in above equations hints us that for a function in real space $rho(bvec(r))$, we can decompose it into radial part and angular part, which is actually spherical expansion, where spherical harmonics are the basis functions and the radial part is the expansion coefficient. This expansion is valid even when the radial part and angular part are coupled.  

Angular momentum operator: 
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

Analytical formular for spherical harmonics: 
$
& Y_l^m (theta, phi.alt) = sqrt((2l+1)/(4pi) ((l-m)!)/((l+m)!)) P_l^m (cos theta) e^(i m phi.alt) \ 

& int_0^(2pi) int_0^pi (Y_(l^prime)^(m^prime))^ast Y_l^m sin theta d theta d phi.alt = delta_(l l^prime) delta_(m m^prime) \ 
$

But spherical harmonics (and the linear combined real spherical harmonics) are never used in both Cartesian and pure form basis functions, instead we use the mathematically equivalent polynomials in the Cartesian coordinate (including the spherical basis function). The corresponding relation between spherical variables and their Cartesian counterparts can be found here #link("https://en.wikipedia.org/wiki/Table_of_spherical_harmonics"). 


Formally: 
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

But in practical calculations the spherical form is not used. 

In the pure basis function we take the norm of a vector $|bvec(r)-bvec(R)|$, while in the Cartesian form we *DO NOT* take the absolute value $times |x-\a_x|$. 

Pure basis function has no redundance for arbitrary angular momentum $l$. For example for $l=2$ there are 5 $d$-orbitals, for $l=3$ there are 7 $f$-orbitals, and so on and so forth. We use spherical harmonics to build up the angular part of pure form basis function, and therefore no redundance exists. 

But for the Cartesian basis function things become different. The angular part is represented as polynomial $x^(l_\x) y^(\l_y) z^(\l_z)$. When the angular momentum $l>=2$, Cartesian form basis functions have redundance, and the higher the angular momentum the more redundances Cartesian form basis functions will have:  

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

The ordering of the angular momentum in Cartesian form basis functions is usually lexicographic, (for pure form basis functions typically two ordering conventions are used, vide infra): ${x\x, x\y, x\z, y\y, y\z, z\z}$, ${x\x\x, x\x\y, x\x\z, x\y\y, x\y\z, x\z\z, y\y\y, y\y\z, y\z\z, z\z\z}$. 
```py
for lx in range(lmax, 0-1, -1):
    for ly in range(lmax-lx, 0-1, -1):
        lz = lmax - lx - ly
        sx = lx * 'x'
        sy = ly * 'y'
        sz = lz * 'z'
        print(f"({lx}, {ly}, {lz})    ({sx}{sy}{sz})")
```



=== About Normalization
Normalization factor for Cartesian form basis function: @obarasaika
$
& phi.alt(x,y,z\; alpha, l_x, l_y, l_z, bold(R)) = N (x-R_x)^(l_x) (y-R_y)^(l_y) (z-R_z)^(l_z) "Exp"{-alpha |bold(r) - bold(R)|^2} \

& N = N(alpha, l_x, l_y, l_z) = ((2 alpha)/pi)^(3/4) (4 alpha)^((l_x+l_y+l_z)/2) 1/sqrt((2l_x-1)!! space (2l_y-1)!! space (2l_z-1)!!) \
& (-1)!! = 1 \ 
$<os>

The normalization factor of pure form basis function: 
$
& phi.alt(r, theta, phi) = N Y_l^m r^l e^(-alpha r^2) \ 
& N = sqrt((2 (2 alpha)^(l + 3/2)) / (Gamma(l + 3/2)))
$

*N.B. Basis functions are NOT orthogonal and NOT normalized. *

Basis functions localize at different centers are not orthogonal, that's why the overlap matrix $S_(i\j)$ is not diagonal. But anti-intuitively yet mathematically allowed is that basis functions are not necessarily normalized. As long as they are linearly independent, they are valid basis functions in mathematics. Different software adopt different normalization conventions. Besides the coefficients hard-coded in basis set files, primitive functions can also have:
1. Each primitive function multiplies its own normalization constant(@os). Leading to all primitive functions are normalized (but normalization of contraction function is not guaranteed). 
2. Only the axis-aligned primitive functions (i.e. xx, yy, xxx, ...) are normalized, and all other primitives are multiplied to the same normalization constant $N = ((2 alpha)/(pi)) (4 alpha)^((l_"tot")/(2)) sqrt(1/(2 l_"tot" - 1)!!)$. 
3. Force all basis functions to be normalized by computing the seld-overlapping integral. 


*The diagonal elements in overlap matrix are not guaranteed to be unity!*

The simplest case is to ignore normalization at all and do calculations directly. But it brings advantages to normalize basis functions/primitive functions in advance. The first is numerical stability, avoiding the occurrence of extremely small/large values (`1e-49, 1e+107`, etc.. Condition number of overlap matrix is better); And the second advantage is better interpretability if basis functions are normalized, which means the atomic orbitals are also normalized and are physically more reasonable. But no matter which convention is used, just to *make sure to use the same convention throughout*, especially in cases that basis function values on the grid are needed (such as computing the exchange-correlation potential numerically). 

Normalization does not affect physical observables (or expectation values), but is does affect the construction of operator matrices. Assume $diracc(phi.alt_mu, phi.alt_mu) equiv 1$, then the un-normalized basis function is normalized basis function times scaling factor: 
$
& phi.alt_mu^prime = a_mu phi.alt_mu \ 
& S_(mu nu)^prime = diracc(phi.alt_mu^prime, phi.alt_mu^prime) = a_mu a_nu S_(mu nu) \ 
& "let " C_(mu i)^prime = C_(mu i) / a_mu \ 
& => f_i = sum_mu C_(mu i) phi.alt_mu = sum_mu C_(mu,i)/a_mu  a_mu phi.alt_mu = sum_mu C_(mu,i)^prime phi.alt_mu^prime \ 
$

viz. the wavefunction is unaffected. 

But operator matrices change, as their constructions depend on the choice of coordinate system and basis vectors. 
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

Coefficient matrix changes accordingly. 

Formalism of Kohn-Sham equation is unchanged, so is its eigenvalues: 
$
& F^prime C^prime = A F A dot.c A^(-1) C = A F C equiv A S C epsilon \ 
& S^prime C^prime epsilon = A S C epsilon \ 
& => F^prime C^prime = S^prime C^prime epsilon \ 
$

Density matrix after transformation becomes ($n_i$ is the occupation number): 
$
& P_(mu nu) = sum_i^"occ" n_i C_(mu i) C_(nu i) = a_mu a_nu sum_i n_i C_(mu i)^prime C_(nu i)^prime = a_mu a_nu P_(mu nu)^prime \ 
& P_(mu nu)^prime = 1/a_(mu) P_(mu nu) 1/a_(nu) \ 
& => P^prime = A^(-1) P A^(-1)
$

Expectation values are invariant because of trace invariance: 
$
& dirac(hat(O)) = "Tr"{P^prime O^prime} = "Tr"{A^(-1) P A^(-1) A O A} <=> "Tr"{O A A^(-1) P} equiv "Tr"{P O}
$

The normalization of wavefunction is retained, which is the constraint of Lagrange multiplier: 
$
& diracc(f_i, f_i) = sum_mu sum_nu diracc(C_(mu i) phi.alt_mu, C_(nu i) phi.alt_nu) \ 
& = sum_(mu nu) C_(mu i) C_(nu i) S_(mu nu) equiv 1 \ 
& <=> C^T S C equiv I \ 
& (C^prime)^T S^prime C^prime = (A^(-1) C)^T A S A A^(-1) C = C^T S C equiv I
$



#linebreak()
== Gauss quadrature
N.B This section is NOT strict mathematical definition, theorem, axiom, etc. 

Numerical integration involves taking the weighted average of the function values of a function $f(x)$ at discrete points within the interval $[a,b]$ as an approximation to the integral (based on the fundamental idea of the mean value theorem for integrals):  
$
& a <= x_0 <= x_1 <= dots.c <= x_n <= b \ 
& int_a^b f(x) d\x approx sum_(i=0)^n A_i f(x_i) \  
$
where ${A_i}$ are coefficients (weights), and ${x_i}$ are quadrature roots. 

The simplest trapezoidal method employs equally spaced nodes with equal weights at each point, which needs a significant number of discrete points to achieve satisfactory numerical accuracy. To attain higher accuracy with fewer nodes, the selection of nodes and the determination of weights are crucial. 

If there exists $x_i in [a,b]$ and coefficients ${A_i}$ such that $int_a^b f(x) d\x approx sum_{i=0}^n A_i f_i$ has algebraic precision of degree $2n+1$, then the nodes ${x_i}$ are called a *Gaussian points*, $A_i$ is called a *Gaussian coefficient*, and the quadrature formula is called a *Gaussian-type quadrature formula*. 

Typically the Gaussian points are the zeros of various orthogonal polynomials in their respective orthogonal intervals. The coefficients (weights) have defined formulae to compute, and the weighting functions of various orthogonal polynomials need to be considered when performing numerical quadrature: 
$
& int_a^b f(x) omega(x) d\x approx sum_(i=1)^n f(x_i) w_i \ 
$




=== Orthogonal polynomials
If the weighted inner product of two arbitrary polynomials in a family of polynomial is zero, this family of polynomial is orthogonal polynomial: 
$
& chevron.l p_i dot.c p_j chevron.r = int_a^b omega(x) p_i (x) p_j (x) d\x = N_(i j) delta_(i j) \ 
$

Trigonometric functions ${1, cos\x, sin\x, cos\2x, sin\2x, dots.c}$ form orthogonal functions in the interval of $[-pi, pi]$ (not polynomials). 

- If ${phi.alt_k}_(k=0)^infinity$ are weighted orthogonal in the interval of $[a,b]$, and $H_n$ is the linear space consisting of all polynomials of degree not exceeding $n$, then ${phi.alt_i | i=0, 1, dots.c, n}$ constitute a basis for $H_n$. And for all $p(x) in H_n$, $p(x) = sum_(j=0)^n c_j phi.alt_j$. 

- Orthogonal polynomial satisfies three-term recurrence relation: 
$
& phi.alt_(-1) = 0 \ 
& phi.alt_0 = 1 \ 
& phi.alt_(n+1) = (x - alpha_n) phi.alt_n - beta_n phi.alt_(n-1), space n = 0, 1, 2, dots.c \ 
& alpha_n = ((x phi.alt_n, phi.alt_n))/((phi.alt_n, phi.alt_n)) \ 
& beta_n = ((phi.alt_n, phi.alt_n))/((phi.alt_(n-1), phi.alt_(n-1)))
$

  Three-term recurrence relation is an efficient way of computing orthogonal polynomial. 

- If ${phi.alt_n}_0^infinity$ are weighted orthogonal in the interval $[a,b]$, then $phi.alt_n$ has $n$ zeros in this interval. 




=== Some commonly used orthogonal polynomials
- Legendre quadrature
Integral interval $[-1,1]$, weighting function $omega(x) = 1$, $P_0=1$, $P_1 = x$, $P_n = 1/(2^n n!) (d^n)/(d x^n) (x^2-1)^n$. 
1. $
& int_(-1)^1 P_n P_m d\x = cases(0\, space m != n, 2/(2n+1)\, space m = n) \ 
$

2. $
& P_n (-x) = (-1)^n P_n (x) \ 
$

3. $
& (n+1) P_(n+1) = (2n+1) x P_n - n P_(n-1), space n=1,2,dots.c \ 
$

4. $P_n (x)$ has $n$ different zeros in the interval $[-1,1]$. 

Quadrature formula: $int_(-1)^1 f(x) d\x approx sum_(i=1)^n f(x_i) w_i$




- Laguerre quadrature
Integral interval $[0, infinity)$, weighting function $omega(x) = e^(-x)$, $L_n = e^x (d^n)/(d x^n) (x^n e^(-x))$. 

1. $
& int_0^infinity e^(-x) L_n L_m d\x = cases(0\, space m!=n , (n!)^2\, space m=n) \ 
$

2. $
& L_0 = 1 \ 
& L_1 = 1 - x \ 
& L_(n+1) = (1+2n-x) L_n - n^2 L_(n-1), space n=1,2,dots.c 
$

Quadrature formula: $int_(0)^infinity f(x) e^(-x) d\x approx sum_(i=1)^n f(x_i) w_i$




- Hermite quadrature 
Interval $(-infinity, infinity)$, weighting function $omega(x) = e^(-x^2)$, $H_n = (-1)^n e^(x^2) (d^n)/(d x^n) (e^(-x^2))$. 

1. $
& int_(-infinity)^infinity e^(-x^2) H_n H_m = cases(0\, space m!=n , 2^n n! sqrt(pi)\, space m = n) \ 
$

2. $
& H_0 = 1 \ 
& H_1 = 2 x \ 
& H_(n+1) = 2x H_n - 2n H_(n-1), space n=1,2,dots.c \ 
$

Quadrature formula: $int_(-infinity)^infinity f(x) e^(-x^2) d\x approx sum_(i=1)^n f(x_i) w_i$



- Chebyshev quadrature of the first kind
Interval $[-1,1]$, weighting funtion $omega(x) = 1/sqrt(1-x^2)$

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

4. Zeros of $T_n$
$
& x_k = cos((2k-1)/(2n) pi) space k=1,2,dots.c, n \ 
$

Quadrature formula: $int_(-1)^1 f(x) / sqrt(1-x^2) d\x approx sum_(i=1)^n f(x_i) w_i$





- Chebyshev quadrature of the second kind
Interval $[-1,1]$, weighting function $omega(x) = sqrt(1-x^2)$

1. $
& U_n = (sin[(n+1) arccos\x])/(sqrt(1-x^2)) \ 
& U_0=1 \ 
& U_1 = 2x \ 
& U_(n+1) = 2x U_n - U_(n-1), space n=1,2,dots.c \ 
$

2. $
& int_(-1)^1 U_n U_m sqrt(1-x^2) = cases(0\, space m!=n , pi/2\, space m=n) \ 
$

Quadrature formula: $int_(-1)^1 f(x)  sqrt(1-x^2) d\x approx sum_(i=1)^n f(x_i) w_i$


*Summary*

#align(center)[
#table(
  columns: 3,
  stroke: black,
  align: left, 
  [Type], [Interval], [Weighting function], 
  [Legendre], [$[-1,1]$],        [$omega(x)=1$], 
  [Laguerre], [$[0, infinity)$], [$omega(x) = e^(-x)$], 
  [Hermite], [$(-infinity, infinity)$], [$omega(x) = e^(-x^2)$], 
  [Chebyshev I], [$(-1,1)$], [$omega(x) = 1/sqrt(1-x^2)$], 
  [Chebyshev II], [$[-1,1]$], [$omega(x) = sqrt(1-x^2)$]
)
]

*Example*

- $int_(-1)^1 1+x+x^2+x^3+x^4+x^5 \dx$


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

When using Gauss-Chebyshev quadrature, it is important to note that if the target integral is $int_(-1)^1 f(x) sqrt(1-x^2) \dx$, then for a $2n-1$ degree function $f(x)$, only $n$ nodes are required to exactly calculate the integral. However, for $int f(x) \dx$, using Chebyshev quadrature requires dividing by the weighting function first: $f(x) -> f(x) / sqrt(1-x^2)$, where the target function $(f(x))/(sqrt(1-x^2))$ is not a polynomial of finite degree! The integral result obtained with only a few nodes will have significant errors. 



=== Change of variable
Many mathematical physical functions are absolute integrable functions defined in the infinite/semi-infinite intervals $f in LL^2$. But usually Laguerre or Hermite quadrature is not used when computing this kind of functions. 

Taking Hermite quadrature as an example, since its integration interval is the entire real domain, its zeros distribute in the interval $-infinity < x_i < infinity$. If a function can be written exactly as a polynomial function multiplied by a Gaussian function, then using Hermite quadrature is a very natural choice: $F(x) = sum_j c_j x^j e^(-x^2)$. However, if it cannot be written in this form, the quadrature requires dividing by the Gaussian weighting function: $int f(x) d\x = int f(x)/e^(-x^2) e^(-x^2) d\x$. Since the zeros of Hermite polynomials are distributed throughout the entire real domain, the denominator $e^(-x^2)$ may be quite small, leading to numerical instability in $f(x)/e^(-x^2)$. Additionally, when we need to "make up" the required quadrature formula by dividing by the Gaussian weighting function, the integral kernel $f(x)/e^(-x^2)$ is undoubtedly no longer a polynomial (an infinite series). Finite-term summation can only obtain approximate results, and Hermite quadrature has lost its superiority. Increasing the number of terms can improve accuracy, but it will encounter numerical stability issues. 

But, it can still be used (and provide pretty good result). For example: $f(x) = 1/(x^2+1) e^(-x^2) + (x-1)^3 e^(-5 (x+1)^2)$

```py
import numpy.polynomial.hermite as H

# exact -43 sqrt(pi/5) / 5 + e pi Erfc[1] = -5.4736295302356037725
f = lambda x: (1/(x**2+1)) * np.exp(-x**2) + (x-1)**3 * np.exp(-5 * (x+1)**2)
res = H.hermgauss(128)
quad = np.sum(f(res[0]) / np.exp(-res[0]**2) * res[1]) # -5.47362953023580089962, AbsErr = 1e-13
# N.B. np.exp(-res[0]**2) min=1e-102
```

Another method for numerically calculating improper integrals is through the use of Legendre quadrature. One of the advantages of Legendre quadrature lies in its simplest weighting function, $omega(x)=1$. Additionally, since its integration interval is $[-1,1]$ and its zeros are also confined within this range, numerically there will be no occurrence of illy small or large values. Similarly, due to the boundedness of its effective interval, variable transformation is required to map the bounded interval to a (pseudo-)unbounded interval. 

For the improper integral of an absolute integrable function on a semi-infinite interval, numerical quadrature cannot handle infinite boundaries. However, considering that $f|_(x=infinity)=0$, the calculation of the improper integral can be approximated as the calculation of a definite integral with a very large upper bound:
```
N[Integrate[Exp[-x^2], {x, 0, Infinity}] - Integrate[Exp[-x^2], {x, 0, 10^4}]]; 
Out[]= 0.
```

To map $(-1,1)$ to $(0,"Big")$, define (one possible) mapping function: 
$
& x in (-1,1) mapsto T(x) in (0, "Big") \ 
& "left" = 1e-3 \ 
& "right" = 1e+4 \ 
& alpha = ln("right"/"left") \ 
& T(x) = "left" (e^(alpha (x+1)/2) - 1) \ 
& T^prime (x) = alpha/2 "left" e^(alpha (x+1)/2) \ 
$

and the calculation of integral $int_0^infinity F(x) d\x$ can be approximated as: 
$
& int_0^infinity F(x) d\x approx int_0^"Big" F[T(x)] d T(x) <=> int_(-1)^1 F[T(x)] T^prime (x) d\x \ 
& approx sum_(i=1)^n F[T(x_i)] underbrace((T^prime (x_i)), "Jacobian") w_i
$

where $T^prime (x)$ is the resulting Jacobian determinant due to the change of variable (1D scalar). 

*Example* 

- Compute $int_0^infinity 1/(x^2+1) e^(-x^2) + (x-1)^3 e^(-5 (x+1)^2)$ via Legendre quadrature + change of variable

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


Similarly, the same trick can be applied to Chebyshev quadrature of the second kind. Different from the Legendre quadrature, performing Chebyshev quadrature needs to think about the weighting function. But its quadrature roots and weights are very easy to compute, and *the roots can be obtained from arithmetic sequence*:  
$
& z = 1,2, dots.c, N space "equally spaced grid" \ 
& x_i = cos(z_i/(N+1) pi) in (1,-1) space "root" \ 
& w_i = pi/(N+1) sin^2(z_i/(N+1) pi) space "weight" \ 
$

Mapping function suggested by Becke: $x mapsto r=r(x) = r_m (1+x)/(1-x) in ("Big", 0) $, where $r_m$ is the element Bragg-Slater radius@becke1988@slaterradii. This design on one hand makes Chebyshev quadrature possible to deal with integrals in semi-infinite interval, and on other hand it is dense near the nucleus and sparse away from it, just like how electron density behaves. 

#figure(
  image("pics/fig1.svg", width: 50%), 
  caption: "Becke mapping function.", 
  supplement: "Figure"
)

$z in [1,N]$, $x in (1,-1)$, $r in (infinity, 0)$. The final variable ${r_i|i=1,2,dots.c,N}$ is a function of the original variable ${z_i}$, so all functions with $r$ as a variable are functions of the variable $z$ through the chain rule. This brings a benefit: if the derivative of $f(r)$ with respect to $r$ needs to be calculated through finite difference, the problem can be transformed into calculating the derivative of $f(r)$ with respect to $z$ through the chain rule, while $Delta z=1$. If Legendre quadrature is used, the original variable $z$ is not on a uniform grid, and finite difference is not easy to operate on non-uniform grids

Chebyshev quadrature requires division by the weighting function. Taking into account the Jacobian determinant brought about by change of variable, these factors can ultimately be incorporated into the weights: 
$
& (d\r)/(d\x)|_(x_i) = 2 r_m 1/(1-x_i)^2 \ 
& int_0^infinity f(r) r^2 d\r approx sum_(i=1)^N (f[r(x_i)] r_i^2) / sqrt(1-x_i^2) pi/(N+1) sin^2(z_i/(N+1) pi) 2 r_m 1/(1-x_i)^2 \ 
& = sum_(i=1)^N f(r_i) r_i^2 pi/(N+1) sin^2(z_i/(N+1) pi) (2 r_m) / (sqrt(1-x_i^2) (1-x_i)^2) \ 
& => w_i = pi/(N+1) sin^2(z_i/(N+1) pi) (2 r_m) / (sqrt(1-x_i^2) (1-x_i)^2) \ 
$

Note that the roots are arranged in a reverse order: ${x_i} in (1,-1)$, mapped variables ${r_i} in ("Big", 0)$ are also in the reverse order. The weights ${w_i}$ need to correspond to this. 

*Example*
- $int_0^infinity 1/(r^2+1) e^(-r^2) + (r-1)^3 e^(-5 (r+1)^2) r^2 d\r$. Compute the radial part integral for a 3D function (1D problem)
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



== Numeric quadrature on the sphere
=== Real spherical harmonics<rsh> 

Spherical harmonics are eigenfunctions of angular momentum operator ($hat(L)^2$) and its $z$-component operator ($hat(L)_z$): 
$
& hat(L)^2 Y_l^m (theta, phi) = l(l+1) hbar^2 Y_l^m (theta, phi) \
& hat(L)_z Y_l^m (theta, phi) = m hbar Y_l^m (theta, phi) \ 
$

Spherical harmonics are usually complex: 
$
& Y_l^m (theta, phi.alt) = (-1)^m N_(l m) P_l^m (cos theta) e^(i m phi.alt) \ 
& N_(l m) = sqrt((2l+1)/(4pi) ((l-m)!)/((l+m)!)) \ 
$

$(-1)^m$ is Condon-Shortley phase and $P_l^m$ is associated-Legendre polynomial. 

Spherical harmonics are orthonormal on unit sphere: 
$
& int_0^pi int_0^(2pi) Y_l^m space Y_(l^prime)^(m^prime) sin theta d theta d phi = delta_(l l^prime) delta_(m m^prime) = chevron.l Y_l^(m)mid(|)Y_(l^prime)^(m^prime) chevron.r_Omega \ 
$

For the convenience of practical calculations, and since the linear superposition of solutions to a differential equation remains a solution to the original equation (superposition principle), it is common in practice to perform linear combinations on spherical harmonics to eliminate the imaginary part, resulting in what are known as real spherical harmonics. Both complex and real spherical harmonics are typically functions of the angles $(theta,phi)$: $Y_l^m (Omega) = Y_l^m (theta,phi)$. The angular variables can be transformed into Cartesian coordinates through inverse trigonometric functions, so $Y_l^m (theta,phi) equiv Y_l^m (x,y,z)$. 

The construction of real spherical harmonics is not unique, but orthonormality must be fulfilled: $& int_0^(2pi) int_0^pi (Y_(l^prime)^(m^prime))^ast Y_l^m sin theta d theta d phi.alt = delta_(l l^prime) delta_(m m^prime)$. 

Below are a few possible constructions: 
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

And possible implementations: 
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


4. C++ std library implemented Legendre polynomial, which is ready to use: 

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



=== Lebedev quadrature
Lebedev quadrature is the most commonly used angular quadrature scheme@lebedev1975values. 

- *Example*@beentjes2015quadrature

$int_Omega f(x,y,z) d Omega = int_Omega 1+x+y^2+x^2 y + x^4 + y^5 + x^2 y^2 z^2$

```py
import numpy as np
from scipy.integrate import lebedev_rule

f = lambda x,y,z: 1 + x + y**2 + x**2 * y + x**4 + y**5 + x**2 * y**2 * z**2

res = lebedev_rule(17) # 110 angular pts
quad = np.sum(f(*res[0]) * res[1]) # 19.388114662154152
# exact 216pi/35 = 19.388114662154152
```

The quadrature result can be verified via Mathematica ($r=1$ on unit sphere): 
```txt
f[x_,y_,z_]:=1+x+y^2+x^2 y+x^4+y^5+x^2 y^2 z^2;
Integrate[(f[x,y,z]/.{x->Sin[\[Theta]] Cos[\[Phi]],y->Sin[\[Theta]] Sin[\[Phi]],z->Cos[\[Theta]]}) Sin[\[Theta]],{\[Theta],0,Pi},{\[Phi],0,2 Pi}]
Out[]=216Pi/35
```





== Numeric solution to Poisson's equation<POISSON>
=== First order functional derivative

Analogous to the derivative of a function: $f^prime (x) = (d f(x))/(d\x)$, $d\f = f^prime d\x$. The first-order functional derivative is: 
$
  & F[f(x)] \
  & F'[f(x)] = (delta F[f(x)]) / (delta f(x)) \
  & delta F = F[f(x) + delta f(x)] - F[f(x)] = underbrace(integral ((delta F) / (delta f)) delta f \dx, "contributions from all x") + cal(O)(delta^2 f) \
$

Unlike the derivative of a function, there is no formula for the functional derivative. The derivative result depends on the specific form of the functional. The coefficient of $delta f$ in the integral is the functional derivative. 

*Example*

- Compute the first-order functional derivative of $F[n(x)] = integral_0^1 n^2 \dx$: 
$
  & delta F = integral_0^1 (n + delta n)^2 - n^2 \dx \
  & delta F = integral_0^1 2n delta n \dx \
  & => (delta F) / (delta n) = 2n \
$

*Example*

- First order derivative of integral-type functional: 

Assume $delta n$ is fixed at the boundaries: $delta n = 0$，therefore $(partial f) / (partial n') delta n mid(|)_"boundary" = 0$ 
$
  & F[n] = integral f(n, n', x) \dx \
  & F[n+delta n] = integral f(n+delta n, d / (\dx) (n + delta n), x) \dx \
  & \
  & "Taylor expansion to the first term" \
  
  & => F[n+delta n] = int f(n, n^prime, x) + {(partial f)/(partial n) delta n + (partial f)/(partial n^prime) (delta n)^prime} space \dx \
  
  & F[n+delta n] - F[n] = delta F = int (partial f)/(partial n) delta n + (partial f)/(partial n^prime) (delta n)^prime space \dx \
  & \

  & "do integration by parts to the second term对第二项进行分部积分" \
  & => delta F = int (partial f)/(partial n) delta n space \dx + int (partial f)/(partial n^prime) d (delta n) = int (partial f)/(partial n) delta n space \dx + (partial f)/(partial n^prime) delta n |_("边界") - int delta n space d ((partial f)/(partial n^prime)) \

  & => delta F = int (partial f)/(partial n) delta n space \dx - int delta n space d/(\dx) ((partial f)/(partial n^prime)) space \dx \

  & delta F = int (delta F[n])/(delta n) space delta n space \dx = int {(partial f)/(partial n) - d/(\dx) ((partial f)/(partial n^prime))} space delta n space \dx \

  & therefore (delta F[n])/(delta n) = (partial f)/(partial n) - d/(\dx) ((partial f)/(partial n^prime))
$

And this is Euler-Lagrange equation. 


The functional derivative of Coulomb energy $E[n(bvec(r))]_"Coul"$ w.r.t. electron density gives Coulomb potential $v[n(bvec(r))]_"Coul" = (delta E[n])/(delta n)$: 

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


=== Poisson's equation

An important identity: 
$ 
& nabla^2_bvec(r) 1/(|bvec(r)-bvec(r)^prime|) = -4 pi delta(bvec(r)-bvec(r)^prime) \
$<green>

where $1/(|bvec(r) - bvec(r)^prime|)$ is Green's function. Subscript in $nabla^2_bvec(r)$ indicates applying Laplacian to $bvec(r)$.  

@hartreepotential is the convolution of the Coulomb kernel function with the electron density: 

$
& (f * g)(t) = int f(tau) g(-tau+t) d tau \
& f = n(bvec(r)) \
& g = 1/(|bvec(r)|) \ 
& => (f * g)(bvec(r)) = int n(bvec(r^prime)) 1/(|-bvec(r^prime) + bvec(r)|) \d^3 bvec(r^prime) = int (n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r)^prime \
$

By plugging @green into @hartreepotential we get: 
$
& nabla_bvec(r)^2 v_"Ha" (bvec(r)) = nabla^2 int (n(bvec(r^prime)))/(|bvec(r-r^prime)|) d^3 bvec(r^prime) = int n(bvec(r^prime)) nabla^2 1/(|bvec(r)-bvec(r)^prime|) d^3 bvec(r^prime) \ 

& = int n(bvec(r^prime)) dot.c -4 pi delta(bvec(r)-bvec(r)^prime) d^3 bvec(r^prime) = -4 pi n(bvec(r)) \
$<hartreepotential2>

And we get the *integral form and derivative form of Coulomb potential*, which are mathematically equivalent: 
$
cases(
v_"Ha" (bvec(r)) = int (n(bvec(r)^prime))/(|bvec(r)-bvec(r)^prime|) \d^3 bvec(r)^prime ,
nabla^2 v_"Ha" (bvec(r)) = -4 pi n(bvec(r))
) \ 
$<poisson>

The derivative form is the *Poisson equation*. Solving Poisson equation is to calculate the Coulomb potential at given electron density.  






=== Spherical expansion

Since (real) spherical harmonics are orthonormal on the unit sphere ($theta in (0, pi)$, $phi.alt in (0, 2 pi)$, $r equiv 1$), they constitute a complete basis set on the unit sphere. Any absolute integrable function on the unit sphere can be expressed as a linear combination of spherical harmonics: 
$
f(theta, phi.alt) = sum_(l=0)^infinity sum_(m=-l)^l a_(l m) Y_l^m (theta, phi), space space l in bb(N) \
$

Basis function in the spherical coordinate is the product of spherical harmonics and Gaussian functions (N.B. primitive gaussian normalized, contracted gaussian normalization not guaranteed)：
$
phi.alt(r,theta,phi) = Y_l^m (theta,phi) |bvec(r)-bvec(R)|^l sum_i c_i n_i e^(-alpha_i |bvec(r)-bvec(R)|^2)
$ 

The spherical form is equivalet to the Cartesian form: 
$
phi.alt(x,y,z) = (x-R_x)^(l_x) (y-R_y)^(l_y) (z-R_z)^(l_z) sum_i c_i n_i e^(-alpha_i |bvec(r)-bvec(R)|^2) \ 
$

The radial part and angular part are coupled and cannot be explicitly separated, but the radial part can be extracted via spherical expansion: 
$
& |phi.alt(r,theta,phi)|^2 equiv |phi.alt(x,y,z)|^2 = sum_(l m) rho_(l m) (r) Y_l^m (theta,phi) \ 
& rho_(l m) (r) = int Y_l^m (theta,phi)^ast |phi.alt(r,theta,phi)|^2 d Omega \ 
& rho_(l m) (r_i) approx sum_j^M |phi.alt(r_i, Omega_j)|^2 Y_l^m (Omega_j)^ast w_j
$<sphericalexpansion>

$Omega_j$ in @sphericalexpansion  is the sampling point on the unit sphere (built via Lebedev quadrature), and $w_j$ is the quadrature weight. 


// *不需要根据原子中心位置进行平移或缩放*，原子网格在构建时已经平移缩放过了，这里出现的$Y_l^m$只是为了通过投影计算展开系数. 

The equation is valid only when the summation in @sphericalexpansion goes to infinity; unless the function to be expanded is itself spherical harmonics or power of spherical harmonics, in which case the expansion is finite and the expansion coefficient can be determined by Clebsch-Gordan coefficient. In common cases, Becke suggested $l$ to be half the Lebedev order. For example, `nleb=29 (302 angular pts)`, `lmax=14`, `nlm=225` @becke19882. 

=== Solving Poisson's equation via spherical expansion and finite difference

Do spherical expansion on electronic Coulomb potential $U(bvec(r))$ and electron density $rho(bvec(r))$: 
$
& U(bvec(r)) = U(r,theta,phi) = sum_(\lm) (u_(\lm)(r))/r Y_l^m (theta, phi) \
& rho(bvec(r)) = rho(r,theta,phi) = sum_(\lm) rho_(\lm)(r) Y_l^m (theta,phi) \
$

plug in @poisson

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

For each pair of $(l,m)$ we get one 2nd-order ODE: 
$
& d^2/(d r^2) u_(\lm) (r) - (l(l+1))/r^2 u_(\lm) (r) = -4pi r rho_(\lm) (r) \
$

The problem is simplified from solving one 2nd-order PDE to solving a set of 2nd-order ODE. 

Uniform grid ${z_i|i=1,2,dots.c,N}$ is converted to the zeros of Chebyshev of the second kind: ${x_i = cos (z_i)/(N+1) pi | i=1,2,dots.c,N} in (1,-1)$, which is then mapped to semi-infinite interval: ${r_i = r_m (1+x_i)/(1-x_i) | i=1,2,dots.c, N} in ("Big",0)$. 

2nd-order derivative function needs at least two boundary conditions, which are determinied by $(z_0, x_0, r_0)$ and $(z_(N+1), x_(N+1), r_(N+1))$ (Dirichlet condition). 

By using chain rule, the radial part of Coulomb potential is a function of scalar variable $r$: $u_(\lm) = u_(\lm) (r)$, and therefore is a function of uniform grid ${z_i}$: $u_(\lm) = u_(\lm)[r(z)]$:

$
& (d^2 u(r))/(d r^2) = d/(\dr) (d u)/(d r) = d/(\dr) ((\du)/(\dz) (\dz)/(\dr)) \
& = d/(\dr) ((\du)/(\dz)) (\dz)/(\dr) + (\du)/(\dz) d/(\dr) ((\dz)/(\dr)) \ 
& = d/(\dz) (\dz)/(\dr) ((\du)/(\dz)) (\dz)/(\dr) + (\du)/(\dz) (d^2z)/(\dr^2) \
& = d/(\dz) ((\du)/(\dz)) ((\dz)/(\dr))^2 + (\du)/(\dz) (d^2z)/(\dr^2) \

& => (d^2 u)/(\dr^2) = (d^2u)/(\dz^2) ((\dz)/(\dr))^2 + (\du)/(\dz) (d^2z)/(\dr^2) \

& (d^2u_(\lm))/(\dz^2) ((\dz)/(\dr))^2 + (\du_(\lm))/(\dz) (d^2z)/(\dr^2) - (l(l+1))/(r^2) u_(\lm) = -4pi r rho_(\lm) \ 

& => {((\dz)/(\dr))^2 (d^2)/(\dz^2) + (d^2z)/(\dr^2) d/(\dz) - (l(l+1))/r^2} space u_(\lm) = -4pi r rho_(\lm) \ 
$

*Equation above is equivalent to $A arrow(x) = arrow(b)$. *

$((\dz)/(\dr))^2$ and $(d^2z)/(\dr^2)$ are known:  
$
& ((\dz)/(\dr))^2 = (N+1)^2 (r_m)/(pi^2 r (r_m+r)^2) \ 
& ((d^2z)/(\dr^2)) = (N+1) (r_m^2 (3r+r_m))/(2 pi (r dot.c r_m)^(3\/2) (r+r_m)^2) \ 
$

${z}$ is uniform grid, and finite difference discretization (fdd) is naturally applied to calculate $(d^2u_(\lm))/(\dz^2)$ and $(\du_(\lm))/(\dz)$. 

As suggested by ref@becke19882, we use 7th-order fdd: 

 - *First order derivative fdd*
#set math.mat(delim: "{")
#align(left)[
  #grid(
    columns: 3,
    column-gutter: 3em,
    row-gutter: 1.5em,
    [Sampling point] , [Derivative formula], [Matrix representation] , 

    [-3,-2,-1,0,1,2,3] , 
    [$(d f)/(d x) approx (-f_(i-3) + 9f_(i-2) - 45f_(i-1) + 45f_(i+1) - 9f_(i+2) + f_(i+3)) / (60 \dx)$] , 
    [$mat(-1, 9, -45, 0, 45, -9, 1)$]
  )
]

- *Second order derivative fdd*
#align(left)[
  #grid(
    columns: 3,
    column-gutter: 3em,
    row-gutter: 1.5em,
    [Sampling point] , [Derivative formula], [Matrix representation] , 
    [-3,-2,-1,0,1,2,3] , 
    [$(d^2 f)/(d x^2) \ approx (2f_(i-3) - 27f_(i-2) + 270f_(i-1) - 490f_i + 270f_(i+1) - 27f_(i+2) + 2f_(i+3)) / (180 \dx^2)$] , 
    [$mat(2, -27, 270, , -490, 270, -27, 2)$]
  )
]

For the discrete function space ${f_i|i=0,1,2,3,4,dots.c}$, the seventh-order central finite difference formula can only calculate up to $f^prime_3$. At the boundary, an asymmetric finite difference formula is used to calculate the differential of the boundary points. Note that constructing a differential operator matrix through finite differences is not a good method for numerically solving differential equations with *initial conditions*, as the constructed differential operator is likely to be singular, leading to the inability to solve the corresponding linear equation system of the differential equation. However, solving differential equations with *boundary conditions* through a differential operator matrix is quite suitable, as it only requires providing two boundary conditions in the first and last rows of the differential operator matrix. 


- *Asymmetric finite difference at boundaries*
For function $arrow(f) = {f_1, f_2, f_3, f_4, dots.c, f_(N-3), f_(N-2), f_(N-1), f_(N)}^T$, the first- and second-order derivative at $f_4 space dots.c space f_(N-3)$ can be calculated by 7th-order central fdd; for derivatives at $(f_2, f_3)$ 与 $(f_(N-2), f_(N-1))$ can be calculated by: 
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

At the boundaries $(f_1, f_N)$ we use boundary conditions: $cases(f(x_1) = y_1, f(x_N) = y_N)$, instead of calculating their derivatives. 

*Example*
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



= *Molecular DFT calculation*

Kohn-Sham equation: 
$
& F S = S C epsilon \ 
& F = underbrace(T + V_"ext", "analytic") + underbrace(U, "analytic\nor\nnumeric") + underbrace(K+C, "numeric")  
$

Gird is used to numerically calculate exchange-correlation functional. Coulomb potential can also be numerically calculated on the grid. 

Notations : 
#table(
  columns: 2, 
  [`natom`], [number of atoms], 
  [`nbf`], [number of basis function], 
  [`nlm`], [number of spherical pair], 
  [`nrad`], [number of radial shells/grids], 
  [`nang`], [number of angular grids], 
  [`natgrid`], [number of atomic grids, `natgrid = nrad x nang`], 
  [`ngrid`], [number of total grids, `ngrid = natom x nrad x nang = natom x natgrid`], 
  [`wi, wj`], [weight of radial point/grid i, angular point/grid j], 
  [`wa`], [Becke weight (atomic weight)], 
  [`i, j, k/lm`], [indices for radial grid, angular grid, spherical pair `((0,0), (1,-1), (1,0), ..., (l,l))`]

)

== Atomic grid

Based on previously mentioned Gauss quadrature, spherical expansion, and quadrature on the unit sphere, we can build atomically centered *Gauss-Chebyshev-Lebedev grid*. 

The angular grid is constructed on the Lebedev sphere sampling points. The radial grid is constructed using the second kind of Chebyshev quadrature. At each layer of the radial grid, the angular grid is scaled by the radius of that layer, and then translated with the atomic position as the origin. This results in a layer of atomic grid. There are ultimately `nrad` layers of atomic grids, with the grid closer to the atomic nucleus being denser. The total number of grids is `natgrid=nrad x nang`. Each grid has its own Cartesian coordinates, and the atomic grid can be represented in matrix form: 

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

3D function values on the grid are organized in matrix with shape `(nrad, nang)` : $f({bvec(r)_(i j)}) = f({x_(i j),y_(i j),z_(i j)}) = A_(i j)$. 

Embedding the radial Jacobian into radial weight, and the vector outer product with angular weight gives another matrix representation: 
$
w_"atgrid" = vec(r_0^2 w^"rad"_0 , r_1^2 w^"rad"_1 , dots.v, r_"nrad-1"^2 w^"rad"_"nrad-1") dot.c mat(w^"ang"_0, w^"ang"_1, dots.c, w^"ang"_"nang-1") \ 
$

Traversing the radial grid, the angular integral is calculated at each layer of the radial grid through the Lebedev quadrature, and then multiplied by the radial weight $w_i$ and the radial Jacobian determinant $r_i^2$ and summed up. This allows for the numerical calculation of the integral of a three-dimensional function over the entire real space:

$
& int_0^infinity r^2 d\r integral.double_Omega f(r, theta, phi) sin theta d Omega approx sum_i^"nrad" w_i r_i^2 sum_j^"nang" A_(i j) w_j  \ 
& int arrow.double.l w_i r_i^2 A_(i j) w_j \ 
$

The last line assumes Einstein summation convention is used. 

A computationally more efficient way is doing element-wise matrix multiplication between function values $A_(i j)$ and weight $w_"atgrid"$, and finally sum up all elements. 

*Example*

- Calculate the carbon `2px` orbital self-overlap integral at `sto-3g` basis set: $phi.alt(x,y,z\; a\x, a\y, a\z) = sum_i^3 c_i n_i (x-a\x) e^(-alpha_i [(x-a\x)^2 + (y-a\y)^2 + (z-a\z)^2])$. 

The atomic position does not affect the calculation results, because for the atomic grid of that specific atom, its position serves as the coordinate origin. Therefore, regardless of where the atom is actually located, it is equivalent to the atom being at the coordinate origin. 

Below is one possible way of doing it: 
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

Operator matrix element on the atomic grid can be calculated using the same trick. 

If operator's values, except for derivative operators, on the atomic grid are known by any means, this operator can be represented in matrix form: $hat(O) = O(bvec(r)) approx O({bvec(r)_(i j)}) = O_(i j)$. 

Basis function values on the atomic grid can be represented in matrix form as well: ${phi.alt_mu (bvec(r)) | mu=1,2,dots.c} approx {phi.alt^mu_(i j) | mu=1,2,dots.c}$. 

And the matrix element is nothing but an integral: 
$
& O_(mu nu) = diraccc(phi.alt_mu, hat(O), phi.alt_nu) = int underbrace(phi.alt_mu (bvec(r)) O(bvec(r)) phi.alt_nu (bvec(r)), f(bvec(r))) d^3 bvec(r) \ 
& approx sum_i^"nrad" r_i^2 w_i sum_j^"nang" underbrace(phi.alt^mu_(i j) phi.alt^nu_(i j) O_(i j), A_(i j)) w_j \ 
& O_(mu nu) arrow.double.l w_i r_i^2 {phi.alt^mu_(i j) phi.alt^nu_(i j) O_(i j)} w_j \ 
$<atomicoperator>

Basis function $phi.alt_mu$ centers on atom C, $phi.alt_nu$ centers on atom B, to compute the matrix element on atom A's grid only the discrete values of $phi.alt_mu$, $phi.alt_nu$, $hat(O)$ on A's grid are needed. The centers of basis functions are not necessarily aligned with the center of grid. 




== Becke molecular grid@becke1988

For polyatomic molecule, all positions in the space are affected by all atoms simultaneously, but some atoms may have larger influence and some atoms' influence are negligible. Namely, atoms have weights in different positions in the space. 

Introducing atomic weighting function (different from $w_i$m and $w_j$): $w_a (bvec(r))$, and requires $sum_a^"natom" w_a (bvec(r)) equiv 1$. So that the summation of weights of all atoms in any position in space gives unity. 

Performing integral in a `natom` system: 
$
& int f(bvec(r)) d^3 bvec(r) = int {sum_a^"natom" w_a (bvec(r)) equiv 1} f(bvec(r)) d^3 bvec(r) \ 
& = sum_a^"natom" int underbrace(w_a (bvec(r)) f(bvec(r)), f_a (bvec(r))) d^3 bvec(r) \ 
$

The integral over the whole molecule is factored into integrals over each atomic component: 
$
& int f(bvec(r)) d^3 bvec(r) = sum_a^"natom" underbrace({sum_i^"nrad" r_i^2 w_i sum_j^"nang" (w_(i j)^a f_(i j)) w_j}, "atgrid") \ 
$

The introduction of atomic weighting function ${w_a (bvec(r)) | a=1,2,dots.c,"natom"}$ brings boundaries in space that are based on atoms, which requires: 
1. Cannot be hard, atomic weights near boundaries should transition smoothly (not step function). 
2. Boundaries cannot be too smooth, atomic contributions should be distinguishable (switching function). 

In short, we need a boundary blurred Voronoi cell. 

Define elliptical coordinate: $mu_(i j) = (r_i - r_j)/(R_(i j)) in [-1,1]$, where $r_i$ is the distance between a point and atom i, $r_j$ is the distance between point and atom j, and $R_(i j)$ is the inter-nuclear distance. 

#figure(
  image("pics/fig2.svg", width: 50%), 
  caption: "Confocal elliptical coordinate system.", 
  supplement: "Figure"
)

Define cutoff profile: $s(mu_(i j))$, which requires：
1. $s(-1) = 1$
2. $s(1) = 0$
3. $(d\s)/(d mu)|_(-1) = (d\s)/(d mu)|_1 = 0$

It should behave like step function but without jumping, and should be continuous near the nucleus (no cusp). Definition of function that meets all these requirements is not unique. 

Define function $f(mu)$:
1. $f(-1) = -1$
2. $f(1) = 1$
3. $(d\f)/(d mu)|_(-1) = (d\f)/(d mu)|_1 = 0$

and $s(mu)$ is obtained from $f(mu)$: $s(mu) = 1/2 [1-f(mu)]$. 

One function that satisfies all the criteria: $p(mu) = f(mu) = 3/2 mu - 1/2 mu^3$. But this function is too smooth and boundaries are not well distinguished. Repeating iteration on $p(mu)$ can systematically increase its separation capability: 
$
& f_1(mu) = p(mu) \ 
& f_2(mu) = p[p(mu)] \ 
& f_3(mu) = p{p[p(mu)]} \ 
& dots.v \ 
& s_k (mu) = 1/2 [1-f_k (mu)] \ 
$

When the iteration is very high, the limit of $s_k$ is step function. 

#table(
  columns: 2,
  align: center,
  [#image("pics/fig3.svg")], 
  [#image("pics/fig4.svg")]
)

And finally the atomic weighting function is obtained: 
$
& w_n (bvec(r)) = P_n (bvec(r)) \/ sum_m^"natom" P_m (bvec(r)) \ 
& P_i (bvec(r)) = product_(j != i) s(mu_(i j))
$

*Core idea of Becke grid: deal with atomic grid first, then loop over all atoms. *


== Building up Coulomb potential
=== Building up Coulomb potential analytically via electron repulsion integral

Same as wavefunction-based method: 
$
  J_(mu nu) = sum_(lambda sigma) P_(lambda sigma) (mu nu|lambda sigma) \ 
$

Coulomb potential is calculated by analytically evaluating tensor contraction between ERI and electron density. 

=== Building up Coulomb potential numerically via solving Poisson's equation@becke19882

Problem to be solved: 
$
& U(bvec(r)) = sum_(l m) (u_(l m) (r)) / r  Y_(l m) (Omega) \ 
& rho(bvec(r)) = sum_(l m) rho_(l m) (r) Y_(l m) (Omega) \ 
& {((\dz)/(\dr))^2 (d^2)/(\dz^2) + (d^2z)/(\dr^2) d/(\dz) - (l(l+1))/r^2} space u_(\lm) = -4pi r rho_(\lm) \ 
$

which is equivalent to solving a set of linear equations: $A arrow(x) = arrow(b)$. Both second-order derivative operator and first-order derivative are heptadiagonal matrices, $(l(l+1)) / r^2$ is diagonal, and $(d z)/(d r)$is constant. 

2nd-order ODE needs two boundary conditions, there are `nrad` radial grids so the shape of matrix $A$ is `(nrad+2, nrad+2)`. $rho_(l m)$ is spherical expansion coefficient, which has shape `(nrad+2, )`, where the first and last elements are boundary conditions for $r->infinity$ and $r->0$ (Dirichlet boundary condition). Note that ${r_i}$ is reversely ordered. 



- *Boundary conditions*

When $r->infinity$，any charge distributions can be seen as point changes no matter how complex they actually are, so $lim_(r->infinity) U(bvec(r)) = U(r) = q/r = sum_(l m) u_(l m)/r Y_(l m)$. The potential now is a function of scalar $r$ only, therefore only the coefficient of $Y_(00)$ is non-zero: $q/r = u_00 / r sqrt(1/(4pi)) => u_00 = sqrt(4 pi) q$. 

When $r->0$, $u_(l m) = 0$. 

N.B. $q$ is partial charge, not atomic number. 

```py
mweights = np.array([np.outer(rrad[iatom]**2 * wrad[iatom], wleb).ravel() * watom[iatom] for iatom in range(natom)]): # (natom, natgrid)
    qn = np.sum(rho[iatom] * mweights[iatom]) # partial charge for atom i
```



- *Spherical basis*

Define spherical basis $Y_(l m, j) = Y_(k j)$ `shape (nlm, nang)`, which is shared by all atoms and the angular sampling points are the same as Lebedev quadrature points. 

```py
def _build_spherical_basis(lm:list[(int,int)], x:np.ndarray, y:np.ndarray, z:np.ndarray) -> np.ndarray:
    r'''
    Arguments
    ---------
    lm : list[(int,int)]
        Angular number. [(0,0), (1,-1), (1,0), ..., (lmax,lmax)]
    x, y, z : np.ndarray:
        Cartesian coordinates of lebedev sampling points on unit sphere

    Return
    ------
    y_kj : np.ndarray (nlm, nang)
    '''
    nlm = len(lm)
    nang = len(x)
    y_kj = np.zeros((nlm,nang))
    for ilm, (l,m) in enumerate(lm):
        y_kj[ilm] = compute_real_spherical_harmonics(l, m, theta=np.arccos(z), phi=np.arctan2(y,x))
    return y_kj
```



- *Spherical expansion on electron density*


Using the orthogonality of spherical harminics: 
$
& rho_(l m) (r) = integral.double_Omega Y_(l m) (Omega) rho(bvec(r)) sin theta d Omega \ 
& rho_(i k) = sum_j underbrace((rho_(i j) w^a_(i j)), rho^a_(i j)) underbrace((Y_(k j) w_j)^T, "row-wise product") \ 
$

$rho_(i j)$ is eletron density values on the grid of atom A in matrix form, multiplied by A's Becke weight on its grid $w^a_(i j)$ we get $rho^a_(i j)$. By doing matrix multiplication with spherical basis the electron density expansion coefficients are obtained $rho_(i k)$ `shape (nrad, nlm)`

```py
rho_ik[iatom] = (watom[iatom] * rho[iatom]).reshape((nrad, nang)) @ (y_kj * wleb).T
```






- *Coulomb potential on atomic grid*

After getting $rho^a_(i k)$, by using fdd and solving linear equations the Coulomb potential expansion coefficients on atomic grid are obtained: $rho^a_(i k) -> u^a_(i k)$. 

But the Coulomb potential obtained at this point is "local". As shown in @fig6, the target point is not on the grid of atom a, but the contribution of atom a to the Coulomb potential at this point is greater than that of atom b, which generates this point. If the contribution of atom a is not considered, the Coulomb potential at this point will be greatly underestimated. Since this point is not a grid point of atom a, the Coulomb potential of atom a at this point cannot be obtained by solving the Poisson equation, and therefore interpolation calculation is required. 

#figure(
  image("pics/fig6.svg", width: 50%), 
  caption: "Numerically build Hartree potential on the grids.", 
  supplement: "Figure"
)<fig6>

$u_(l m) (r)$ is a function of scalar $r$, intuitively we would like to do interpolation on $r$. However $r$ is non-uniform grid, and the large the $r$ the larger the spacing, which leads to very inaccurate interpolated results at large $r$. 

#figure(
  image("pics/fig7.svg", width: 50%), 
  supplement: "Figure", 
  caption: "Cubic spline interpolation."
)

- *Construction Coulomb potential on the whole grid via interpolation*

As suggested in @becke19882, the interpolation is performed on uniform grid ${z}$. Note that $z$ is the same for all atoms while $r$ is scaled differently for different elements, this is another advantage of using Chebyshev radial grid. 

```py
for ilm, (l, _) in ennumerate(lm):
    u_ik[iatom] = poisson_solver(rho_ik[iatom, :, ilm], l, rrad[iatom], qn) # (nrad, nlm)
    u_ik_splines[iatom].append(CubicSpline(z, u_ik[iatom][:, ilm]))
```

After constructing an interpolator for each pair of $(l,m)$, interpolation needs to be performed on *all* grids. Thus, the constructed Coulomb potential is no longer "local" but "global". Therefore, it is necessary to calculate the distances from all atoms to all grids:

```py
def _get_distance(grid_global:np.ndarray, xyz:np.ndarray) -> np.ndarray:
'''
    Arguments
    ---------
    grid_global : np.ndarray (natom,nradxnang,3)
    xyz : np.ndarray (natom,3)

    Return
    ------
    dist : np.ndarray (natom,natomxnradxnang)
    '''
    grid_total = np.concatenate(grid_global) # (natomxnradxnang,3)
    natom = len(xyz)
    ngrid = len(grid_total)
    dist  = np.zeros((natom,ngrid))
    for iatom in range(natom):
        _disp = grid_total[:,:3] - xyz[iatom]
        _dist = np.linalg.norm(_disp, axis=1)
        dist[iatom]  = _dist

    return dist

 zz = [(nrad + 1) / np.pi * np.arccos((dist[iatom] - rms[iatom]) / (dist[iatom] + rms[iatom])) for iatom in range(natom)] # r = r(z) z = z(r) for interpolation
```

and calculate the contribution of atom i to all grids via cubic spline interpolation: 
```py
for ilm, (l, _) in enumerate(lm):
  u_ik[iatom] = poisson_solver(rho_ik[iatom], ...) # (nrad, nlm)
  u_ik_splines[iatom].append(CubicSpline(z, u_ik[iatom][:,ilm])) # (nlm, )
  u_ik_interp[iatom][ilm] = u_ik_splines[iatom][ilm](zz[iatom]) # (nlm, ngrid)
```

Re-constructing the global Coulomb potential generated by atom i by inverse spherical expansion: $U^a_"ngrid" = sum_(k) u^"a interp"_(k,"ngrid") \/ "dist"^a dot.c y^a_(k,"ngrid")$, where $y^a_(k, "ngrid")$ needs to compute solid angles from atom i to all grids, and then compute the spherical harmonics values of all these solid angles at all $(l,m)$ pairs: 
```py
def _get_solid_angle(grid_global:np.ndarray, xyz:np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    '''
    Arguments
    ---------
    grid_global : np.ndarray (natom,nradxnang,3)
    xyz : np.ndarray (natom,3)

    Return
    ------
    theta : np.ndarray (natom,natomxnradxnang)
    phi : np.ndarray (natom,natomxnradxnang)
    '''
    grid_total = np.concatenate(grid_global) # (natomxnradxnang,3)
    natom = len(xyz)
    ngrid = len(grid_total)
    theta = np.zeros((natom,ngrid))
    phi   = np.zeros((natom,ngrid))
    for iatom in range(natom):
        _disp = grid_total[:,:3] - xyz[iatom]
        _disp_normalized = _disp / _dist[:,None]
        _theta = np.arccos(_disp_normalized[:,2])
        _phi = np.arctan2(_disp_normalized[:,1], _disp_normalized[:,0])
        theta[iatom] = _theta
        phi[iatom]   = _phi

    return theta, phi

def _get_spherical_vals(lm:list[(int,int)], theta:np.ndarray, phi:np.ndarray) -> np.ndarray:
    '''
    Arguments
    ---------
    lm : list[(int,int)]
    theta : np.ndarray (natom, natomxnradxnang)
    phi : np.ndarray (natom, natomxnradxnang)

    Return
    ------
    ylm : np.ndarray (natom, nlm, natomxnradxnang)
    '''
    natom, ngrid = theta.shape
    nlm = len(lm)
    ylm = np.zeros((natom,nlm,ngrid)) # (natom,natomxnradxnang)

    for iatom in range(natom):
        for ilm, (l,m) in enumerate(lm):
            ylm[iatom,ilm] = rsh(l, m, theta=theta[iatom], phi=phi[iatom])
    
    return ylm
```

And finally: 
```py 
u_ij = np.zeros((ngrid,))

for iatom in range(natom):
    for ilm, (l, _) in enumerate(lm):
        u_ik[iatom] = poisson_solver(rho_ik[iatom], ...) # (nrad, nlm)
        u_ik_splines[iatom].append(CubicSpline(z, u_ik[iatom][:,ilm])) # (nlm, )
        u_ik_interp[iatom][ilm] = u_ik_splines[iatom][ilm](dist) # (nlm, ngrid)
    u_ij += np.sum((u_ik_interp[iatom] / dist[iatom]) * ylm[iatom], axis=0)
```

The interpolated Coulomb potential loops over all atoms, and takes contributions from all atoms on the whole grid into consideration. Thus what we get here is the *global Coulomb potential*. 

`u_ij shape (ngrid = natomxnradxnang, )`. Do so to facilitate the subsequent evaluation of operator matrix element.

- *Coulomb potential matrix element*

@atomicoperator computes matrix element on atomic grid. To conveniently calculate matrix elements on the whole molecular grid, flatten all quantities represented in matrix form: `u_ij -> shape (natomxnradxnang,)`. 

The calculation of matrix element involves calculation of basis function values on the grid. Organize basis function values as multi-dimensional array: ${phi.alt^mu_(a,i j) | mu=1,2,dots.c,"nbf" | a=1,2,dots.c,"natom"}$. Flatten it atom-wisely as a matrix with shape `(nbf x ngrid)` ：${phi.alt_(mu, "ngrid") | mu=1,2,dots.c,"nbf"}$. 

And combine all kinds of weights and Jacobian and flatten them: 
```py
# watom (natom, natgrid)
mweights = [np.outer(rrad[iatom]**2 * wrad[iatom], wleb).ravel() * watom[iatom] for iatom in range(natom)] # (natom, natgrid)
mweights = np.flatten(mweights) # (natomxnradxnang, )
```

And the Coulomb potential matrix can be evaluated by: 
$
& J_(mu nu) = diraccc(phi.alt_mu, U, phi.alt_nu) \ 
& cases(phi.alt_mu -> "(ngrid, )" , 
        phi.alt_nu -> "(ngridm, )",
        U -> "u_ij (ngrid, )" ) \ 

& diraccc(phi.alt_mu, U, phi.alt_nu) = integral r^2 d\r integral.double phi.alt_mu (bvec(r)) U(bvec(r)) phi.alt_nu (bvec(r)) sin theta d Omega \ 
& approx sum_a^"natom" w_a sum_i^"nrad" (r^a_i)^2 w^a_i sum_j^"nang" phi.alt_mu (r^a_i, Omega_j) phi.alt_nu (r^a_i, Omega_j) U(r^a_i, Omega_j) w_j \ 
& = sum_a sum_i sum_j underbrace({w_a w^a_i (r^a_i)^2 w_j}, "mweights") phi.alt^mu_(a,i j) phi.alt^nu_(a, i j) U_(a, i j) \ 
& = sum_a sum_(i j)^"natgrid" W_(a, "natgrid") phi.alt^mu_(a, "natgrid") phi.alt^nu_(a, "natgrid") U_(a, "natgrid") \ 
& = sum_(a i j)^"ngrid" W_"ngrid" phi.alt^mu_"ngrid" phi.alt^nu_"ngrid" U_"ngrid" => J_(mu nu) \ 
$

One possible implementation: 
```py
for mu in range(nbf):
    for nu in range(nbf):
        Juv[mu,nu]  = np.sum(aos_vals[mu].ravel() * u_ij * aos_vals[nu].ravel() * mweights)
```

And we have numerically calculated the Coulomb potential operator matrix. 





== Numeric quadrature of exchange-correlation functional

Exchange-correlation operator matrix elements are usually evaluated numerically via quadrature, as the mathematical expression for exchange-correlation functional is complicated and the integral over real space is pretty hard to compute analytically. 

Exchange-correlation functional is functional of electron density, so we calculate the electron density on the grid: 
$
& rho(bold(r)) = sum_i^"nocc" 2 |f_i|^2 \ 
& = sum_i^"nocc" 2 sum_mu C_(mu i) phi.alt_mu sum_nu C_(nu i) phi.alt_nu \ 
& = sum_(mu nu) (sum_i^"nocc" 2 C_(mu i) C_(nu i)) phi.alt_mu phi.alt_nu \ 
& = sum_(mu nu) P_(mu nu) phi.alt_mu phi.alt_nu \ 
$

in matrix form: 
```py
rho = np.zeros((ngrid,))
for mu in range(nbf):
    for nu in range(nbf):
        rho += 2 * aos_vals[mu] * Puv[mu,nu] * aos_vals[nu]
```

Integrating electron density over the whole real space reproduce the total number of electrons: $"ne" = int rho(bold(r)) d^3 bold(r) = $ ```py np.sum(rho * mweights)```. 

By using the same trick as before, the exchange-correlation potential matrix element and energy density matrix element are obtained. Once electron density on the grid is calculated, the exchange-correlation potential and energy density can be computed by calling `libxc`: 
```cpp
#include "xtensor.hpp"
#include "xc.h"
// link library `-lxc`

xc_func_type funcx, funcc;
xc_func_init(&funcx, X_id, XC_UNPOLARIZED); // init exchange
xc_func_init(&funcc, C_id, XC_UNPOLARIZED); // init correlation

auto rho = compute_rho();

xt::xtensor<double, 1> vx = xt::zeros<double>({ngrid}); // exchange potential
xt::xtensor<double, 1> ex = xt::zeros<double>({ngrid}); // exchange energy density
xt::xtensor<double, 1> vc = xt::zeros<double>({ngrid});
xt::xtensor<double, 1> ec = xt::zeros<double>({ngrid});

xc_lda_vxc(&funcx, ngrid, rho.data(), vx.data());
xc_lda_exc(&funcx, ngrid, rho.data(), ex.data());
xc_lda_vxc(&funcc, ngrid, rho.data(), vc.data());
xc_lda_exc(&funcc, ngrid, rho.data(), ec.data());

for (auto ibf=0; ibf < nbf; ++ibf) {
    for (auto jbf=0; jbf < nbf; ++jbf) {
        Kuv(ibf,  jbf) = xt::sum(aos_vals_mu * vx * aos_vals_nu * mweights)();
        Cuv(ibf,  jbf) = xt::sum(aos_vals_mu * vc * aos_vals_nu * mweights)();
        Exuv(ibf, jbf) = xt::sum(aos_vals_mu * ex * aos_vals_nu * mweights)();
        Ecuv(ibf, jbf) = xt::sum(aos_vals_mu * ec * aos_vals_nu * mweights)();
        }
    }
}

xc_func_end(&funcx); // free exchange
xc_func_end(&funcc); // free correlation
```


== Fock matrix

Summing up all operator matrices gives Fock matrix: 
$
  F = underbrace(T + V_"ext", "analytic") + underbrace(U, "analytic\nor\nnumeric") + underbrace(K + C, "numeric")
$

The energy components can be calculated by using density matrix: $dirac(hat(O)) = "Tr"{P O}$. And the energy decomposition: 
```py
kin_e         = np.trace(Puv @ tij.T)
ext_e         = np.trace(Puv @ vij.T)
hartree_e     = np.trace(Puv @ Juv.T) / 2 # double counting
exchange_e    = np.trace(Puv @ Exuv.T) # energy density matrix DOES NOT appear in Fock
correlation_e = np.trace(Puv @ Ecuv.T)
etot = kin_e + ext_e + hartree_e + exchange_e + correlation_e + nuc_repulsion
```



== Transformation from Cartesian basis function to pure basis function

The ordering of Cartesian basis function is typically lexicographic, for details see @Cartbf. While there are two commonly used ordering conventions for the pure basis functions: 
1. HORTON/ORCA convention. see #link("https://theochem.github.io/horton/2.1.1/tech_ref_gaussian_basis.html?highlight=basis")[this link]

  $m = 0, 1, -1, 2, -2, dots.c, l, -l$. $C$ stands for cos $because m > 0$, $S$ stands for sin $because m < 0$. 

  #table(
align: left,
columns: 2, 
[$l$], [Ordering], 
[0], [$C_00$], 
[1], [$C_10 =z, C_11 = x, S_11 = y$], 
[2], [$C_20, C_21, S_21, C_22, S_22$], 
[3], [$C_30, C_31, S_31, C_32, S_32, C_33, S_33$], 
[4], [$C_40, C_41, S_41, C_42, S_42, C_43, S_43, C_44, S_44$], 
)


2. libint/CCA convention
  $m = -l, -l+1, dots.c, l-1, l$. 

The pure basis function is linear combination of Cartesian basis function. In practical calculations, the results (such as integrals, function values, etc.) are first obtained in Cartesian form, and then by matrix multiplication with transform matrix the results in pure form are obtained. 

Transformation matrices under the HORTON convention can be found at  #link("https://github.com/theochem/iodata/blob/main/tools/harmonics.py"). 

For example, matrix form of basis function in Cartesian form on the grid: `aos_vals (nbf, ngrid)`. The transformation matrix is block diagonal: 

#figure(
  image("pics/fig8.svg", width: 50%), 
  supplement: "Figure", 
  caption: "Cartesian to pure transformation matrix."
)

Matrix multiplying with `aos_vals` gives pure basis function values on the grid: `spsdf 21 -> 17`. 

The same idea applies to the operator matrices. Note that the transformation matrices are highly sparse. 

Below lists some transformation matrices up to `lmax=3` in HORTON convention: 
$
& C_00 = 1 \

& vec(
    C_10,  
    C_11,  
    S_11) = mat(
    dot.c , dot.c ,     1 ;
        1 , dot.c , dot.c ;
    dot.c ,     1 , dot.c ;
     ) dot.c 
vec(
    x, 
    y, 
    z
) \

& vec(
    C_20 , 
    C_21 , 
    S_21 , 
    C_22 , 
    S_22) = mat(
    - 1/2     , dot.c , dot.c , - 1/2        , dot.c ,     1 ;
    dot.c     , dot.c , 1     , dot.c        , dot.c , dot.c ;
    dot.c     , dot.c , dot.c , dot.c        , 1     , dot.c ;
    sqrt(3)/2 , dot.c , dot.c , - sqrt(3)/2  , dot.c , dot.c ;
    dot.c     , 1     , dot.c , dot.c        , dot.c , dot.c ;
) dot.c
vec(
    x\x , 
    x\y , 
    x\z , 
    y\y , 
    y\z , 
    z\z
) \

& vec(
    C_30 , 
    C_31 , 
    S_31 , 
    C_32 , 
    S_32 , 
    C_33 , 
    S_33
) = mat(
    dot.c       , dot.c         , - (3 sqrt(5))/10 , dot.c                  , dot.c , dot.c      , dot.c        , - (3 sqrt(5))/10 , dot.c      ,     1 ;
    - sqrt(6)/4 , dot.c         , dot.c            , - sqrt(30)/20          , dot.c , sqrt(30)/5 , dot.c        , dot.c            , dot.c      , dot.c ;
    dot.c       , - sqrt(30)/20 , dot.c            , dot.c                  , dot.c , dot.c      , - sqrt(6)/4  , dot.c            , sqrt(30)/5 , dot.c ;
    dot.c       , dot.c         , sqrt(3)/2        , dot.c                  , dot.c , dot.c      , dot.c        , - sqrt(3)/2      , dot.c      , dot.c ;
    dot.c       , dot.c         , dot.c            , dot.c                  , 1     , dot.c      , dot.c        , dot.c            , dot.c      , dot.c ;
    sqrt(10)/4  , dot.c         , dot.c            , - (3 sqrt(2))/4        , dot.c , dot.c      , dot.c        , dot.c            , dot.c      , dot.c ;
    dot.c       , (3 sqrt(2))/4 , dot.c            , dot.c                  , dot.c , dot.c      , - sqrt(10)/4 , dot.c            , dot.c      , dot.c ;
) dot.c
vec(
    x\x\x , 
    x\x\y , 
    x\x\z , 
    x\y\y , 
    x\y\z , 
    x\z\z , 
    y\y\y , 
    y\y\z , 
    y\z\z , 
    z\z\z
)
$










= *Appendix*

- *密度矩阵*

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



$
& "ne" = sum_i^"occ" 2 int |psi_i (bold(r))|^2 d^3 bold(r) \ 
& = sum_i^"occ" int 2 (sum_mu |C_(mu,i) phi.alt_mu|) (sum_nu C_(nu,i)  phi.alt_nu) d^3 bold(r) \ 
& = sum_mu sum_nu 2 sum_i^"occ" C_(mu,i) C_(nu,i) int phi.alt_mu^dagger phi.alt_nu d^3 bold(r) \ 
& = sum_(mu nu) P_(mu,nu) chevron.l phi.alt_mu | phi.alt_nu chevron.r \ 
& <=> sum_(mu nu) P_(mu,nu) S_(mu,nu) \ 
& equiv "Tr"{P S^T} \ 
$

直接計算"密度矩陣" $P_(mu,nu)$ 的跡不會給出電子數, 因爲*基函數不正交*. 





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
  ...


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
