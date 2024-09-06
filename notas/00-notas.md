# Aula 1 -- 13-08-2024

## Problema estacionário 1D (Diferenças Finitas)

* Problema de Valor de Contorno (Formulação Forte)

$$ \begin{cases}
  - \alpha \; u_{xx}(x) + \beta \; u(x) = f(x)
    \quad\text{, em } \Omega = (0, 1) \\
  u(0) = u(1) = 0 \\
  \alpha > 0, \beta \ge 0
\end{cases} $$

Forma discreta:

$$ \begin{cases}
  - \alpha \; u_i'' + \beta \; u_i = f_i
    \quad\text{, } i = 1, 2, \dots, N \\
  \alpha > 0, \beta \ge 0
\end{cases} $$

Notação:
$$ \begin{cases}
  u_0 = u_1 = 0 \\
  u_i'' = \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2} \\
  f_i = f(x_i) \\
  u_i = u(x_i) \\
  h = x_{i+1} - x_i = \frac{x_{N+1} - x_0}{N+1}
\end{cases} $$

Com isso temos um sistema linear:

$$ \begin{cases}
  - \alpha \; \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2} + \beta \; u_i = f_i
    \quad\text{, } i = 1, 2, \dots, N \\
  u_0 = u_1 = 0 \\
  \alpha > 0, \beta \ge 0
\end{cases} $$

Contas:
$$
  - \alpha \; \frac{u_{i+1} - 2u_i + u_{i-1}}{h^2} + \beta \; u_i = f_i
$$ $$
  - \alpha \; u_{i+1} + \alpha \; 2u_i - \alpha \; u_{i-1} + \beta \; u_i \; h^2 = f_i \; h^2
$$ $$
  - \alpha \; u_{i+1} + (2\alpha + \beta \; h^2) \; u_i - \alpha \; u_{i-1} = f_i \; h^2
$$
Matrix template:
$$ \begin{bmatrix}
  2\alpha + \beta \; h^2 & -\alpha & 0 & 0 & 0 & \cdots & 0 \\
  -\alpha & 2\alpha + \beta \; h^2 & -\alpha & 0 & 0 & \cdots & 0 \\
  0 & -\alpha & 2\alpha + \beta \; h^2 & -\alpha & 0 & \cdots & 0 \\
  0 & \ddots & \ddots & \ddots & \ddots & \cdots & \vdots \\
  0 & \cdots & 0 & -\alpha & 2\alpha + \beta \; h^2 & -\alpha & 0 \\
  0 & \cdots & 0 & 0 & -\alpha & 2\alpha + \beta \; h^2 & -\alpha \\
  0 & \cdots & 0 & 0 & 0 & -\alpha & 2\alpha + \beta \; h^2 \\
\end{bmatrix} \cdot u = \begin{bmatrix}
  f_1 \; h^2 + \alpha \; u_0 \\ f_2 \; h^2 \\ \\ \vdots \\ \\ f_{N-1} \; h^2 \\ f_N \; h^2 + \alpha \; u_{N+1} \\
\end{bmatrix} $$

### Exemplos

* Seja $u(x) = x \; (x-1)$:

$$
  - \alpha \; u_{xx}(x) + \beta \; u(x) = f(x)
  \Rightarrow
  f(x) = - 2 \; \alpha + \beta \; x \; (x-1)
$$
Logo,
$$ \begin{cases}
  - \alpha \; u_{xx}(x) + \beta \; u(x) = - 2 \; \alpha + \beta \; x \; (x-1)
    \quad\text{, em } \Omega = (0, 1) \\
  u(0) = u(1) = 0
\end{cases} $$

* Seja $u(x) = sen (\pi \; x)$:
$$
  f(x) = \pi^2 \; \alpha \; sen(\pi \; x) + \beta \; sen(\pi \; x)
$$
Logo,

$$ \begin{cases}
  - \alpha \; u_{xx}(x) + \beta \; u(x) = \pi^2 \; \alpha \; sen(\pi \; x) + \beta \; sen(\pi \; x)
    \quad\text{, em } \Omega = (0, 1) \\
  u(0) = u(1) = 0
\end{cases} $$

* Cálculo do Erro Relativo

$$
  E_{rel} = \frac{\|u - u^h\|}{\|u\|}
$$
onde $u$ e $u^h$ são vetores,
$u_i = u(x_i)$ e $u^h_i = u^h(x_i)$,
$u^h$ é a solução aproximada.

#### Branch (Condicionamento)

Condicionamento está relacionado com
"se mexer um pouco nos valores (adicionar erros),
quanto o resultado muda?"
$$
  cond(A) = \|A\| \; \|A^{-1}\|
$$
É provado que $cond(A) \ge 1$.

Se $cond(A)$ é perto de $1$ é melhor
(muda menos com os erros).

## Formulação Fraca

Espaço das soluções
$$
  H = \{
    u \text{ é suficientemente suave } : u(0) = u(1) = 0
  \}
$$
Espaço das funções de teste (funções peso)
$$
  V = \{
    v \text{ é suficientemente suave } : v(0) = v(1) = 0
  \}
$$

> **Nota**: Nesse caso $H = V$

Seja $v \in V$ qualquer:
$$
  - \alpha \; \int_0^1{ u_{xx}(x) \; v(x) \; dx } + \beta \; \int_0^1{ u(x) \; v(x) \; dx } = \int_0^1{ f(x) \; v(x) \; dx }
$$ $$
  - \alpha \; \left[ \left( v(x) \; u_x(x) \right|_0^1 - \int_0^1{ u_x(x) \; v_x(x) \; dx } \right] + \beta \; \int_0^1{ u(x) \; v(x) \; dx } = \int_0^1{ f(x) \; v(x) \; dx }
$$ $$
  - \alpha \; \left[ - \int_0^1{ u_x(x) \; v_x(x) \; dx } \right] + \beta \; \int_0^1{ u(x) \; v(x) \; dx } = \int_0^1{ f(x) \; v(x) \; dx }
$$
Logo
$$ \begin{cases}
  \alpha \; \int_0^1{ u_x(x) \; v_x(x) \; dx } + \beta \; \int_0^1{ u(x) \; v(x) \; dx } = \int_0^1{ f(x) \; v(x) \; dx } \\
  u(0) = u(1) = 0 \\
  v(0) = v(1) = 0 \\
  \alpha > 0, \beta \ge 0
\end{cases} $$

Notação:

* $$
  a(u, v) =
  \alpha \; \int_0^1{ u_x(x) \; v_x(x) \; dx } + \beta \; \int_0^1{ u(x) \; v(x) \; dx }
$$
* $$
  (u, v) =
  \int_0^1{ u(x) \; v(x) \; dx }
$$
Com essa notação o sistema fica:
$$ \begin{cases}
  a(u, v) = (f, v) \\
  u \in U \\
  v \in V \\
  \alpha > 0, \beta \ge 0
\end{cases} $$

* Forma discreta:

Espaço das soluções:
$$
  H^h = \{
    u \text{ é suficientemente suave } : u(0) = u(1) = 0
  \}
$$
Espaço das funções de teste (funções peso)
$$
  V^h = \{
    v \text{ é suficientemente suave } : v(0) = v(1) = 0
  \}
$$
Ambos tem dimensão finitas $m$ e base $\{ \Phi_1, \Phi_2, \dots, \Phi_m\}$

Problema aproximado:
Determinar $u^h \in H^h$ tal que
$$ \begin{cases}
  a(u^h, v^h) = (f, v) \\
  v^h \in V^h
\end{cases} $$

Logo, $u^h(x) = \sum_{j=1}^m c_j \; \Phi_j(x)$
$$ \begin{cases}
  a(\sum_{j=1}^m{ c_j \; \Phi_j }, v^h) = (f, v^h) \\
  v^h \in V^h
\end{cases} $$ $$ \begin{cases}
  \sum_{j=1}^m{ c_j \; a(\Phi_j, v^h) } = (f, v^h) \\
  v^h \in V^h
\end{cases} $$

Para cada $v^h = \Phi_i$ com $i = 1, 2, \dots, m$:
$$ \begin{cases}
  \sum_{j=1}^m{ c_j \; a(\Phi_j, \Phi_i) } = (f, \Phi_i) \\
\end{cases} $$

Notação:

* $K_{i,j} = a(\Phi_j, \Phi_i)$
* $F_i = (f, \Phi_i)$

Então temos o sistema linear $m \times m$:
$$
  \mathbb{K} \; \mathbb{C} = \mathbb{F}
$$

## Resumo da aula

$$
S \Rightarrow W \Rightarrow G \Rightarrow M
$$

* $S$: Strong   System
* $W$: Weak     System
* $G$: Galerkin System (Discretização)
* $M$: Matrix   System

---
# Aula 2 -- 15-08-2024

## Definindo base das funções $\Phi$

Usando os $\Phi_i$, definido da seguinte forma:
$$
  \Phi_i(x) := \begin{cases}
    0 &\quad\text{, if } 0 \le x < x_{i-1} \\
    \frac{x - x_{i-1}}{x_i - x_{i-1}} &\quad\text{, if } x_{i-1} \le x < x_i \\
    1 - \frac{x - x_i}{x_{i+1} - x_i} &\quad\text{, if } x_i \le x < x_{i+1} \\
    0 &\quad\text{, if } x \le x_{i+1} \le 1 \\
  \end{cases}
$$
Ou usando a notação $h_i = x_i - x_{i-1}$
$$
  \Phi_i(x) := \begin{cases}
    0 &\quad\text{, if } 0 \le x < x_{i-1} \\
    \frac{x - x_{i-1}}{h_i} &\quad\text{, if } x_{i-1} \le x < x_i \\
    \frac{x_{i+1} - x}{h_{i+1}} &\quad\text{, if } x_i \le x < x_{i+1} \\
    0 &\quad\text{, if } x \le x_{i+1} \le 1 \\
  \end{cases}
$$

Temos:

* Como (quando $|i - j| \ge 2$
  $\Phi_i(x) \; \Phi_j(x) = 0$),
  $\mathbb{K}$ é tridiagonal
* Se $$
  \Phi_j(x_i) = \begin{cases}
    1 &\quad\text{, if } i = j \\
    0 &\quad\text{, if } i \ne j \\
  \end{cases}
$$ então $$
  u^h(x_i)
  = \sum_{j = 1}^m c_j \; \Phi_j(x_i)
  = c_i
$$

Contas para $K_{i, i}$:

$$
  K_{i, i}
$$ $$
  a(\Phi_i, \Phi_i)
$$ $$
  \alpha \; \int_0^1{ \frac{d \Phi_i}{dx}(x) \; \frac{d \Phi_i}{dx}(x) \; dx }
    + \beta \; \int_0^1{ \Phi_i(x) \; \Phi_i(x) \; dx }
$$ $$
  \alpha \; \int_{x_{i-1}}^{x_{i+1}}{ \frac{d \Phi_i}{dx}(x) \; \frac{d \Phi_i}{dx}(x) \; dx }
    + \beta \; \int_{x_{i-1}}^{x_{i+1}}{ \Phi_i(x) \; \Phi_i(x) \; dx }
$$ $$
  \alpha \; \int_{x_{i-1}}^{x_i}{ \frac{d \Phi_i}{dx}(x) \; \frac{d \Phi_i}{dx}(x) \; dx }
    + \alpha \; \int_{x_i}^{x_{i+1}}{ \frac{d \Phi_i}{dx}(x) \; \frac{d \Phi_i}{dx}(x) \; dx }
    + \beta \; \int_{x_{i-1}}^{x_i}{ \Phi_i(x) \; \Phi_i(x) \; dx }
    + \beta \; \int_{x_i}^{x_{i+1}}{ \Phi_i(x) \; \Phi_i(x) \; dx }
$$ $$
  \alpha \; \int_{x_{i-1}}^{x_i}{ \frac{1}{h_i} \; \frac{1}{h_i} \; dx }
    + \alpha \; \int_{x_i}^{x_{i+1}}{ \frac{-1}{h_{i+1}} \; \frac{-1}{h_{i+1}} \; dx }
    + \beta \; \int_{x_{i-1}}^{x_i}{ \frac{x - x_{i-1}}{h_i} \; \frac{x - x_{i-1}}{h_i} \; dx }
    + \beta \; \int_{x_i}^{x_{i+1}}{ \frac{x_{i+1} - x}{h_{i+1}} \; \frac{x_{i+1} - x}{h_{i+1}} \; dx }
$$ $$
  \frac{\alpha}{{h_i}^2} \; h_i
    + \frac{\alpha}{{h_{i+1}}^2} \; h_{i+1}
    + \frac{\beta}{{h_i}^2} \; \left( \frac{(x - x_{i-1})^3}{3} \right|_{x_i}^{x_{i-1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \left( \frac{(-1) \; (x_{i+1} - x)^3}{3} \right|_{x_{i+1}}^{x_i}
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
    + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
$$

Contas para $K_{i, i+1}$ e $K_{i+1, i}$:
$$
  K_{i, i+1} = K_{i+1, i}
$$ $$
  a(\Phi_i, \Phi_{i+1})
$$ $$
  \alpha \; \int_0^1{ \frac{d \Phi_i}{dx}(x) \; \frac{d \Phi_{i+1}}{dx}(x) \; dx }
    + \beta \; \int_0^1{ \Phi_i(x) \; \Phi_{i+1}(x) \; dx }
$$ $$
  \alpha \; \int_{x_i}^{x_{i+1}}{ \frac{d \Phi_i}{dx}(x) \; \frac{d \Phi_{i+1}}{dx}(x) \; dx }
    + \beta \; \int_{x_i}^{x_{i+1}}{ \Phi_i(x) \; \Phi_{i+1}(x) \; dx }
$$ $$
  \alpha \; \int_{x_i}^{x_{i+1}}{ \left( \frac{-1}{h_{i+1}} \right) \; \left( \frac{1}{h_{i+1}} \right) \; dx }
    + \beta \; \int_{x_i}^{x_{i+1}}{ \left( \frac{x_{i+1} - x}{h_{i+1}} \right) \; \left( \frac{x - x_i}{h_{i+1}} \right) \; dx }
$$ $$
  \frac{-\alpha}{{h_{i+1}}^2} \; h_{i+1}
    + \beta \; \int_{x_i}^{x_{i+1}}{ \left( \frac{x_{i+1} - x}{h_{i+1}} \right) \; \left( \frac{x - x_i}{h_{i+1}} \right) \; dx }
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \int_{x_i}^{x_{i+1}}{ (x_{i+1} - x) \; (x - x_i) \; dx }
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \int_{x_i}^{x_{i+1}}{ (-x^2 + (x_{i+1} + x_i)x - x_{i+1}x_i) \; dx }
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \left( -\frac{x^3}{3} + (x_{i+1} + x_i)\frac{x^2}{2} - x_{i+1}x_ix \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \left( -\frac{{x_{i+1}}^3 - {x_i}^3}{3} + (x_{i+1} + x_i)\frac{{x_{i+1}}^2 - {x_i}^2}{2} - x_{i+1}x_i(x_{i+1} - x_i) \right)
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \left( -\frac{{x_{i+1}}^3 - {x_i}^3}{3} + \frac{{x_{i+1}}^3 + {x_{i+1}}^2x_i - x_{i+1}{x_i}^2 - {x_i}^3}{2} - {x_{i+1}}^2x_i + x_{i+1}{x_i}^2 \right)
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{{h_{i+1}}^2} \; \left( -\frac{{x_{i+1}}^3 - {x_i}^3}{3} + \frac{{x_{i+1}}^3 + {x_{i+1}}^2x_i - x_{i+1}{x_i}^2 - {x_i}^3}{2} - {x_{i+1}}^2x_i + x_{i+1}{x_i}^2 \right)
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{6{h_{i+1}}^2} \; ( -2{x_{i+1}}^3 + 2{x_i}^3 + 3{x_{i+1}}^3 + 3{x_{i+1}}^2x_i - 3x_{i+1}{x_i}^2 - 3{x_i}^3 - 6{x_{i+1}}^2x_i + 6x_{i+1}{x_i}^2 )
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{6{h_{i+1}}^2} \; ( {x_{i+1}}^3 - {x_i}^3 - 3{x_{i+1}}^2x_i + 3x_{i+1}{x_i}^2 )
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{6{h_{i+1}}^2} \; (x_{i+1} - x_i)^3
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{6{h_{i+1}}^2} \; {h_{i+1}}^3
$$ $$
  \frac{-\alpha}{h_{i+1}}
    + \frac{\beta}{6} \; h_{i+1}
$$

Contas para $F_i$:
$$
  F_i
$$ $$
  \int_0^1{ f(x) \; \Phi_i(x) dx }
$$ $$
  \int_{x_{i-1}}^{x_{i+1}}{ f(x) \; \Phi_i(x) dx }
$$ $$
  \int_{x_{i-1}}^{x_i}{ f(x) \; \Phi_i(x) dx }
    + \int_{x_i}^{x_{i+1}}{ f(x) \; \Phi_i(x) dx }
$$ $$
  \int_{x_{i-1}}^{x_i}{ f(x) \; \frac{x - x_{i-1}}{h_i} dx }
    + \int_{x_i}^{x_{i+1}}{ f(x) \; \frac{x_i - x}{h_{i+1}} dx }
$$ $$
  \frac{1}{h_i} \; \int_{x_{i-1}}^{x_i}{ f(x) \; (x - x_{i-1}) dx }
    + \frac{1}{h_{i+1}} \int_{x_i}^{x_{i+1}}{ f(x) \; (x_i - x) dx }
$$ $$
  \frac{1}{h_i} \; \left[ \left( f^*(x) \; (x - x_{i-1}) \right|_{x_{i-1}}^{x_i} - \int_{x_{i-1}}^{x_i}{ f^*(x) dx } \right]
    + \frac{1}{h_{i+1}} \left[ \left( f^*(x) \; (x_i - x) \right|_{x_i}^{x_{i+1}} + \int_{x_i}^{x_{i+1}}{ f^*(x) dx } \right]
$$ $$
  \frac{1}{h_i} \; \left[ (f^*(x_i) \; (x_i - x_{i-1}) - f^*(x_{i-1}) \; (x_{i-1} - x_{i-1})) - \left( f^{*2}(x) \right|_{x_{i-1}}^{x_i} \right]
    + \frac{1}{h_{i+1}} \left[ (f^*(x_{i+1}) \; (x_i - x_{i+1}) - f^*(x_i) \; (x_i - x_i)) + \left( f^{*2}(x) \right|_{x_i}^{x_{i+1}} \right]
$$ $$
  \frac{1}{h_i} \; \left[ f^*(x_i) \; h_i - ( f^{*2}(x_{i-1}) - f^{*2}(x_i) ) \right]
    + \frac{1}{h_{i+1}} \left[ - f^*(x_{i+1}) \; h_{i+1} + ( f^{*2}(x_{i+1}) - f^{*2}(x_i) ) \right]
$$ $$
  \frac{1}{h_i} \; \left[ f^*(x_i) \; h_i + ( f^{*2}(x_i) - f^{*2}(x_{i-1}) ) \right]
    + \frac{1}{h_{i+1}} \left[ - f^*(x_{i+1}) \; h_{i+1} + ( f^{*2}(x_{i+1}) - f^{*2}(x_i) ) \right]
$$ $$
  f^*(x_i) - f^*(x_{i+1})
    + \frac{f^{*2}(x_i) - f^{*2}(x_{i-1})}{h_i}
    + \frac{f^{*2}(x_{i+1}) - f^{*2}(x_i)}{h_{i+1}}
$$

Então:
$$
  K_{i, i} =
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
    + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
$$ $$
  K_{i, i+1} = K_{i, i+1} =
  -\frac{\alpha}{h_{i+1}} + \frac{\beta}{6} \; h_{i+1}
$$ $$
  F_i =
  ???
$$

Passo a passo:

1. Montar $\mathbb{K}$
2. Montar $\mathbb{F}$
3. Resolver $\mathbb{K} \; \mathbb{C} = \mathbb{F}$

### Exemplos

* Seja $u(x) = x + \frac{e^{-x} - e^x}{e^1 - e^{-1}}$ e
  $f(x) = x$:

$$
  F_i
$$ $$
  \int_{x_{i-1}}^{x_i}{ x \; \frac{x - x_{i-1}}{h_i} dx }
    + \int_{x_i}^{x_{i+1}}{ x \; \frac{x_i - x}{h_{i+1}} dx }
$$ $$
  \frac{1}{h_i} \; \int_{x_{i-1}}^{x_i}{ x^2 - x \; x_{i-1} dx }
    + \frac{1}{h_{i+1}} \; \int_{x_i}^{x_{i+1}}{ x_i \; x - x^2 dx }
$$ $$
  \frac{1}{h_i} \; \left( \frac{x^3}{3} - \frac{x^2 \; x_{i-1}}{2} \right|_{x_{i-1}}^{x_i}
    + \frac{1}{h_{i+1}} \; \left( \frac{x_i \; x^2}{2} - \frac{x^3}{3} \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{1}{6 \; h_i} \; \left( 2 \; x^3 - 3 \; x^2 \; x_{i-1} \right|_{x_{i-1}}^{x_i}
    + \frac{1}{6 \; h_{i+1}} \; \left( 3 \; x_i \; x^2 - 2 \; x^3 \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{1}{6 \; h_i} \; ( 2 \; {x_i}^3 - 3 \; {x_i}^2 \; x_{i-1} - 2 \; {x_{i-1}}^3 + 3 \; {x_{i-1}}^3 )
    + \frac{1}{6 \; h_{i+1}} \; ( 3 \; x_i \; {x_{i+1}}^2 - 2 \; {x_{i+1}}^3 - 3 \; {x_i}^3 + 2 \; {x_i}^3 )
$$ $$
  \frac{1}{6 \; h_i} \; ( 2 \; {x_i}^3 - 3 \; {x_i}^2 \; x_{i-1} + \; {x_{i-1}}^3 )
    - \frac{1}{6 \; h_{i+1}} \; ( 2 \; {x_{i+1}}^3 - 3 \; x_i \; {x_{i+1}}^2 + {x_i}^3 )
$$ TODO $$
  \frac{1}{6 \; h_i} \; (2 \; x_i + x_{i-1}) \; {h_i}^2
    + \frac{1}{6 \; h_{i+1}} \; ({x_{i+1}}^3 - 3 \; x_{i+1} \; {x_i}^2 + 2 \; {x_i}^3)
$$ TODO $$
  \frac{1}{6 \; h_i} \; (2 \; x_i + x_{i-1}) \; {h_i}^2
    + \frac{1}{6 \; h_{i+1}} \; (x_{i+1} + 2 \; x_i) \; {h_{i+1}}^2
$$ $$
  \frac{h_i}{6} \; (2 \; x_i + x_{i-1})
    + \frac{h_{i+1}}{6} \; (x_{i+1} + 2 \; x_i)
$$

Se $h = h_i = h_{i+1}$:
$$
  \frac{h}{6} \; (2 \; x_i + x_i - h)
    + \frac{h}{6} \; (x_i + h + 2 \; x_i)
$$ $$
  \frac{h}{6} \; (x_i + x_i - h + x_i + h + 2 \; x_i)
$$ $$
  \frac{h}{6} \; (6 \; x_i)
$$ $$
  x_i \; h
$$

* Seja $u(x) = x \; (x-1)$ e
  $f(x) = - 2 \; \alpha + \beta \; x \; (x-1)$:

TODO

* Seja $u(x) = sen (\pi \; x)$ e
  $f(x) = \pi^2 \; \alpha \; sen(\pi \; x) + \beta \; sen(\pi \; x)$:

TODO

---
# Aula 5 -- 27-08-2024

Cálculo do erro:
$$
  \|u - u_h\|^2
$$ $$
  \int_0^1{ (u(x) - u_h(x))^2 \; dx}
$$ $$
  \int_0^1{ \left(u(x) - \sum_{i=1}^N{ c_i \; \varphi_i(x) }\right)^2 \; dx}
$$ $$
  \sum_{i=1}^{N+1}{ \int_{(i-1)\;h}^{i\;h}{ \left(u(x) - \sum_{j=1}^N{ c_j \; \varphi_j(x) }\right)^2 \; dx} }
$$ $$
  \sum_{i=1}^{N+1}{ \int_{(i-1)\;h}^{i\;h}{ (u(x) - (c_{i-1} \; \varphi_{i-1}(x) + c_i \; \varphi_i(x)))^2 \; dx} }
$$

Trocando variável:
$$ \begin{cases}
  (\frac{2}{h} \; x - 1 \Leftrightarrow \xi)
    \leftrightarrow
    (x \Leftrightarrow \frac{h}{2} \; (\xi + 1)) \\
  \frac{2}{h} \; dx \Leftrightarrow d\xi \\
  c_0 = c_{N+1} = 0 \\
  \phi_1(\xi) = \frac{1 - \xi}{2} \\
  \phi_2(\xi) = \frac{1 + \xi}{2} \\
  x(\xi, i) = h \; \left(\frac{\xi + 1}{2} + (i-1)\right) \\
\end{cases} $$

Temos:
$$
  \sum_{i=1}^{N+1}{ \frac{h}{2} \; \int_{-1}^1{ (u(x(\xi, i)) - (c_{i-1} \; \phi_1(\xi) + c_i \; \phi_2(\xi)))^2 \; d\xi} }
$$ $$
  \frac{h}{2} \; \sum_{i=1}^{N+1}{ \int_{-1}^1{ (u(x(\xi, i)) - c_{i-1} \; \phi_1(\xi) - c_i \; \phi_2(\xi))^2 \; d\xi} }
$$

Aplicando quadratura gaussiana:
$$
  \frac{h}{2} \; \sum_{i=1}^{N+1}{ \sum_{j=0}^{n_{pg}}{ W_j \; (u(x(P_j, i)) - c_{i-1} \; \phi_1(P_j) - c_i \; \phi_2(P_j))^2 } }
$$

---
# Aula 6 -- 29-08-2024

* Lembretes
$$ \begin{array}{lcl}
  h = x^e_2 - x^e_1 \\
  \varphi_1^e(x) = \frac{x^e_2 - x}{h} &\quad&
  \varphi_2^e(x) = \frac{x - x^e_1}{h} \\
  x(\xi, e) = \frac{h}{2} \; (\xi + 1) + x^e_1 &&
  dx = \frac{h}{2} \; d\xi \\
  \phi_1(\xi) = \frac{1 - \xi}{2} &&
  \phi_2(\xi) = \frac{1 + \xi}{2} \\
\end{array} $$

* Matriz local $K^e$, para $a, b \in \{1, 2\}$

$$
  K^e_{a,b}
  =
  \alpha \; \int_{x^e_1}^{x^e_2}{ \frac{d \; \varphi^e_a}{dx}(x) \; \frac{d \; \varphi^e_b}{dx}(x) \; dx}
  + \beta \; \int_{x^e_1}^{x^e_2}{ \varphi^e_a(x) \; \varphi^e_b(x) \; dx}
$$

Realizando troca de variável para $\xi \in [-1, 1]$:
$$
  \alpha \; \int_{-1}^1{ \frac{d \; \varphi^e_a}{dx}(x(\xi, e)) \; \frac{d \; \varphi^e_b}{dx}(x(\xi, e)) \; \frac{h}{2} \; d\xi}
  + \beta \; \int_{-1}^1{ \varphi^e_a(x(\xi, e)) \; \varphi^e_b(x(\xi, e)) \; \frac{h}{2} \; d\xi}
$$ $$
  \alpha \; \frac{h}{2} \; \int_{-1}^1{ \frac{d \; \varphi^e_a}{dx}(x(\xi, e)) \; \frac{d \; \varphi^e_b}{dx}(x(\xi, e)) \; d\xi}
  + \beta \; \frac{h}{2} \; \int_{-1}^1{ \varphi^e_a(x(\xi, e)) \; \varphi^e_b(x(\xi, e)) \; d\xi}
$$ $$
  \alpha \; \frac{h}{2} \; \int_{-1}^1{ \frac{d \; \varphi^e_a}{dx}(x(\xi, e)) \; \frac{d \; \varphi^e_b}{dx}(x(\xi, e)) \; d\xi}
  + \beta \; \frac{h}{2} \; \int_{-1}^1{ \phi_a(\xi) \; \phi_b(\xi) \; d\xi}
$$

Para trocar variável de $\frac{d \; \varphi_a^e}{dx}(x(\xi, e))$ em função de $\phi_a(\xi)$:

$$
  \frac{d \; \varphi_a^e}{dx}(x(\xi, e)) \; \frac{x(\xi, e)}{dx} = \frac{d \; \phi_a}{d\xi}(\xi)
$$ $$
  \frac{d \; \varphi_a^e}{dx}(x(\xi, e)) \; \frac{h}{2} = \frac{d \; \phi_a}{d\xi}(\xi)
$$ $$
  \frac{d \; \varphi_a^e}{dx}(x(\xi, e)) = \frac{2}{h} \; \frac{d \; \phi_a}{d\xi}(\xi)
$$

Usando isso em $K^e_{a, b}$:
$$
  \alpha \; \frac{h}{2} \; \int_{-1}^1{ \frac{2}{h} \; \frac{d \; \phi_a}{d\xi}(\xi) \; \frac{2}{h} \; \frac{d \; \phi_b}{d\xi}(\xi) \; d\xi}
  + \beta \; \frac{h}{2} \; \int_{-1}^1{ \phi_a(\xi) \; \phi_b(\xi) \; d\xi}
$$ $$
  \alpha \; \frac{2}{h} \; \int_{-1}^1{ \frac{d \; \phi_a}{d\xi}(\xi) \; \frac{d \; \phi_b}{d\xi}(\xi) \; d\xi}
  + \beta \; \frac{h}{2} \; \int_{-1}^1{ \phi_a(\xi) \; \phi_b(\xi) \; d\xi}
$$

* Vetor local $F^e$, para $a, b \in \{1, 2\}$

$$
  F^e_a
  = \int_{x^e_1}^{x^e_2}{ f(x) \; \varphi^e_a(x) \; dx}
  = \int_{-1}^1{ f(x(\xi, e)) \; \varphi^e_a(x(\xi, e)) \; \frac{h}{2} \; d\xi}
$$ $$
  F^e_a =
  \frac{h}{2} \; \int_{-1}^1{ f(x(\xi, e)) \; \phi_a(\xi) \; d\xi}
$$
