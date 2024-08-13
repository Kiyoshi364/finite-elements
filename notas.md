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
  f_1 \; h^2 \\ \\ \\ \vdots \\ \\ \\ f_N \; h^2 \\
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
  a(\sum_{j=1}^m{ c_j \; \Phi_j }, v^h) = (f, v) \\
  v^h \in V^h
\end{cases} $$ $$ \begin{cases}
  \sum_{j=1}^m{ c_j \; a(\Phi_j, v^h) } = (f, v) \\
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

* $S$: Strong System
* $W$: Weak   System
* $G$: Galerk System (Discretização)
* $M$: Matrix System
