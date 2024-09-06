# Aula 8 -- 05-09-2024

## Nova equação

### Formulação Forte

$$ \begin{cases}
  - \alpha \; u_{xx}(x) + \beta \; u(x) + \gamma \; u_x(x) = f(x)
    &\qquad x \in (0, 1) \\
  u(0) = u(1) = 0 \\
\end{cases} $$

### Formulação Fraca

Dado um $v$ num "espaço infinito":

$$
  - \alpha \; \int_0^1{ u_{xx}(x) \; v(x) \; dx}
  + \beta \; \int_0^1{ u(x) \; v(x) \; dx}
  + \gamma \; \int_0^1{ u_x(x) \; v(x) \; dx}
  = \int_0^1{ f(x) \; v(x) \; dx}
$$ $$
  - \alpha \; \left[ \left( u_x(x) \; v(x) \right|_0^1 - \int_0^1{ u_x(x) \; v_x(x) \; dx} \right]
  + \beta \; \int_0^1{ u(x) \; v(x) \; dx}
  + \gamma \; \int_0^1{ u_x(x) \; v(x) \; dx}
  = \int_0^1{ f(x) \; v(x) \; dx}
$$ $$
  - \alpha \; \left[ - \int_0^1{ u_x(x) \; v_x(x) \; dx} \right]
  + \beta \; \int_0^1{ u(x) \; v(x) \; dx}
  + \gamma \; \int_0^1{ u_x(x) \; v(x) \; dx}
  = \int_0^1{ f(x) \; v(x) \; dx}
$$ $$
  \alpha \; \int_0^1{ u_x(x) \; v_x(x) \; dx}
  + \beta \; \int_0^1{ u(x) \; v(x) \; dx}
  + \gamma \; \int_0^1{ u_x(x) \; v(x) \; dx}
  = \int_0^1{ f(x) \; v(x) \; dx}
$$

Usando notação mágica:

* $$
  \kappa(u, v) =
    \alpha \; \int_0^1{ u_x(x) \; v_x(x) \; dx} + \beta \; \int_0^1{ u(x) \; v(x) \; dx} + \gamma \; \int_0^1{ u_x(x) \; v(x) \; dx}
$$

* $$
  (u, v) = \int_0^1{ f(x) \; v(x) \; dx}
$$

O sistema fica
$$ \begin{cases}
  \kappa(u, v) = (f, v)
\end{cases} $$

### Problema aproximado (Galerkin)

Reduzindo o espaço de funções para dimensão $m$.
Temos
$$ \begin{cases}
  \kappa(u^h, v^h) = (f, v^h)
\end{cases} $$

### Forma matriz-vetor

Podemos usar as bases
$\{ \varphi_1, \varphi_2, \dots, \varphi_m \}$,
temos que $u^h(x) = \sum_{j=1}^m{ c_j \; \varphi_j(x) }$

E adicionar mais linhas no sistema.
Uma linha $i \in \{ 1, \dots, m \}$,
usamos $v^h_i = \varphi_i$:

$$ \begin{cases}
  \kappa(u^h, v^h) = (f, v^h)
\end{cases} $$ $$ \begin{cases}
  \kappa(\sum_{j=1}^m{ c_j \; \varphi_j(x) }, \varphi_i) = (f, \varphi_i)
\end{cases} $$ $$ \begin{cases}
  \sum_{j=1}^m{ c_j \; \kappa(\varphi_j, \varphi_i) } = (f, \varphi_i)
\end{cases} $$

Notação:

* $K_{i,j} = \kappa(\varphi_j, \varphi_i)$
* $F_i = (f, \varphi_i)$
  (igual ao antigo)

Então, temos um sistema $m \times m$:
$$
  \mathbb{K} \; \mathbb{C} = \mathbb{F}
$$

* Calculando $K_{i,i}$

$$
  K_{i,i}
$$ $$
  \kappa(\varphi_i, \varphi_i)
$$ $$
  \alpha \; \int_0^1{ {\varphi_i}_x(x) \; {\varphi_i}_x(x) \; dx}
  + \beta \; \int_0^1{ \varphi_i(x) \; \varphi_i(x) \; dx}
  + \gamma \; \int_0^1{ {\varphi_i}_x(x) \; \varphi_i(x) \; dx}
$$ $$
  \alpha \; \left[
    \int_{x_{i-1}}^{x_i}{ {\varphi_i}_x(x) \; {\varphi_i}_x(x) \; dx}
    + \int_{x_i}^{x_{i+1}}{ {\varphi_i}_x(x) \; {\varphi_i}_x(x) \; dx}
  \right]
  + \beta \; \left[
    \int_{x_{i-1}}^{x_i}{ \varphi_i(x) \; \varphi_i(x) \; dx}
    + \int_{x_i}^{x_{i+1}}{ \varphi_i(x) \; \varphi_i(x) \; dx}
  \right]
  + \gamma \; \left[
    \int_{x_{i-1}}^{x_i}{ {\varphi_i}_x(x) \; \varphi_i(x) \; dx}
    + \int_{x_i}^{x_{i+1}}{ {\varphi_i}_x(x) \; \varphi_i(x) \; dx}
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \int_{x_{i-1}}^{x_i}{ \frac{1}{h_i} \; \frac{x - x_{i-1}}{h_i} \; dx}
    + \int_{x_i}^{x_{i+1}}{ \frac{-1}{h_{i+1}} \; \frac{x_i - x}{h_{i+1}} \; dx}
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{{h_i}^2} \; \int_{x_{i-1}}^{x_i}{ (x - x_{i-1}) \; dx}
    + \frac{1}{{h_{i+1}}^2} \; \int_{x_i}^{x_{i+1}}{ (x - x_i) \; dx}
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{{h_i}^2} \; \left( \frac{x^2}{2} - x_{i-1} \; x \right|_{x_{i-1}}^{x_i}
    + \frac{1}{{h_{i+1}}^2} \; \left( \frac{x^2}{2} - x_i \; x \right|_{x_i}^{x_{i+1}}
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{2 \; {h_i}^2} \; \left( x^2 - 2 \; x_{i-1} \; x \right|_{x_{i-1}}^{x_i}
    + \frac{1}{2 \; {h_{i+1}}^2} \; \left( x^2 - 2 \; x_i \; x \right|_{x_i}^{x_{i+1}}
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{2 \; {h_i}^2} \; \left( ({x_i}^2 - 2 \; x_{i-1} \; x_i) - ({x_{i-1}}^2 - 2 \; x_{i-1} \; x_{i-1}) \right)
    + \frac{1}{2 \; {h_{i+1}}^2} \; \left( ({x_{i+1}}^2 - 2 \; x_i \; x_{i+1}) - ({x_i}^2 - 2 \; x_i \; x_i) \right)
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{2 \; {h_i}^2} \; \left( {x_i}^2 - 2 \; x_{i-1} \; x_i + {x_{i-1}}^2 \right)
    + \frac{1}{2 \; {h_{i+1}}^2} \; \left( {x_{i+1}}^2 - 2 \; x_i \; x_{i+1} + {x_i}^2 \right)
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{2 \; {h_i}^2} \; (x_i - x_{i-1})^2
    + \frac{1}{2 \; {h_{i+1}}^2} \; (x_{i+1} - x_i)^2
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[
    \frac{1}{2 \; {h_i}^2} \; {h_i}^2
    + \frac{1}{2 \; {h_{i+1}}^2} \; {h_{i+1}}^2
  \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma \; \left[ \frac{1}{2} + \frac{1}{2} \right]
$$ $$
  \alpha \; \left( \frac{1}{h_i} + \frac{1}{h_{i+1}} \right)
  + \beta \; \left( \frac{h_i}{3} + \frac{h_{i+1}}{3} \right)
  + \gamma
$$

* Calculando $K_{i,i+1}$

$$
  K_{i,i+1}
$$ $$
  \alpha \; \int_0^1{ {\varphi_i}_x(x) \; {\varphi_{i+1}}_x(x) \; dx}
  + \beta \; \int_0^1{ \varphi_i(x) \; \varphi_{i+1}(x) \; dx}
  + \gamma \; \int_0^1{ {\varphi_i}_x(x) \; \varphi_{i+1}(x) \; dx}
$$ $$
  \alpha \; \int_{x_i}^{x_{i+1}}{ {\varphi_i}_x(x) \; {\varphi_{i+1}}_x(x) \; dx}
  + \beta \; \int_{x_i}^{x_{i+1}}{ \varphi_i(x) \; \varphi_{i+1}(x) \; dx}
  + \gamma \; \int_{x_i}^{x_{i+1}}{ {\varphi_i}_x(x) \; \varphi_{i+1}(x) \; dx}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \int_{x_i}^{x_{i+1}}{ \frac{-1}{h_{i+1}} \; \frac{x - x_i}{h_{i+1}} \; dx}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{{h_{i+1}}^2} \; \int_{x_i}^{x_{i+1}}{ (x - x_i) \; dx}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{{h_{i+1}}^2} \; \left( \frac{x^2}{2} - x_i \; x \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; \left( x^2 - 2 \; x_i \; x \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; ( ({x_{i+1}}^2 - 2 \; x_i \; x_{i+1}) - ({x_i}^2 - 2 \; x_i \; x_i)
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; ({x_{i+1}}^2 - 2 \; x_i \; x_{i+1} + {x_i}^2)
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; (h_{i+1})^2
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  - \frac{\gamma}{2}
$$

* Calculando $K_{i+1,i}$

$$
  K_{i,i+1}
$$ $$
  \alpha \; \int_0^1{ {\varphi_i}_x(x) \; {\varphi_{i+1}}_x(x) \; dx}
  + \beta \; \int_0^1{ \varphi_i(x) \; \varphi_{i+1}(x) \; dx}
  + \gamma \; \int_0^1{ {\varphi_{i+1}}_x(x) \; \varphi_i(x) \; dx}
$$ $$
  \alpha \; \int_{x_i}^{x_{i+1}}{ {\varphi_i}_x(x) \; {\varphi_{i+1}}_x(x) \; dx}
  + \beta \; \int_{x_i}^{x_{i+1}}{ \varphi_i(x) \; \varphi_{i+1}(x) \; dx}
  + \gamma \; \int_{x_i}^{x_{i+1}}{ {\varphi_{i+1}}_x(x) \; \varphi_i(x) \; dx}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \int_{x_i}^{x_{i+1}}{ \frac{1}{h_{i+1}} \; \frac{x_{i+1} - x}{h_{i+1}} \; dx}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{{h_{i+1}}^2} \; \int_{x_i}^{x_{i+1}}{ (x - x_{i+1}) \; dx}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{{h_{i+1}}^2} \; \left( \frac{(x - x_{i+1})^2}{2} \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; \left( (x - x_{i+1})^2 \right|_{x_i}^{x_{i+1}}
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; \left( (x_{i+1} - x_{i+1})^2 - (x_i - x_{i+1})^2 \right)
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; \left( -(x_i - x_{i+1})^2 \right)
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \gamma \; \frac{-1}{2 \; {h_{i+1}}^2} \; (- {h_{i+1}}^2)
$$ $$
  \frac{- \alpha}{h_{i+1}}
  + \frac{\beta}{6} \; h_{i+1}
  + \frac{\gamma}{2}
$$
