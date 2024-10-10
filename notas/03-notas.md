# Aula 18 -- 10-10-2024

## Espaço 2D Mudança de Variável

$$
  x_1(\xi) = \frac{h_1}{2}(\xi_1 + 1) + p_1^e
$$ $$
  x_2(\xi) = \frac{h_2}{2}(\xi_2 + 1) + p_2^e
$$ $$
  |J| = \left|\begin{matrix}
    \frac{\partial}{\partial \xi_1} x_1(\xi)
    &
    \frac{\partial}{\partial \xi_2} x_1(\xi)
    \\
    \frac{\partial}{\partial \xi_1} x_2(\xi)
    &
    \frac{\partial}{\partial \xi_2} x_2(\xi)
  \end{matrix}\right|
  = \left|\begin{matrix}
    \frac{h_1}{2} & 0 \\
    0 & \frac{h_2}{2}
  \end{matrix}\right|
  = \frac{h_1 \; h_2}{4}
$$

### Relação entre $\frac{\partial}{\partial x_i} \varphi^e_a(x(\xi))$ e $\frac{\partial}{\partial \xi_i} \phi_a(\xi)$

$$
    \varphi^e_a(x(\xi)) = \phi_a(\xi)
$$ $$
    \frac{\partial}{\partial \xi_i} \varphi^e_a(x(\xi))
    = \frac{\partial}{\partial \xi_i} \phi_a(\xi)
$$ $$
    \frac{\partial}{\partial x_1} \varphi^e_a(x(\xi))
    \; \frac{\partial}{\partial \xi_i} x_1(\xi)
    + \frac{\partial}{\partial x_2} \varphi^e_a(x(\xi))
    \; \frac{\partial}{\partial \xi_i} x_2(\xi)
    = \frac{\partial}{\partial \xi_i} \phi_a(\xi)
$$ $$
    \frac{\partial}{\partial x_i} \varphi^e_a(x(\xi))
    \; \frac{\partial}{\partial \xi_i} x_i(\xi)
    = \frac{\partial}{\partial \xi_i} \phi_a(\xi)
$$ $$
    \frac{\partial}{\partial x_i} \varphi^e_a(x(\xi))
    \; \frac{h_i}{2}
    = \frac{\partial}{\partial \xi_i} \phi_a(\xi)
$$ $$
    \frac{\partial}{\partial x_i} \varphi^e_a(x(\xi))
    = \frac{2}{h_i}
    \; \frac{\partial}{\partial \xi_i} \phi_a(\xi)
$$

### Definição da matriz global $K$

$$
    K_{i,j} =
    \alpha \; \int_\Omega{
      \frac{\partial}{\partial x_1}\varphi_j(x)
      \; \frac{\partial}{\partial x_1}\varphi_i(x)
      \; dx}
    + \alpha \; \int_\Omega{
      \frac{\partial}{\partial x_2}\varphi_j(x)
      \; \frac{\partial}{\partial x_2}\varphi_i(x)
      \; dx}
    + \beta \; \int_\Omega{ \varphi_j(x) \; \varphi_i(x) \; dx}
$$

### Definição da matriz global $K^e$

$$
    K^e_{a,b}
$$ $$
    \alpha \; \int_{\Omega^e}{
      \frac{\partial}{\partial x_1}\varphi^e_b(x)
      \; \frac{\partial}{\partial x_1}\varphi^e_a(x)
      \; dx}
    + \alpha \; \int_{\Omega^e}{
      \frac{\partial}{\partial x_2}\varphi^e_b(x)
      \; \frac{\partial}{\partial x_2}\varphi^e_a(x)
      \; dx}
    + \beta \; \int_{\Omega^e}{
      \varphi^e_b(x) \; \varphi^e_a(x)
      \; dx}
$$ $$
    \alpha \; \int_R{
      \frac{\partial}{\partial x_1}\varphi^e_b(x(\xi))
      \; \frac{\partial}{\partial x_1}\varphi^e_a(x(\xi))
      \; |J|
      \; d\xi}
    + \alpha \; \int_R{
      \frac{\partial}{\partial x_2}\varphi^e_b(x(\xi))
      \; \frac{\partial}{\partial x_2}\varphi^e_a(x(\xi))
      \; |J|
      \; d\xi}
    + \beta \; \int_R{
      \varphi^e_b(x(\xi)) \; \varphi^e_a(x(\xi))
      \; |J|
      \; d\xi}
$$ $$
    \alpha \; \int_R{
      \frac{2}{h_1} \; \frac{\partial}{\partial \xi_1}\phi_b(\xi)
      \; \frac{2}{h_1} \frac{\partial}{\partial \xi_1}\phi_a(\xi)
      \; \frac{h_1 \; h_2}{4}
      \; d\xi}
    + \alpha \; \int_R{
      \frac{2}{h_2} \; \frac{\partial}{\partial \xi_2}\phi_b(\xi)
      \; \frac{2}{h_2} \; \frac{\partial}{\partial \xi_2}\phi_a(\xi)
      \; \frac{h_1 \; h_2}{4}
      \; d\xi}
    + \beta \; \int_R{
      \phi_b(\xi) \; \phi_a(\xi)
      \; \frac{h_1 \; h_2}{4}
      \; d\xi}
$$ $$
    \alpha \; \frac{h_2}{h_1} \; \int_R{
      \frac{\partial}{\partial \xi_1}\phi_b(\xi)
      \; \frac{\partial}{\partial \xi_1}\phi_a(\xi)
      \; d\xi}
    + \alpha \; \frac{h_1}{h_2} \; \int_R{
      \frac{\partial}{\partial \xi_2}\phi_b(\xi)
      \; \frac{\partial}{\partial \xi_2}\phi_a(\xi)
      \; d\xi}
    + \beta \; \frac{h_1 \; h_2}{4} \; \int_R{
      \phi_b(\xi) \; \phi_a(\xi)
      \; d\xi}
$$
