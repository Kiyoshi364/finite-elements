# Aula 13 -- 24-09-2024

## (2) Cálculo de $U^0$ como a projeção de $L^2$ de $u_0$

Seja $U^0 \in V_m$, com
$$
  (U^0 - u_0, v_h) = 0, \qquad \forall v_h \in V_m
$$

Fazendo
$U^0 = \sum_{j=0}^m C^0_j \; \varphi_j$
e
$v_h = \varphi_i$ para $i = 1, \dots, m$,
temos:
$$
  (\sum_{j=0}^m{ C^0_j \; \varphi_j - u_0 }, \varphi_i) = 0
$$ $$
  \sum_{j=0}^m{ (C^0_j \; \varphi_j - u_0, \varphi_i) } = 0
$$ $$
  \sum_{j=0}^m{ (C^0_j \; \varphi_j, \varphi_i) - (u_0, \varphi_i) } = 0
$$ $$
  \sum_{j=0}^m{ C^0_j \; (\varphi_j, \varphi_i) } - (u_0, \varphi_i) = 0
$$ $$
  \sum_{j=0}^m{ C^0_j \; (\varphi_j, \varphi_i) } = (u_0, \varphi_i)
$$

Então:
$$
  \mathbb{M} \; C^0 = \mathbb{F}
$$
com
$$
  \mathbb{M}_{i,j} = (\varphi_j, \varphi_i)
$$
e
$$
  \mathbb{F}_i = (u_0, \varphi_i)
$$

## (3) Cálculo de $U^0$ como a projeção de $H^1_0$ de $u_0$

Seja $U^0 \in V_m$, com
$$
  (\frac{d}{dx}(U^0 - u_0), \frac{d}{dx}v_h) = 0, \qquad \forall v_h \in V_m
$$

Fazendo
$U^0 = \sum_{j=0}^m C^0_j \; \varphi_j$
e
$v_h = \varphi_i$ para $i = 1, \dots, m$,
temos:
$$
  (\frac{d}{dx}(\sum_{j=0}^m{ C^0_j \; \varphi_j - u_0 }), \varphi_i) = 0
$$ $$
  (\frac{d}{dx}(\sum_{j=0}^m{ C^0_j \; \varphi_j }) - \frac{d}{dx}u_0, \varphi_i) = 0
$$ $$
  (\frac{d}{dx}(\sum_{j=0}^m{ C^0_j \; \varphi_j }), \varphi_i) - (\frac{d}{dx}u_0, \varphi_i) = 0
$$ $$
  \sum_{j=0}^m{ C^0_j \; (\frac{d}{dx}\varphi_j, \varphi_i) } - (\frac{d}{dx}u_0, \varphi_i) = 0
$$ $$
  \sum_{j=0}^m{ C^0_j \; (\frac{d}{dx}\varphi_j, \varphi_i) } = (\frac{d}{dx}u_0, \varphi_i)
$$

Então:
$$
  \mathbb{A} \; C^0 = \mathbb{B}
$$
com
$$
  \mathbb{A}_{i,j} = (\frac{d}{dx}\varphi_j, \varphi_i)
$$
e
$$
  \mathbb{B}_i = (\frac{d}{dx}u_0, \varphi_i)
$$

## (4) Cálculo de $U^0$ usando o operador $\kappa$

Seja $U^0 \in V_m$, com
$$
  \kappa(U^0 - u_0, v_h) = 0, \qquad \forall v_h \in V_m
$$

Fazendo
$U^0 = \sum_{j=0}^m C^0_j \; \varphi_j$
e
$v_h = \varphi_i$ para $i = 1, \dots, m$,
temos:
$$
  \kappa(\sum_{j=0}^m{ C^0_j \; \varphi_j - u_0 }, \varphi_i) = 0
$$ $$
  \kappa(\sum_{j=0}^m{ C^0_j \; \varphi_j }, \varphi_i) = \kappa(u_0, \varphi_i) = 0
$$ $$
  \sum_{j=0}^m{ C^0_j \; \kappa(\varphi_j, \varphi_i) } = \kappa(u_0, \varphi_i)
$$

Então:
$$
  \mathbb{K} \; C^0 = \mathbb{B}
$$
com
$$
  \mathbb{K}_{i,j} = \kappa(\varphi_j, \varphi_i)
$$
e
$$
  \mathbb{B}_i = \kappa(u_0, \varphi_i)
$$
