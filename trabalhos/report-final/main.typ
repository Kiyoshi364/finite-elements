#import "template.typ": template, fakepar

#let qquad = math.wide
#let varphi = sym.phi.alt
#let mathblock(numbered: true, body) = math.equation(
    block: true,
    numbering: if numbered { "(1)" } else { none },
    body,
)

#let images-folder = "images/"

#let title = [
    Usando iteradores e valores pre-computados
    #linebreak()
    para economizar recusos
    na construção de matrizes
    #linebreak()
    para Método de Elementos Finitos
    #linebreak()
    com a Linguagem de Programação Julia
]
#let author = [Daniel K. Hashimoto V. de Andrade --- 124259224]

#let problema(body) = {
    let indent = 0em
    set par(
        first-line-indent: indent,
    )
    set text(style: "italic")
    body
}

#let todo(body) = {
    text(fill: rgb("#ff0000"))[TODO:]
    body
}

#show: template.with(
    indent-first: false,
    title: title,
    author: author,
)

= Introdução

- elementos finitos
- performance
- implementaçoes
- nossa proposta
- sumário

= Especificação do Problema e o Método de Elementos Finitos

Nessa seção, vamos
formalizar o problema teórico;
revisar como
o método de elementos finitos
aproxima a solução,
chegando a formulação matricial;
e depois
discutir como as matrizes são construídas.
O objetivo dessa seção
é criar um contexto para
que seja possível compreender melhor
a melhoria proposta.

== Especificação do Problema

Nesse trabalho estamos aplicando elementos finitos
no seguinte problema:

#problema[
    Dados
    uma função $f : accent(Omega, macron) -> RR$
    e constantes $alpha > 0, beta >= 0 : RR$,
    encontre
    uma função $u : accent(Omega, macron) -> RR$
    tal que
    #mathblock()[$
        cases(
            - alpha Delta u(arrow(x)) + beta u(arrow(x))
                = f(arrow(x))
                \,&qquad x in Omega,
            u(arrow(x)) = 0
                \,&qquad x in Gamma,
        )
    $] <problema-forte>
    para algum sub-espaço $Omega subset RR^2$,
    com $Gamma$ fronteira de $Omega$ e
    $accent(Omega, macron) = Omega union Gamma$.
]

Obtemos a versão fraca do @problema-forte[Problema]
substituindo a equação principal pela seguinte:
#mathblock()[$
    alpha integral_Omega
        nabla u(arrow(x)) dot nabla v(arrow(x))
        d Omega
    + beta integral_Omega
        u(arrow(x)) v(arrow(x))
        d Omega
    =  integral_Omega
        f(arrow(x)) v(arrow(x))
        d Omega
$] <eq:fraca-ext>
onde $v : accent(Omega, macron) -> RR$
é qualquer função suficentemente suave
tal que $forall arrow(y) in Gamma, v(arrow(y)) = 0$.
Comumente, introduzimos operadores
para esconder a complexidade da equação
#mathblock()[$
    kappa(u, v) = (f, v)
$] <eq:fraca>
#fakepar()

O Método de Elementos Finitos
pode ser aplicado depois de fazermos duas discretizações:
uma no domínio $accent(Omega, macron)$,
e outra na dimensão do espaço das funções.
A primeira divide o domíno $accent(Omega, macron)$
em $N$ regiões contínuas $r_1, dots.h, r_N$.
Cada região contínua
não compartilha nenhum ponto com as outras regiões,
ou seja,
$forall 1 <= i, j <= N, i != j <=> r_i sect r_j = emptyset$.
Chamamos cada região $r_e$ ($1 <= e <= N$) de elemento $e$.

A segunda discretização
reduz o espaço de funções para dimensão $m$,
isso significa que vamos ter $m$
funções da base $phi_1, dots.h, phi_m$.
Com isso podemos representar a nossa função solução $u^h$
usando um vetor $c$:
$u^h (arrow(x)) = sum_(j=1)^m c_j phi_j (x)$.

Agora, para que @eq:fraca se mantenha verdade
para todas as funções do espaço discretizado,
basta que seja verdade
para as funções da base $phi_1, dots.h, phi_m$.
Transformando o problema infinito
em um sistema finito de equações da seguinte forma:
#mathblock()[$
    cases(
        c_1 kappa(phi_1, phi_1) &+ c_2 kappa(phi_2, phi_1)
            &+ dots.h.c + c_m kappa(phi_m, phi_1) &= (f, phi_1),
        c_1 kappa(phi_1, phi_2) &+ c_2 kappa(phi_2, phi_2)
            &+ dots.h.c + c_m kappa(phi_m, phi_2) &= (f, phi_2),
        med dots.v,
        c_1 kappa(phi_1, phi_m) &+ c_2 kappa(phi_2, phi_m)
            &+ dots.h.c + c_m kappa(phi_m, phi_m) &= (f, phi_m),
    )
$] <eq:sis:discreto>
#fakepar()

O @eq:sis:discreto[Sistema]
pode ser reescrito
em uma única equação matricial:
#mathblock()[$
    KK med c = FF
$] <eq:matriz>
onde $KK$ e $FF$ são definidas como:
$
    KK_(a,b) = kappa(phi_b, phi_a)
    qquad "e" qquad
    FF_a = (f, phi_a)
    , qquad qquad "com" 1 <= a, b <= m
$

== Sistema Matricial e Construção das Matrizes

No caso geral,
resolver um sistema matricial
com $r$ linhas e $c$ colunas é custoso
tanto para tempo quanto para memória.
Entretanto é comum fazer escolhas espertas
de funções da base,
que torne $0$ muitos coeficientes de $KK$,
tornando $KK$ uma matriz esparsa
e reduzindo a complexidade da solução do sistema matricial.

Com esse truque em mente,
a estratégia
para montar $KK$ e $FF$
já estabelecida na literatura
é calcular a contribuição de cada elemento
para essas matrizes.
Para isso,
fazemos todos os cálculos
de um certo elemento $e$
em um mesmo mundo canônico
refletimos sua contribuição nas $KK$ e $FF$.
Nós chamamos as _coisas_
relacionadas ao sistema matricial
de _coisas_ globais,
enquanto chamamos as _coisas_
relacionadas ao mundo canônico
de _coisas_ locais.

#todo[continue from here]

$
phi quad varphi
$

= A Melhoria Proposta e Implementação Realizada

#todo[
descrever a
proposta de melhoria na implementação,
mais especificamente na montagem das matrizes.
]

= Resultados

== Resultados locais

#let readcsv = name => csv(images-folder + name + ".csv")
#let csv_as_table(it) = align(
    center,
    grid(
        columns: (auto,) + ((1fr,) * (it.at(0).len() - 1)),
        row-gutter: 0.5em,
        align: (left,) + ((center,) * (it.at(0).len() - 1)),
        ..it.at(0).map(it => align(center, [*#it*])),
        ..it.slice(1).map(row => (
            row.at(0),
            ..row.slice(1).map(s => {
                let ss = s.split(".")
                if ss.len() == 1 {
                    ss.at(0)
                } else {
                    let len = ss.at(1).len()
                    ss.at(0)
                    [.]
                    ss.at(1).slice(0,calc.min(3,len))
                }
            })
        )).flatten()
    )
)

#csv_as_table(readcsv("smallvec"))
#csv_as_table(readcsv("smallmat"))

== Resultados Globais

#let magic_grid(names, files) = {
    let types = ("vec", "mat", "both")
    let len = names.len()
    grid(
        columns: (auto, 1fr, 1fr, 1fr),
        column-gutter: (1em),
        align: horizon,
        ..([], [Vetor], [Matriz], [Ambos]).map(body => align(center, body)),
        ..(array.range(len)).map(i => (
            names.at(i),
            types.map(type =>
                image(images-folder + type + "-" + files.at(i) + ".svg")
            ),
        )).flatten()
    )
}

#let graph_names = ("Tempo", "Memória", "Alocações")
#let graph_files = ("min_time", "min_memory", "min_alloc")

#magic_grid(graph_names, graph_files)

== Resultados Globais para cada Implementação

#let impls_names = ("Baseline", "Ref", "Iter", "Iter and Ref", "Bruno Carmo")
#let impls_files = ("baseline", "ref", "iter", "iter_ref", "bacarmo")

#magic_grid(impls_names, impls_files)

= Conclusão

== Trabalhos Futuros

- Aproveitar memória na LG
- Entender porque a construção da matriz está lenta

#lorem(200)
