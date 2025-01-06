#import "template.typ": template, fakepar

#let qquad = math.wide
#let varphi = sym.phi.alt
#let mathblock(numbered: true, body) = math.equation(
    block: true,
    numbering: if numbered { "(1)" } else { none },
    body,
)
#let julia(body) = raw(
    block: true,
    lang: "julia",
    body,
)
#let compare-julia(name1, name2, code1, code2) = figure(
    grid(
        columns: (1fr, 1fr),
        align: left,
        ..(name1, name2).map(it => align(center, it)),
        julia(code1), julia(code2),
    ),
    caption: [#(name1) and #(name2)],
    placement: auto,
    supplement: [Comparação de Código],
    kind: "julia-compare",
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
    text(fill: rgb("#ff0000"))[TODO#(if body != [] [: ] else [])]
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
em $N$ regiões contínuas $Omega_1, dots.h, Omega_N$.
Cada região contínua
não compartilha nenhum ponto com as outras regiões,
ou seja,
$forall 1 <= i, j <= N, i != j <=> Omega_i sect Omega_j = emptyset$.
Chamamos cada região $Omega_e$ ($1 <= e <= N$) de elemento $e$.

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
para cada componente dessas matrizes.
Para isso,
fazemos todos os cálculos
de um certo elemento $e$
em um mesmo mundo canônico,
gerando matrizes menores $KK^e$ e $FF^e$.
E então,
refletimos suas contribuição nas $KK$ e $FF$ originais.
Nós chamamos as _coisas_
relacionadas ao sistema matricial
de _coisas_ globais,
enquanto chamamos as _coisas_
relacionadas ao mundo canônico
de _coisas_ locais.

No mundo local,
ao invés de usar a base $phi_1, dots.h, phi_m$
vamos usar funções $varphi_1, dots.h, varphi_d : cal(R) -> RR$,
sendo $d$ o número de funções
em que um elemento tem contribuições,
ou seja,
o número de funções que retorna um valor diferente de $0$
quando recebe algum ponto do elemento.
$cal(R)$ é o espaço canônico,
um quadrado no intervalo aberto
$]-1, 1[ times ]-1, 1[$.
As funções $varphi_i$ não variam
de acordo ao elemento escolhido.

Dessa forma,
para um elemento $e$,
calculamos $KK^e$ e $FF^e$
da seguinte forma:
#mathblock()[$
    KK^e_(a,b) =
        alpha integral_(cal(R))
            frac(1, |J(arrow(xi))|)
            nabla varphi_b (arrow(xi))^T
            dot.c
            (H (arrow(xi))^T dot.c H (arrow(xi)))
            dot.c
            nabla varphi_a (arrow(xi))
            d arrow(xi)
        + beta integral_(cal(R))
            |J(arrow(xi))|
            varphi_b (arrow(xi))
            varphi_a (arrow(xi))
            d arrow(xi)
$] <Ke>
#mathblock()[$
    FF^e_a =
        integral_(cal(R))
            |J(arrow(xi))|
            f(x_e (arrow(xi)))
            varphi_a (arrow(xi))
            d arrow(xi)
$] <Fe>
onde
$x_e : cal(R) -> Omega_e$ é a função que mapeia
$arrow(xi)$
de volta o ponto global relativo que está no elemento $e$,
$x_(e 1)$ e $x_(e 2)$ mapeiam apenas
para a primeira e segunda coordenada, respectivamente,
do resultado de $x_e$.
$J(arrow(xi))$ é o determinante da matriz Jacobiana de $x_e$,
$H(arrow(xi))$ é uma matriz.
$J(arrow(xi))$ e $H(arrow(xi))$ são definidas como
$
    #let dd(i, j) = $frac(partial x_(e#i), partial arrow(xi)_#j)(arrow(xi))$
    J(arrow(xi)) = dd(1, 1) dd(2, 2) - dd(1, 2) dd(2, 1)
    qquad e qquad
    H(arrow(xi)) = mat(
        delim: "[",
        dd(2, 2), -dd(2, 1);
        -dd(1, 2), dd(1, 1);
    )
$

= A Melhoria Proposta e Implementação Realizada

Nessa seção,
descrevemos a proposta de melhoria
na implementação da montagem das matrizes,
os cálculos que podem ser reaproveitados,
e depois discutimos a implementação realizada.

== Aproveitamento dos Cálculos

As integrações nos cálculos de $KK^e$ (@Ke) e $FF^e$ (@Fe)
(e por sua vez $KK$ e $FF$)
geralmente
são feitas usando os pontos e pesos de gauss.
Mas todos os cálculos são feitos
no mesmo mundo local,
ou seja,
com os mesmos pontos de gauss e com as mesmas $varphi_i$'s.
Por isso, acabamos recalculando o mesmo valor diversas vezes.
Por exemplo,
os vetores $nabla varphi_i (arrow(xi))$ ($1 <= i <= d$)
são usados em todas as matrizes $KK^e$,
enquanto
os valores $varphi_i (arrow(xi))$ ($1 <= i <= d$)
são usados em todas matrizes $KK^e$ e vetores $FF^e$.

Os outros elementos da fórmula
dependem apenas de $x_e$.
Entretanto, dependendo do mapeamento escolhido,
podemos definir $x_e$ em termos de $varphi_i$.
Por exemplo, no caso específico de
uma discretização do espaço em quadrados
e base linear;
se $X^e$ é o vetor com as primeiras componentes
do quadrilátero $e$
e $Y^e$ é o vetor com as segundas componentes
do quadrilátero $e$,
podemos definir
#mathblock[$
    x_e (arrow(xi)) = vec(
        sum_(i=0)^d (X^e)_i dot.c varphi_i (arrow(xi)),
        sum_(i=0)^d (Y^e)_i dot.c varphi_i (arrow(xi))
    )
$] <x2xis>
#fakepar()

Considerando isso,
propormos pré-calcular os valores de
$varphi_i (arrow(xi))$ e
$nabla varphi_i (arrow(xi))$
para cada ponto de gauss em uma tabela.
Então durante cálculo das matrizes locais
podemos apenas olhar o valor na tabela
evitanto cálculos já realizados.

== Implementação

Para essa implementação,
usamos uma discretização de espaço
em quadriláteros possívelmente não regulares
e base linear.
Dessa forma,
temos $4$ funções locais da base ($d = 4$) e
definimos as funções $varphi_i$ e $nabla varphi_i$
da seguinte forma
#mathblock()[$
    #let def(q, w) = $frac((1 #q arrow(xi)_1) (1 #w arrow(xi)_2), 4)$
    varphi_1 (arrow(xi)) = def(-, -) qquad , qquad
    varphi_2 (arrow(xi)) = def(+, -) ,\
    varphi_3 (arrow(xi)) = def(+, +) qquad "e" qquad
    varphi_4 (arrow(xi)) = def(-, +)
$] <varphis>
#mathblock()[$
    #let de(q, w, i) = $#(if w == [+] [] else [#w]) frac(1 #q arrow(xi)_#i, 4)$
    #let def(q, w) = $vec(de(#w, #q, 2), de(#q, #w, 1))$
    nabla varphi_1 (arrow(xi)) = def(-, -) qquad , qquad
    nabla varphi_2 (arrow(xi)) = def(+, -) ,\
    nabla varphi_3 (arrow(xi)) = def(+, +) qquad "e" qquad
    nabla varphi_4 (arrow(xi)) = def(-, +)
$] <dvarphis>
Por dividir cada elemento em um quadrilátero,
usamos o mesmo mapeamento local-global
descrito na @x2xis.

Fizemos 2 implementações de construtores de $KK^e$ e $FF^e$
separadamente,
totalizando em 4.
Todas utilizando os valores pré-calculados em uma tabela.
As implementações do primeiro grupo de construtores locais
(chamadas de `Baseline`)
alocam a matriz internamente e retornam ela,
enquanto as do segundo grupo
(chamadas de `Ref`)
recebem memória já alocada e sobreescrevem ela
com o resultado.

Fizemos 4 implementações de construtores de $KK^e$ e $FF^e$
separadamente,
e depois as mesmas 4 implementações
calculando ambas as matrizes juntas,
totalizando em 12 implementações.
Essas implementações são divididas em 2 critérios:
reutilizar a memória das matrizes locais
entre elementos; e
utilizar uma abstração de iterador da linguagem Julia.
O grupo de implementações `Baseline`
não cumpre nenhum dos dois critérios.
`Ref` cumpre apenas o primeiro;
`Iter` cumpre apenas o segundo;
e finalmente,
`Iter and Ref` cumpre os dois.
Vamos mostrar o corpo da implementação
da construção de $FF$
para exemplificar as diferentes versões.

O primeiro critério
indica que a implementação reserva toda a memória necessária
antes de iterar os elementos,
reutiliza elas durate as iterações
e utiliza as versões `Ref` das implementações locais.
Esse critério está aí para ver
se o compilador de Julia
consegue perceber
que poderíamos estar reutilizando
a mesma memória.
Em @comparacao:base-ref,
percebemos que
as diferenças mais notáveis são
a alocação explícita antes do loop principal
e o uso de funções de mutação.

#compare-julia(
    `Baseline`, `Ref`,
    "
local F = fill(0.0, (m+1,))





for e in 1:N_e
    local LGe = view(LG, :, e)
    local Xe = view(X, LGe)
    local Ye = view(Y, LGe)

    local x2xis = x2xis_f(phis, Xe, Ye)
    local dx2xis = dx2xis_f(phi_derivs, Xe, Ye)

    local F_e = build_small_vec_2d(

        f,
        x2xis, dx2xis,
        phis, ws, gauss_n
    )

    local EQoLG_ = view(EQoLG, :, e)
    for i in 1:dim
        F[EQoLG_[i]] += F_e[i]
    end
end
F[begin:end-1]
    ", "
local F = fill(0.0, (m+1,))
local Fe = Ref(fill(0.0, (dim,)))
local x2xis = Ref(fill(0.0, (size(phis)[1:2]..., sdim)))
local dx2xis = Ref(fill(0.0, (size(phi_derivs)[1], size(phi_derivs)[1], dim)))
for e in 1:N_e
    local LGe = view(LG, :, e)
    local Xe = view(X, LGe)
    local Ye = view(Y, LGe)

    x2xis_f_ref!(x2xis, phis, Xe, Ye)
    dx2xis_f_ref!(dx2xis, phi_derivs, Xe, Ye)

    build_small_vec_2d_ref!(
        Fe,
        f,
        x2xis[], dx2xis[],
        phis, ws, gauss_n
    )

    local EQoLG_ = view(EQoLG, :, e)
    for i in 1:dim
        F[EQoLG_[i]] += Fe[][i]
    end
end
F[begin:end-1]
    ",
) <comparacao:base-ref>
#fakepar()

Já o segundo critério
indica que a implementação utiliza iteradores
para fazer os pré-calculos de $J(arrow(xi))$ e $H(arrow(xi))$
ao invés de pré-calcular explicitamente dentro do loop.
A ideia desse critério é tentar deixar corpo do loop mais compacto,
simples e com menos variáveis intermediárias visíveis;
ela também testa o quanto o Julia consegue otimizar
a abstração de iterador.
A principal diferença entre as implementações
(em @comparacao:base-iter)
é o uso do iterador `FiniteElementIter` e
o cálculo implícito das variáveis
`LGe`, `Xe`, `Ye`, `x2xis` e `dx2xis`.

#compare-julia(
    `Baseline`, `Iter`,
    "





local F = fill(0.0, (m+1,))
for e in 1:N_e
    local LGe = view(LG, :, e)
    local Xe = view(X, LGe)
    local Ye = view(Y, LGe)

    local x2xis = x2xis_f(phis, Xe, Ye)
    local dx2xis = dx2xis_f(dphis, Xe, Ye)

    local F_e = build_small_vec_2d(
        f,
        x2xis, dx2xis,
        phis, ws, gauss_n
    )

    local EQoLG_ = view(EQoLG, :, e)
    for i in 1:dim
        F[EQoLG_[i]] += F_e[i]
    end
end
F[begin:end-1]
    ", "
local iter = FiniteElementIter(
    X, Y,
    N_e, LG,
    phis, dphis,
)
local F = fill(0.0, (m+1,))
for (e, x2xis, dx2xis) in iter







    local F_e = build_small_vec_2d(
        f,
        x2xis, dx2xis,
        phis, ws, gauss_n
    )

    local EQoLG_ = view(EQoLG, :, e)
    for i in 1:dim
        F[EQoLG_[i]] += F_e[i]
    end
end
F[begin:end-1]
    ",
) <comparacao:base-iter>

A implementação de `Iter and Ref`
procura esconder o reaproveitamento de memória
dentro da abstração de iterador.
Então, em @comparacao:ref-refiter,
vemos que temos bem mais linhas de código implícitas.

#compare-julia(
    `Ref`, `Iter and Ref`,
    "






local F = fill(0.0, (m+1,))
local Fe = Ref(fill(0.0, (dim,)))
local x2xis = Ref(fill(0.0, (size(phis)[1:2]..., sdim)))
local dx2xis = Ref(fill(0.0, (size(phi_derivs)[1], size(phi_derivs)[1], dim)))
for e in 1:N_e
    local LGe = view(LG, :, e)
    local Xe = view(X, LGe)
    local Ye = view(Y, LGe)

    x2xis_f_ref!(x2xis, phis, Xe, Ye)
    dx2xis_f_ref!(dx2xis, phi_derivs, Xe, Ye)

    build_small_vec_2d_ref!(
        Fe,
        f,
        x2xis[], dx2xis[],
        phis, ws, gauss_n
    )

    local EQoLG_ = view(EQoLG, :, e)
    for i in 1:dim
        F[EQoLG_[i]] += Fe[][i]
    end
end
F[begin:end-1]
    ", "
local iter = FiniteElementIterRef(
    X, Y,
    N_e, LG,
    phis, phi_derivs,
)

local F = fill(0.0, (m+1,))





for (e, ref_x2xis, ref_dx2xis, ref_Ke, ref_Fe) in iter
    local x2xis = ref_x2xis[]
    local dx2xis = ref_dx2xis[]




    build_small_vec_2d_ref!(
        ref_Fe,
        f,
        x2xis, dx2xis,
        phis, ws, gauss_n
    )

    local EQoLG_ = view(EQoLG, :, e)
    for i in 1:dim
        F[EQoLG_[i]] += ref_Fe[][i]
    end
end
F[begin:end-1]
    ",
) <comparacao:ref-refiter>
#fakepar()

== Artefatos Entregues

#let commit = "06321afce88290a35a9b5d822c021d58d88fefb4"
#let github-repo = "https://github.com/Kiyoshi364/finite-elements/"
#let github-commit(commit: commit) = github-repo + "tree/" + commit + "/"
#let src-folder = "src/09-finite-elements-any-quadrilateral/"
#let link-raw(the-link, body) = {
    link(the-link, box(raw(body)))
}
#let link-file(file-name, commit: commit, ..args) = {
    let args_pos = args.pos()
    let body = if args_pos.len() == 1 {
        args_pos.at(0)
    } else if args_pos.len() > 1 {
        panic("More than one extra pos argument to link-file:", ..args_pos)
    } else {
        file-name
    }
    link-raw(github-commit(commit: commit) + src-folder + file-name, body)
}

A implementação desenvolvida durante esse trabalho
está disponível em
#link-raw(github-repo, github-repo).
O commit usado para produzir os resultados é
#link-raw(github-commit(), commit).
O código está na pasta
#link-file(src-folder).
O arquivo
#link-file("finite-elements.jl")
tem todas as implementações.
O script
#link-file("tests.jl")
testa a corretude
dos bloquinhos usados pelas implementações.
O arquivo
#link-file("bacarmo.jl")
é uma implementação adaptada
de Bruno Carmo
#footnote[implementação original disponível em
#{
let l = "https://github.com/bacarmo/Elementos_Finitos/blob/dba5361b423a32933118b3fca65edb93fe1cc425/2024_02/Estacionario_2D_equacao1_Pluto.jl"
link-raw(l, l)
}].
O intuito de incluir essa implementação
nos benchmarks é
ter uma comparação com outras implementações.
O script
#link-file("make_benchmarks.jl")
faz as benchmarks
e gera o arquivo `results.json`.
O script
#link-file("compile_benchmarks.jl")
analiza o `results.json` e gera
tabelas e gráficos mostrados na @sec:Resultados.
Os resultados do benchmark
que usamos para gerar as tabelas e gráficos
da @sec:Resultados
estão em
#link-file("results-" + commit + ".json").
Finalmente,
o script
#link-file("errors.jl")
gera o gráfico de convergência de erro
da implementação padrão (`Ref`).

Em uma versão futura de
#{
let future-commit = "beae74b93459cc28e7f57e05e2400d773a463ba4"
link-file("errors.jl")
[
    (no commit
    #link-raw(github-commit(commit: future-commit), future-commit)),
]
}
o script gera
os gráficos de convergência do erro
para todas as implementações.
Na @conv-erros,
mostramos que todas as implementações funcionam
fazendo gráficos de convergência de erro
para cada implementação.

#let error-grid(names, files) = figure(
    grid(
        columns: (1fr, 1fr, 1fr, 1fr),
        column-gutter: (1em,),
        align: center,
        ..array.range(files.len()).map(i => {
            let file = files.at(i)
            let name = names.at(i)
            image(images-folder + "error-" + file + ".svg")
            name
        }),
    ),
    caption: [Convergência de Erros],
    placement: auto,
    supplement: [Convergência de Erros],
    kind: "errors-onvergence",
)

#let impls_names = (`Baseline`, `Ref`, `Iter`, `Iter and Ref`, `Bruno Carmo`)
#let impls_files = ("baseline", "ref", "iter", "iter_ref", "bacarmo")

#error-grid(impls_names.slice(0, impls_names.len() - 1), impls_files.slice(0, impls_files.len() - 1))
<conv-erros>

= Resultados <sec:Resultados>

Nessa seção,
mostramos e discutimos os resultados
dos benchmarks realizados sobre as implementações.

== Resultados locais

#let readcsv = name => csv(images-folder + name + ".csv")
#let csv_as_table(name, it) = figure(
    grid(
        columns: (auto,) + ((1fr,) * (it.at(0).len() - 1)),
        row-gutter: 0.5em,
        align: (left,) + ((center,) * (it.at(0).len() - 1)),
        ..it.at(0).map(it => align(center, [*#it*])),
        ..it.slice(1).map(row => (
            raw(row.at(0)),
            ..row.slice(1).map(s => {
                let ss = s.split(".")
                if ss.len() == 1 {
                    ss.at(0)
                } else {
                    let len = ss.at(1).len()
                    let acc = 3
                    let min = calc.min(acc,len)
                    ss.at(0)
                    [.]
                    ss.at(1).slice(0,min)
                    [0] * (acc - min)
                }
            })
        )).flatten()
    ),
    caption: [#name],
    placement: auto,
    supplement: [Benchmarks],
    kind: "benchmarks",
)

#csv_as_table([Vetores Locais ($FF^e$)], readcsv("smallvec")) <bench-smallvec>
#csv_as_table([Matrizes Locais ($KK^e$)], readcsv("smallmat")) <bench-smallmat>

Nos @bench-smallvec e @bench-smallmat,
mostramos os tabelas com os resultados de benchmark
para a construção de $FF^e$ e $KK^e$,
respectivamente.
Em ambas as tabelas,
aparecem uma versão `with precalculation`.
A diferença entre a versão normal e a `with precalculation`
é que
a versão normal inclui o tempo do pré-calculo
de valores auxiliares
(que seriam compartilhados
na construção simultânea de $FF$ e $KK$)
e valores auxilares para a conversão de coordenada
(apenas no benchmark da $FF^e$).

Em ambos os resultados,
percebemos que,
tanto para construção de $FF^e$ quanto $KK^e$,
as melhores implementações foram
as `Ref`, seguidas pelas `Baseline`
e então as do `Bruno Carmo`.
Também que as versões `Ref`
apenas realizam alocações durante a fase de pré-calculo.

== Resultados Globais

#let magic_grid(name, names, files) = {
    let types = ("vec", "mat", "both")
    let len = names.len()
    let g = grid(
        columns: (auto, 1fr, 1fr, 1fr),
        column-gutter: (1em,),
        row-gutter: (0.5em,),
        align: horizon,
        ..([], [Vetor ($FF$)], [Matriz ($KK$)], [Ambos ($FF$ e $KK$)]).map(body => align(center, body)),
        ..(array.range(len)).map(i => (
            names.at(i),
            types.map(type =>
                image(images-folder + type + "-" + files.at(i) + ".svg")
            ),
        )).flatten()
    )
    figure(
        g,
        caption: [#name],
        placement: auto,
        supplement: [Benchmarks],
        kind: "benchmarks",
    )
}

#let graph_names = ("Tempo", "Memória", "Alocações")
#let graph_files = ("min_time", "min_memory", "min_alloc")

Nos @bench-compareimpl e @bench-eachimpl,
mostramos os resultados dos benchmarks
para as implementações globais.
@bench-compareimpl
mostra todas as implementações
de forma comparativa entre si,
enquanto
@bench-eachimpl
mostra cada implementação isoladamente.
Todos os gráficos,
tem em sua coordenada x
os números de elementos finitos da matriz construída.
Todos os eixos estão em escala logarítmica.

Nas colunas de construção de $FF$ percebemos
um ganho em todas as implementações
comparados com `Bruno Carmo`.
Entretanto esse ganho
desaparece na construção de $KK$.
Na construção de $KK$,
todas as implementações que não são `Bruno Carmo`
aparentam usar tempo assintoticamente maior que exponencial;
e por isso,
acaba dominando o tempo
da construção de ambas matrizes $FF$ e $KK$.

Sobre as implementações,
o reuso de memória em `Ref`
teve uma melhora nas três métricas
tempo, memória e alocações
comparando com a sua versão sem o reuso `Baseline`;
enquanto o reuso de memória em `Iter and Ref`
apenas reduziu o uso de memória na construção de $FF$
comparando com `Iter`.
Já a utilização de um iterador
como uma conveniência não contribuiu
para melhoras em execução:
`Ref` sempre foi melhor que `Iter and Ref`
e `Baseline` sempre foi melhor que `Iter`.
`Ref` aparenta ser a melhor implementação.

#magic_grid([Comparação de Implementações], graph_names, graph_files)
<bench-compareimpl>

#magic_grid([Resultados para cada Implementação], impls_names, impls_files)
<bench-eachimpl>

= Conclusão

Nesse trabalho implementamos
construção de matrizes
para o método de elementos finitos
na linguagem de programação Julia
de quatro formas diferentes
buscando economizar recursos.
Então, comparamos as implementações entre si
e com uma implementação externa.
Todas as quatro implementações buscaram
pre-calcular valores,
com destaque em $varphi$ e $nabla varphi$
aplicados nos pontos de gauss.

As principal diferença entre as implementações
foi o uso ou não de reaproveitamento de memória
e uso ou não de um iterador para conveniência da implementação.
Parece que
as quatro implementações tiveram algum problema
na montagem da matriz $KK$:
o tempo utilizado para construir a matriz
cresce assintoticamente mais rápido que o esperado.

Com isso, concluímos que,
mesmo com o defeito de implementação,
o pre-cálculo desses valores
contribui para a redução do tempo.
Também que reutilizar a memória
contribui tanto para utilizar e alocar menos memória
quanto para construír a matriz em menos tempo.
E finalmente que
essa implementação de iterador,
focada em dar conveniência ao programador,
gera gastos de tempo, memória e alocação.
Acreditamos que esse gasto extra,
da-se por o compilador não ser capaz
de desfazer dessa abstração,
então impedindo otimizações que seria capaz
nas outras implementações.

== Trabalhos Futuros

Entre os trabahos futuros temos três pontos:
entender porque a construção de $KK$ está lenta;
usar a memória de $KK^e$ para construir $FF^e$; e
reaproveitar memória na construção parcial da LG.
Sobre o primeiro ponto,
temos que o código de construção de $KK$
é análogo ao código de construção de $FF$ e
à implementação de Bruno Carmo.
Isso se faz surpreendente que a implementação
demore consideravelmente mais.
Talvez isso tenha acontecido por algum motivo
de otimização ou "pessimização"
do compilador.

O segundo ponto é mais simples.
Na implementação `Ref` da construção das matrizes locais,
passamos uma região de memória para armazenar o resultado.
Atualmente reservamos
uma região para $KK^e$ e outra para $FF^e$,
mas podemos reutilizar a mesma região para as duas
se realizarmos as operações na seguinte ordem:
construir matriz local 1; mover para matriz global 1;
construir matriz local 2; e então mover para matriz global 2.
Chamo as matrizes $FF$ e $KK$ (locais e globais)
de matrizes 1 e 2,
pois a ordem não afeta o algorítmo.
Não imaginamos que essa otimização
vai trazer grandes contribuições,
pois estamos apenas removendo
1 alocação de $n times ("size of float")$,
onde $n$ é o número de elementos finitos
influenciados por uma função da base.
Em uma base linear e um sub-espaço de $RR^2$,
economizaríamos 32 bytes para 4 floats de 64 bits.

Finalmente,
o último ponto é sobre economizar memória na LG.
O mapa LG é utilizado para traduzir
índices locais para índices globais,
muito importante para a construção de $FF$ e $KK$.
Geralmente ela é implementada como uma matriz,
sendo uma forma de "tabela de lookup".
Em muitos casos de particionamento do espaço,
essa tabela é completamente previsível
dados poucos parâmetros.
Julia possui um grande sistema de interfaces
e é possível fazer uma implementação de uma `AbstractMatrix`:
um objeto que se comporta como uma matriz normal,
mas pode ter implementações mais espertas por baixo dos panos.
Um exemplo bem conhecido de objeto
que implementa essa interface é `SparseMatrixCSC`
que considera que muitos elementos da matriz
são iguais (geralmente a 0)
e apenas guarda os a posição e valor dos elementos diferentes.
