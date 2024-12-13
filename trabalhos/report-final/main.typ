#import "template.typ": template

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

#show: template.with(
    title: title,
    author: author,
)

= Introdução

- elementos finitos
- performance
- implementaçoes
- sumário

= Especificação do Problema

Nesse trabalho estamos aplicando elementos finitos
no seguinte problema:

#problema[
    Dados
    uma função $f : accent(Omega, macron) -> RR$
    e constantes $alpha > 0, beta >= 0 : RR$,
    encontre
    uma função $u : accent(Omega, macron) -> RR$
    tal que
    #math.equation(
        block: true,
        numbering: "(1)"
    )[$
        cases(
            - alpha Delta u(arrow(x)) + beta u(arrow(x))
                = f(arrow(x))
                \,&#h(3em) x in Omega,
            u(arrow(x)) = 0
                \,&#h(3em) x in Gamma
        )
    $] <problema-forte>
    para algum sub-espaço $Omega subset RR^2$,
    com $Gamma$ fronteira de $Omega$ e
    $accent(Omega, macron) = Omega union Gamma$.
]

A versão fraca do @problema-forte[Problema]
é dado pela seguinte equação:
#math.equation(
    block: true,
    numbering: "(1)"
)[$
    alpha integral_Omega
        nabla u(arrow(x)) dot nabla v(arrow(x))
        d Omega
    + beta integral_Omega
        u(arrow(x)) v(arrow(x))
        d Omega
    =  integral_Omega
        f(arrow(x)) v(arrow(x))
        d Omega
$]

= Resultados

= Conclusão

== Trabalhos Futuros

- Aproveitar memória na LG
- Entender porque a construção da matriz está lenta

#lorem(200)
