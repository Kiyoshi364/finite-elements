#let template(
  // Text
  font-type: "New Computer Modern",
  lang: "pt",
  region: "br",

  // Par
  indent: 2em,

  // Me
  title: none,
  author: none,
  mainbody,
) = {
  set page(
    paper: "a4",
    margin: (left: 2.0cm, right: 2.0cm, top: 2.5cm, bottom: 2.5cm),
    numbering: "1/1",
  )

  set text(
    font: font-type,
    lang: lang,
    region: region,
  )

  set par(
    justify: true,
    first-line-indent: indent,
  )

  set heading(
    numbering: "1.",
  )

  show heading: c => {
    c
    par(text(fill: white, [a #v(-1.5em)]))
  }

  align(center)[
    #text(
      weight: "bold",
      size: 18pt,
      title,
    )

    #author
  ]
  v(1em)

  mainbody
}
