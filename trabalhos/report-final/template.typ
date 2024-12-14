#let fakepar() = par(text(fill: white, [~ #v(-1.75em)]))

#let template(
  // Text
  text-size: 11pt,
  font-type: "New Computer Modern",
  lang: "pt",
  region: "br",

  // Par
  indent-first: true,
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
    size: text-size,
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
    if indent-first { fakepar() }
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
