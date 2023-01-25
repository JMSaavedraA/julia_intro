using Markdown, Weave

weave(joinpath("", "introduccion.jmd"),
      informat="markdown",
      out_path = :pwd,
      doctype = "md2html")
