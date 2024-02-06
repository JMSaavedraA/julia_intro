# script para instalar Plots en Julia
using Pkg
println("Instalando paquetes necesarios para graficar")
Pkg.update()
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
@warn("Packages installed!")