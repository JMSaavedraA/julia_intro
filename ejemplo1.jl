include("elementoFinito.jl")
include("metodosTridiagonales.jl")
include("Solvers.jl")
using Random

Random.seed!(2022);
m=2001
k=501
n=101
fuzzy=1.0
x=LinRange(-5,5,k)
y=sinpi.(x) + fuzzy*rand(k) .- fuzzy/2 #y_k es sin(pi*x_k) + u_k, donde u_k es aleatorio
z=LinRange(-5,5,m)
fz=sinpi.(z)
w=LinRange(-5,5,n)

#Primero, consideramos λ=0.5

lambda=0.5
fz1=elementoFinitoLineal(z,x,y,w,lambda)

lambda=1.5
fz2=elementoFinitoLineal(z,x,y,w,lambda)

lambda=3.5
fz3=elementoFinitoLineal(z,x,y,w,lambda)