using Random
using Statistics
using Distributions
using Plots
using LaTeXStrings
using Printf
using KernelDensity

Random.seed!(2022)
k=10
n=3000
N=10000
x=zeros(Float64,n+1)
x[1]=2.5
alpha=2
r=100
estimador=zeros(Float64,r)

for j=1:r
    for i=1:n
        lambda=sqrt(x[i])+log(x[i])/x[i] 
        u=rand(Exponential(lambda),k)
        m=mean(u)
        x[i+1]=x[i]+3*log(i)*(alpha-m)/(i*log(N))
    end
    estimador[j]=x[n+1]
end
media=mean(estimador)
desv=sqrt(var(estimador))
estimador=(estimador.-media)./desv
U=kde(estimador)

plt1=histogram(estimador,normalize = :probability,bins=:sturges,label=L"\hat{x}",size=(1600,1200),legend=:topright,palette=:Oranges_3); plot!(U.x,pdf.(Normal(mean(estimador), sqrt(var(estimador))),U.x),lw=5,label="Densidad Normal");plot!(U.x,U.density,lw=5,label="Estimación por Kernel Normal");
savefig(plt1,@sprintf("histograma2-%i.png",n)) #Guardamos la gráfica