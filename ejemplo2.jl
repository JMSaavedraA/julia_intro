using Random
using Statistics
using Distributions
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