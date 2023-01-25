include("Solvers.jl")
include("metodosTridiagonales.jl")
include("sitnikov.jl")
using Plots
using LaTeXStrings
using Printf

#Parámetros
p=1
c=2
exc=0.20
n=1000
T=10
t=LinRange(0,T,n+1)
z0=1.25
z0dot=0.0
tol=1e-8
maxIter=1000
h=t[2]-t[1]
theta=zeros(Float64,n+1)
#Resolvemos θ'=4c(1+ϵcos(θ))^2/p^2
@time eulerTheta(n,theta,h,c,exc,p)
#Calculamos r=p/(2(1+ϵcos(θ)))
R=(p/2)./(1 .+ exc*cos.(theta))
#Resolvemos z"=-z/(z^2+r^2)^()
@time Z=sitnikov(z0,z0dot,R,n,tol,maxIter,h)

#Límites de los ejes para la animación
X=R.*cos.(theta)
Y=R.*sin.(theta)
xlim=(-1,1)
ylim=(-1,1)
zlim=extrema(Z)

#creamos la animación
anim = @animate for i=2:n
    inicio=maximum([1,i-10])
    @views x1, y1, z1 = X[inicio:i], Y[inicio:i], Z[inicio:i].*0
    @views x2, y2, z2 = -X[inicio:i], -Y[inicio:i], Z[inicio:i].*0
    @views x3, y3, z3 = X[i:i]*0, Y[i:i]*0, Z[i:i]
    line=range(0, 10, length = i-inicio+1)
    A = plot(x1, y1, z1; label=L"q_1(t)", lw=line, xlim, ylim, zlim, camera=(36, 30),size=(800,600));plot!(x2, y2, z2; label=L"q_2(t)", lw=line, xlim, ylim, zlim, camera=(36, 30)); plot!(x3, y3, z3; label=L"z(t)", lw=0.3, xlim, ylim, zlim, camera=(36, 30),markersize = 2, seriestype=:scatter,markerstrokewidth=1)
end
gif(anim, @sprintf("sitnikov%f.gif",z0),fps=30)