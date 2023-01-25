include("Solvers.jl")
include("metodosTridiagonales.jl")


function f(x,t)
    #Función del problema de Sitnikov
    y=-x/(x^2+t^2)^(1.5)
    return y
end

function fPrime(x,r)
    #Derivada de f
    y=(2*x^2-r^2)/(x^2+r^2)^(2.5)
    return y
end

function Jacobiana(x::AbstractVector,r::AbstractVector,n::Int)
    #Jacobiana de F
    J=zeros(Float64,n+1,5)
    J[3:n+1,2]=h*h*fPrime.(x[2:n],r[2:n])
    return J
end

function ladoDerecho(x::AbstractVector,r::AbstractVector,n::Int,z0::Number,z0dot::Number,h::AbstractFloat)
    #F(z) del problema
    d=zeros(Float64,n+1)
    d[3:n+1]=h*h*f.(x[2:n],r[2:n])
    d[1]=z0
    d[2]=h*z0dot
    return d
end

function sitnikov(z0::Number,z0dot::Number,r::AbstractVector,n::Int,tol::AbstractFloat,maxIter::Int,h::AbstractFloat)
    #Diferencias finitas para el problema de Sitnikov
    #Construcción de A
    A=zeros(Float64,n+1,5); A[:,1].=11/12; A[:,2].=-5/3; A[:,3].=1/2; A[:,4].=1/3; A[:,5].=-1/12;
    A[1,:].= 0; A[1,3]=1;
    A[2,1]=0; A[2,2]=-11/6; A[2,3]=3; A[2,4]=-3/2; A[2,5]=1/3;
    A[n,1]=1; A[n,2]=-2; A[n,3]=1; A[n,4]=0; A[n,5]=0;
    A[n+1,1]=1; A[n+1,2]=-2; A[n+1,3]=1; A[n+1,4]=0; A[n+1,5]=0;
    #inicialización de z
    z=zeros(Float64,n+1).+z0
    #T(z)=Az-F(z)
    condicion=multiplicaNdiagonalVector(A,z,n+1,2)-ladoDerecho(z,r,n,z0,z0dot,h)
    #inicialización de delta
    delta=zeros(Float64,n+1)
    #Norma de T(z)
    nonStop=norma2(condicion,n+1)
    i=0
    while nonStop>tol && i<maxIter
        J=Jacobiana(z,r,n)
        #LDU de la Jacobiana (Especial para J pentadiagonal)
        L=LDUpentadiagonal(A-J,n+1)
        #Resolver para delta
        resuelveLDUPentadiagonalDada(L,condicion,n+1,delta)
        #Actualización por Newton
        z=z-delta
        condicion=multiplicaNdiagonalVector(A,z,n+1,2)-ladoDerecho(z,r,n,z0,z0dot,h)
        nonStop=norma2(condicion,n+1)
        i+=1
    end
    return z
end

function eulerTheta(n::Int,theta::AbstractVector,h::AbstractFloat,c::Number,exc::Number,p::Number)
    #Solución por diferencias finitas de θ'=4c(1+ϵcos(θ))^2/p^2
    theta[2]=theta[1]+ h*c*(1+exc*cos(theta[1]))^2/(p^2)
    for i=1:n-1
        theta[i+2]=theta[i]+ 2*h*c*(1+exc*cos(theta[i+1]))^2/(p^2)
    end
end