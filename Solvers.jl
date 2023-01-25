using DelimitedFiles
using Printf


function resuelveDiagonal(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    x=zeros(Float64,m)
    for i in 1:m
        try
            x[i]=baux[i]/Aaux[i,i]
        catch
            print("Error, matriz singular")
        end
    end
    return x
end
function resuelveDiagonal(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    for i in 1:m
        try
            x[i]=baux[i]/Aaux[i,i]
        catch
            print("Error, matriz singular")
        end
    end
end

function resuelveDiagonal(Aaux::AbstractVector,baux::AbstractVector,m::Int)
    x=zeros(Float64,m)
    for i in 1:m
        try
            x[i]=baux[i]/Aaux[i]
        catch
            print("Error, matriz singular")
        end
    end
    return x
end

function resuelveDiagonal(Aaux::AbstractVector,baux::AbstractVector,m::Int,x::AbstractVector)
    for i in 1:m
        try
            x[i]=baux[i]/Aaux[i]
        catch
            print("Error, matriz singular")
        end
    end
end

function Identidad(m::Int)
    I=diagonal(zeros(Float64,m).+1.0,m)
    return I
end

function invierteDiagonal(A::AbstractVector,m::Int)
    Aaux=zeros(Float64,m,m)
    for i in 1:m
        try
            Aaux[i,i]=1.0/A[i]
        catch
            print("Error, matriz singular")
        end
    end
    return Aaux
end

function invierteDiagonal(A::AbstractMatrix,m::Int)
    Aaux=zeros(Float64,m,m)
    for i in 1:m
        try
            Aaux[i,i]=1.0/A[i,i]
        catch
            print("Error, matriz singular")
        end
    end
    return Aaux
end

function invierteDiagonalVector(A::AbstractVector,m::Int)
    Aaux=zeros(Float64,m)
    for i in 1:m
        try
            Aaux[i]=1.0/A[i]
        catch
            print("Error, matriz singular")
        end
    end
    return Aaux
end

function invierteDiagonalVector(A::AbstractMatrix,m::Int)
    Aaux=zeros(Float64,m)
    for i in 1:m
        try
            Aaux[i]=1.0/A[i,i]
        catch
            print("Error, matriz singular")
        end
    end
    return Aaux
end

function determinanteTriangular(Aaux::AbstractMatrix,m::Int)
    l=1.0
    for i in 1:m
        l*=Aaux[i,i]
    end
    return l
end

function determinanteTriangular(Aaux::AbstractVector,m::Int)
    l=1.0
    for i in 1:m
        l*=Aaux[i]
    end
    return l
end

function resuelveTriangularInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    x=zeros(Float64,m)
    for i in 1:m
        x[i]=(baux[i]-sum(x[1:i-1].*Aaux[i,1:i-1]))/Aaux[i,i]
    end
    return x
end


function resuelveTriangularSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    x=zeros(Float64,m)
    for i in 1:m
        j=m-i+1
        x[j]=(baux[j]-sum(x[j+1:end].*Aaux[j,j+1:end]))/Aaux[j,j]
    end
    return x
end

function resuelveTriangularInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    for i in 1:m
        x[i]=(baux[i]-sum(x[1:i-1].*Aaux[i,1:i-1]))/Aaux[i,i]
    end
end

function resuelveTriangularSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    for i in 1:m
        j=m-i+1
        x[j]=(baux[j]-sum(x[j+1:end].*Aaux[j,j+1:end]))/Aaux[j,j]
    end
end

function resuelveGaussJordan(A::AbstractMatrix,b::AbstractVector,m::Int)
    Aaux=copy(A)
    baux=copy(b)
    orden=collect(1:m)
    for i in 1:m
        K=argmax(abs.(Aaux[i:end,i:end]))
        k=K[1]+i-1
        l=K[2]+i-1
        if k!=i
            aux=Aaux[k,:]
            Aaux[k,:]=Aaux[i,:]
            Aaux[i,:]=aux
            aux=baux[k]
            baux[k]=baux[i]
            baux[i]=aux
        end
        if l!=i
            aux=Aaux[:,l]
            Aaux[:,l]=Aaux[:,i]
            Aaux[:,i]=aux
            aux=orden[l]
            orden[l]=orden[i]
            orden[i]=aux
        end
        aux=Aaux[i,:]/Aaux[i,i]
        baux[i]=baux[i]/Aaux[i,i]
        Aaux[i,:]=aux
        for j=i+1:m
            baux[j]-=baux[i]*Aaux[j,i]
            Aaux[j,:]=Aaux[j,:]-Aaux[j,i]*aux
        end
    end
    z=resuelveTriangularSuperior(Aaux,baux,m)
    x=zeros(Float64,m)
    for i in 1:m
        x[orden[i]]=z[i]
    end
    return x
end


function resuelveGaussJordan(A::AbstractMatrix,b::AbstractVector,m::Int,x::AbstractVector)
    Aaux=copy(A)
    baux=copy(b)
    orden=collect(1:m)
    for i in 1:m
        K=argmax(abs.(Aaux[i:end,i:end]))
        k=K[1]+i-1
        l=K[2]+i-1
        if k!=i
            aux=Aaux[k,:]
            Aaux[k,:]=Aaux[i,:]
            Aaux[i,:]=aux
            aux=baux[k]
            baux[k]=baux[i]
            baux[i]=aux
        end
        if l!=i
            aux=Aaux[:,l]
            Aaux[:,l]=Aaux[:,i]
            Aaux[:,i]=aux
            aux=orden[l]
            orden[l]=orden[i]
            orden[i]=aux
        end
        aux=Aaux[i,:]/Aaux[i,i]
        baux[i]=baux[i]/Aaux[i,i]
        Aaux[i,:]=aux
        for j=i+1:m
            baux[j]-=baux[i]*Aaux[j,i]
            Aaux[j,:]=Aaux[j,:]-Aaux[j,i]*aux
        end
    end
    z=resuelveTriangularSuperior(Aaux,baux,m)
    for i in 1:m
        x[orden[i]]=z[i]
    end
end

function croutLU(A::AbstractMatrix,m::Int)
    L=zeros(Float64,m,m)
    for i=1:m
        for j=1:i-1
            L[i,j]=A[i,j]-sum(L[i,1:j-1].*L[1:j-1,j])
        end
        L[i, i] = A[i, i] - sum(L[i,1:i-1].*L[1:i-1, i])
        L[i,i+1:end]=(transpose(A[i,i+1:end])-transpose(L[i,1:i-1])*L[1:i-1,i+1:end])./L[i,i]
    end
    return L
end

function doolittleLU(A::AbstractMatrix,m::Int)
    L=zeros(Float64,m,m)
    for i=1:m
        for j=1:i-1
            L[i,j]=(A[i,j]-sum(L[i,1:j-1].*L[1:j-1,j]))/L[j,j]
        end
        L[i, i] = A[i, i] - sum(L[i,1:i-1].*L[1:i-1, i])
        L[i,i+1:end]=transpose(A[i,i+1:end])-transpose(L[i,1:i-1])*L[1:i-1,i+1:end]
    end
    return L
end

function factorizarLDU(A::AbstractMatrix,m::Int)
    L=zeros(Float64,m,m)
    for i=1:m
        for j=1:i-1
            suma=0;
            for k=1:j-1
                suma+=L[i,k]*L[k,k]*L[k,j];
            end
            L[i,j]=(A[i,j]-suma)./L[j,j];
        end
        suma=0;
        for k=1:i-1
            suma+=L[i,k]*L[k,k]*L[k,i];
        end
        L[i, i] = A[i, i] - suma;
        for j=i+1:m
            suma=0;
            for k=1:j-1
                suma+=L[i,k]*L[k,k]*L[k,j];
            end
            L[i,j]=(A[i,j]-suma)/L[i,i];
        end
    end
    return L
end

function lower(X::AbstractMatrix,m::Int)
    Xaux=zeros(Float64,m,m)
    for i=2:m
        Xaux[i,1:i-1]=X[i,1:i-1]
    end
    return Xaux
end

function upper(X::AbstractMatrix,m::Int)
    Xaux=zeros(Float64,m,m)
    for i=1:m-1
        Xaux[i,i+1:end]=X[i,i+1:end]
    end
    return Xaux
end

function diagonal(X::AbstractMatrix,m::Int)
    Xaux=zeros(Float64,m,m)
    for i=1:m
        Xaux[i,i]=X[i,i]
    end
    return Xaux
end

function diagonal(X::AbstractVector,m::Int)
    Xaux=zeros(Float64, m, m)
    for i=1:m
        Xaux[i,i]=X[i]
    end
    return Xaux
end

function vectorDiagonal(X::AbstractMatrix,m::Int)
    Xaux=zeros(Float64, m)
    for i=1:m
        Xaux[i]=X[i,i]
    end
    return Xaux
end

function resuelveCrout(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    L=croutLU(Aaux,m)
    z=resuelveTriangularInferior(L,baux,m)
    for i=1:m
        L[i,i]=1
    end
    x=resuelveTriangularSuperior(L,z,m)
    return x
end

function resuelveDoolittle(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    L=doolittleLU(Aaux,m)
    d=vectorDiagonal(L,m)
    for i=1:m
        L[i,i]=1.0
    end
    z=resuelveTriangularInferior(L,baux,m)
    for i=1:m
        L[i,i]=d[i]
    end
    x=resuelveTriangularSuperior(L,z,m)
    return x
end

function resuelveLDU(A::AbstractMatrix,b::AbstractVector,m::Int)
    L=factorizarLDU(A,m)
    d=vectorDiagonal(L,m)
    for i=1:m
        L[i,i]=1.0
    end
    z=resuelveTriangularInferior(L,b,m)
    y=resuelveDiagonal(d,z,m)
    x=resuelveTriangularSuperior(L,y,m)
    return x
end

function resuelveCrout(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    L=croutLU(Aaux,m)
    z=resuelveTriangularInferior(L,baux,m)
    for i=1:m
        L[i,i]=1
    end
    resuelveTriangularSuperior(L,z,m,x)
end

function resuelveDoolittle(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    L=doolittleLU(Aaux,m)
    d=vectorDiagonal(L,m)
    for i=1:m
        L[i,i]=1.0
    end
    z=resuelveTriangularInferior(L,baux,m)
    for i=1:m
        L[i,i]=d[i]
    end
    resuelveTriangularSuperior(L,z,m,x)
end

function resuelveLDU(A::AbstractMatrix,b::AbstractVector,m::Int,x::AbstractVector)
    L=factorizarLDU(A,m)
    d=vectorDiagonal(L,m)
    for i=1:m
        L[i,i]=1.0
    end
    z=resuelveTriangularInferior(L,b,m)
    y=resuelveDiagonal(d,z,m)
    resuelveTriangularSuperior(L,y,m,x)
    return x
end

function guardarVector(x0::AbstractVector,m::Int,nombre::String)
    open(nombre, "w") do f # "w" for writing JULIA
    write(f, @sprintf("%i 1\n",m)) # \n for newline
    for i=1:m-1
        write(f, @sprintf("%f\n",x0[i]))
    end
    write(f, @sprintf("%f",x0[m]))
    end
end

function guardarMatriz(A::AbstractMatrix,m::Int,n::Int,nombre::String)
    open(nombre, "w") do f # "w" for writing JULIA
    write(f, @sprintf("%i %i\n",m,n)) # \n for newline
    for i=1:m-1
        for j=1:n-1
            write(f, @sprintf("%.20f ",A[i,j]))
        end
        write(f, @sprintf("%.20f\n",A[i,n]))
    end
    for j=1:n-1
        write(f, @sprintf("%.20f ",A[n,j]))
    end
    write(f, @sprintf("%.20f",A[n,n]))
    end
end

function leerMatriz(nombre::String)
    A = convert(Array{Float64,2},readdlm(nombre,' ')[2:end,1:end]);
    M = convert(Array{Int64,1},readdlm(nombre,' ')[1,1:2]);
    m=M[1]
    n=M[2]
    return A, m, n
end

function leerVector(nombre::String)
    b = convert(Array{Float64,1},readdlm(nombre,'\n')[2:end]);
    return b
end

function norma2(x::AbstractVector,m::Int)
    #La norma del vector x de tamaño m
    suma=0
    for i=1:m
        suma+=x[i]*x[i]
    end
    n=sqrt(suma)
    return n
end

function CholeskyLDL(Aaux::AbstractMatrix,m::Int)
    #Esta es la factorizacion de Cholesky A=LDL^t para A simetrica de tamaño mxm
    L=zeros(Float64,m,m)
    for i=1:m
        sumaExterna=0.0
        for k=1:i-1
            vAux=0.0
            suma=0.0
            for j=1:k-1
                suma+=L[i,j]*L[j,j]*L[k,j]
            end
            vAux=(Aaux[i,k]-suma)
            L[i,k]=vAux/L[k,k]
            sumaExterna+=vAux*vAux/L[k,k]
        end
        L[i,i]=Aaux[i,i]-sumaExterna
    end
    return L
end

function metodoJacobi(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Jacobi para un problema en general Ax=y, A de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld.=xNew
        xNew=similar(xOld)
        for j=1:m
            suma=0
            for k=1:j-1
                suma+=A[j,k]*xOld[k]
            end
            for k=j+1:m
                suma+=A[j,k]*xOld[k]
            end
            xNew[j]=(y[j] - suma)/A[j,j]
        end
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i #Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoGaussSeidel(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Gauss-Seidel para un problema en general Ax=y, A de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld.=xNew
        for j=1:m
            suma=0
            for k=1:j-1
                suma+=A[j,k]*xNew[k]
            end
            for k=j+1:m
                suma+=A[j,k]*xNew[k]
            end
            xNew[j]=(y[j] - suma)/A[j,j]
        end
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i #Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoJacobiThreads(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Jacobi para un problema en general Ax=y, A de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld.=xNew
        xNew=similar(xOld)
        Threads.@threads for j=1:m
            suma=0
            for k=1:j-1
                suma+=A[j,k]*xOld[k]
            end
            for k=j+1:m
                suma+=A[j,k]*xOld[k]
            end
            xNew[j]=(y[j] - suma)/A[j,j]
        end
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i #Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoGaussSeidelThreads(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Gauss-Seidel para un problema en general Ax=y, A de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld.=xNew
         Threads.@threads for j=1:m
            suma=0
            for k=1:j-1
                suma+=A[j,k]*xNew[k]
            end
            for k=j+1:m
                suma+=A[j,k]*xNew[k]
            end
            xNew[j]=(y[j] - suma)/A[j,j]
        end
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i #Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end
using LinearAlgebra



function multiplicaMatricesThreads(A::AbstractMatrix,B::AbstractMatrix,C::AbstractMatrix,m::Int,n::Int,p::Int)
    #A=mxn,B=nxp,C=mxp
    Threads.@threads for i=1:m
        for j=1:p
            Cij = zero(eltype(C))
            for k=1:n
                Cij += A[i,k] * B[k,j]
            end
            C[i,j] = Cij
        end
    end
end

function multiplicaMatrizVectorThreads(A::AbstractMatrix,x::AbstractVector,c::AbstractVector,m::Int,n::Int)
    #A=mxn
    Threads.@threads for i=1:m
        cij = zero(eltype(c))
        for j=1:n
            cij += A[i,j] * x[j]
        end
        c[i,j] = cij
    end
end

function maxAbsolutoFueraDiagonal(A::AbstractMatrix,m::Int)
    v=0
    i=0
    j=0
    for k=1:m
        for l=1:k-1
            if v<abs(A[k,l])
                v=abs(A[k,l])
                i=k
                j=l
            end
        end
        for l=k+1:m
            if v<abs(A[k,l])
                v=abs(A[k,l])
                i=k
                j=l
            end
        end
    end
    return v,i,j
end

function multiplicaMatrices(A::AbstractMatrix,B::AbstractMatrix,C::AbstractMatrix,m::Int,n::Int,p::Int)
    #A=mxn,B=nxp,C=mxp
    for i=1:m
        for j=1:p
            Cij = zero(eltype(C))
            for k=1:n
                Cij += A[i,k] * B[k,j]
            end
            C[i,j] = Cij
        end
    end
end
function multiplicaMatrizVector(A::AbstractMatrix,x::AbstractVector,c::AbstractVector,m::Int,n::Int)
    #A=mxn
    for i=1:m
        cij = zero(eltype(c))
        for j=1:n
            cij += A[i,j] * x[j]
        end
        c[i,j] = cij
    end
end

function factorizaQR(A::AbstractMatrix,m::Int)
    Q=zeros(Float64,m,m)
    R=zeros(Float64,m,m)
    for i=1:m
        for j=1:i-1
            aAux=0.0
            for l=1:m
                aAux+=Q[l,j]*A[l,i]
            end
            R[j,i]=aAux
        end
        nAux=0.0
        for l=1:m
            aAux=A[l,i]
            for j=1:i-1
                aAux-=R[j,i]*Q[l,j]
            end
            Q[l,i]=aAux
            nAux+=aAux*aAux
        end
        nAux=sqrt(nAux)
        R[i,i]=nAux
        for l=1:m
            Q[l,i]=Q[l,i]/nAux
        end
    end
    return Q,R
end

function resuelveQR(A::AbstractMatrix,b::AbstractVector,m::Int)
    Q,R=factorizaQR(A,m)
    y=zeros(Float64,m)
    for i=1:m
        for j=1:m
            y[i]+=Q[j,i]*b[j]
        end
    end
    x=resuelveTriangularSuperior(R,y,m)
    return x
end

function resuelveQR(A::AbstractMatrix,b::AbstractVector,m::Int,x::AbstractVector)
    Q,R=factorizaQR(A,m)
    y=zeros(Float64,m)
    for i=1:m
        for j=1:m
            y[i]+=Q[j,i]*b[j]
        end
    end
    resuelveTriangularSuperior(R,y,m,x)
end
