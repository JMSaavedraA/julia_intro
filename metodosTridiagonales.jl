function formaTridiagonal(A::AbstractMatrix,m::Int)
    #Función que nos regresa la diagonal de la matriz A bandada de tamaño de banda n y tamaño de A mxm
    Aaux=zeros(Float64, m, 3)
    Aaux[1,2]=A[1,1]
    Aaux[1,3]=A[1,2]
    for i=2:m-1
        Aaux[i,1]=A[i,i-1]
        Aaux[i,2]=A[i,i]
        Aaux[i,3]=A[i,i+1]
    end
    Aaux[m,1]=A[m,m-1]
    Aaux[m,2]=A[m,m]
    return Aaux
end

function reconstruyeTridiagonal(A::AbstractMatrix,m::Int)
    Aaux=zeros(Float64, m, m)
    Aaux[1,1]=A[1,2]
    Aaux[1,2]=A[1,3]
    for i=2:m-1
        Aaux[i,i-1]=A[i,1]
        Aaux[i,i]=A[i,2]
        Aaux[i,i+1]=A[i,3]
    end
    Aaux[m,m-1]=A[m,1]
    Aaux[m,m]=A[m,2]
    return Aaux
end

function CholeskyLDLTridiagonal(Aaux::AbstractMatrix,m::Int)
    #Esta es la factorizacion de Cholesky A=LDL^t para A bandada simetrica de tamaño de banda n y A de tamaño mxm
    L=zeros(Float64,m,3) #Guardamos en una matriz mxn solo la parte inferior de L, igualmente así se espera que sea A
    L[1,2]=Aaux[1,2]
    for i=2:m
        vAux=Aaux[i,1];
        L[i,1]=vAux/L[i-1,2]
        L[i-1,3]=L[i,1]
        L[i,2]=Aaux[i,2]-(vAux*vAux/L[i-1,2])
    end
    return L
end


function CholeskyLDLPentadiagonal(Aaux::AbstractMatrix,m::Int)
    #Esta es la factorizacion de Cholesky A=LDL^t para A bandada simetrica de tamaño de banda n y A de tamaño mxm
    L=zeros(Float64,m,5) #Guardamos en una matriz mxn solo la parte inferior de L, igualmente así se espera que sea A
    d2=Aaux[1,3]
    L[1,3]=d2
    a2=Aaux[2,2];
    L[2,2]=a2/d2
    L[1,4]=L[2,2]
    L[2,3]=Aaux[2,3]-(a2^2/d2)
    for i=3:m
        d1=L[i-2,3]
        a1=Aaux[i,1];
        d2=L[i-1,3]
        a2=Aaux[i,2]-a1*L[i-1,2];
        L[i,1]=a1/d1
        L[i-2,5]=L[i,1]
        L[i,2]=a2/d2
        L[i-1,4]=L[i,2]
        L[i,3]=Aaux[i,3]-(a1^2/d1 + a2^2/d2)
    end
    return L
end


function CholeskyLDLNdiagonal(Aaux::AbstractMatrix,m::Int,n::Int)
    #Esta es la factorizacion de Cholesky A=LDL^t para A bandada simetrica de tamaño de banda n y A de tamaño mxm
    L=zeros(Float64,m,n) #Guardamos en una matriz mxn solo la parte inferior de L, igualmente así se espera que sea A
    for i=1:m
        iShift=minimum([n-i,0])#Calculamos que entrada de L estamos calculando
        sumaExterna=0.0
        for k=1-iShift:i-1
            kShift=minimum([n-k,0])
            vAux=0.0
            suma=0.0
            for j=1-iShift:k-1
                jShift=minimum([n-j,0])
                suma+=L[i,j+iShift]*L[j,j+jShift]*L[k,j+kShift]
            end
            vAux=(Aaux[i,k+iShift]-suma)
            L[i,k+iShift]=vAux/L[k,k+kShift]
            sumaExterna+=vAux*vAux/L[k,k+kShift]
        end
        L[i,i+iShift]=Aaux[i,i+iShift]-sumaExterna
    end
    return L
end


function resuelveTridiagonalInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    #Eliminación gaussiana para Aaux triangular inferior bandada de banda n y tamaño m
    x=zeros(Float64, m)
    x[1]=baux[1]/Aaux[1,2]
    for i=2:m
        x[i]=(baux[i]-Aaux[i,1]*x[i-1])/Aaux[i,2]
    end
    return x
end

function resuelveTridiagonalInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    #Eliminación gaussiana para Aaux triangular inferior bandada de banda n y tamaño m
    x[1]=baux[1]/Aaux[1,2]
    for i=2:m
        x[i]=(baux[i]-Aaux[i,1]*x[i-1])/Aaux[i,2]
    end
end

function resuelveTridiagonalSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    #Eliminación gaussiana para Aaux triangular superior bandada de banda n y tamaño m
    x=zeros(Float64, m)
    x[m]=baux[m]/Aaux[m,2]
    for i=m-1:-1:1
        x[i]=(baux[i]-Aaux[i,3]*x[i+1])/Aaux[i,2]
    end
    return x
end

function resuelveTridiagonalSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    #Eliminación gaussiana para Aaux triangular superior bandada de banda n y tamaño m
    x[m]=baux[m]/Aaux[m,2]
    for i=m-1:-1:1
        x[i]=(baux[i]-Aaux[i,3]*x[i+1])/Aaux[i,2]
    end
end

function resuelveDiagonalTridiagonal(A::AbstractMatrix,b::AbstractVector,m::Int)
    #Solucion de Dx=b, donde D es la diagonal de A, una matriz bandada de banda n y tamaño mxm
    x=similar(b)
    for i=1:m
        x[i]=b[i]/A[i,2]
    end
    return x
end

function resuelveDiagonalTridiagonal(A::AbstractMatrix,b::AbstractVector,m::Int,x::AbstractVector)
    #Solucion de Dx=b, donde D es la diagonal de A, una matriz bandada de banda n y tamaño mxm
    for i=1:m
        x[i]=b[i]/A[i,2]
    end
end

function metodoJacobiTridiagonal(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Jacobi para un problema en general Ax=y, A tridiagonal simétrica de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld=copy(xNew)
        xNew=similar(xOld)
        xNew[1]=(y[1]-A[1,3]*xOld[2])/A[1,2]
        for j=2:m-1
            xNew[j]=(y[j]-A[j,1]*xOld[j-1]-A[j,3]*xOld[j+1])/A[j,2]
        end
        xNew[m]=(y[m] - A[m,1]*xOld[m-1])/A[m,2]
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i#Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function metodoGaussSeidelTridiagonal(A::AbstractMatrix,y::AbstractVector,m::Int,maxIter::Int,tol::AbstractFloat,xInit::AbstractVector)
    #El método iterativo de Gauss-Seidel para un problema en general Ax=y, A tridiagonal simétrica de tamaño mxm, maxIter el máximo de iteraciones, tol la tolerancia del error relativo, xInit el punto inicial
    xNew=copy(xInit)
    error=10
    i=0
    Cercania=zeros(Float64,maxIter)
    while error>tol && i<maxIter
        xOld=copy(xNew)
        xNew[1]=(y[1]-A[1,3]*xNew[2])/A[1,2]
        for j=2:m-1
            xNew[j]=(y[j]-A[j,1]*xNew[j-1]-A[j,3]*xNew[j+1])/A[j,2]
        end
        xNew[m]=(y[m] - A[m,1]*xNew[m-1])/A[m,2]
        error=norma2(xNew-xOld,m)/norma2(xNew,m)
        i+=1
        Cercania[i]=error
    end
    Cercania=Cercania[1:i]
    return xNew, Cercania, i#Regresamos el vector x aproximado, la cadena de errores de aproximación y la iteración en que encuentra x
end

function croutLUtridiagonal(A::AbstractMatrix,m::Int)
    L=similar(A)
    L[1,2]=A[1,2]
    L[1,3]=A[1,3]/L[1,2]
    for i=1:m
        L[i,1]=A[i,1]
        L[i,2]=A[i,2] - L[i,1]/L[i,2]
        L[i,3]=A[i,3]/L[i,2]
    end
    L[m,1]=A[m,1]
    L[m,2]=A[m,2] - L[m,1]/L[m,2]
    return L
end

function formaNdiagonal(A::AbstractMatrix,m::Int,n::Int)
    #Función que nos regresa la diagonal de la matriz A bandada de tamaño de banda n y tamaño de A mxm
    r=2*n+1
    Aaux=zeros(Float64, m, r)
    for i=1:n
        Aaux[i,n+2-i:r].=A[i,1:n+i]
    end
    for i=n+1:m-n
        Aaux[i,:]=A[i,i-n:i+n]
    end
    for i=n-1:-1:0
        Aaux[m-i,1:n+1+i].=A[m-i,m-n-i:m]
    end
    return Aaux
end

function reconstruyeNdiagonal(A::AbstractMatrix,m::Int,n::Int)
    r=2*n+1
    Aaux=zeros(eltype(A), m, m)
    for i=1:n
        Aaux[i,1:n+i].=A[i,n+2-i:r]
    end
    for i=n+1:m-n
        Aaux[i,i-n:i+n].=A[i,:]
    end
    for i=n-1:-1:0
        Aaux[m-i,m-n-i:m].=A[m-i,1:n+1+i]
    end
    return Aaux
end

function multiplicaNdiagonalVector(A::AbstractMatrix,x::AbstractVector,c::AbstractVector,m::Int,n::Int)
    r=2*n+1
    for i=1:n
        c[i]=A[i,n+2-i:r]'*x[1:n+i]
    end
    for i=n+1:m-n
        c[i]=A[i,:]'* x[i-n:i+n]
    end
    for i=n-1:-1:0
        c[m-i]=A[m-i,1:n+1+i]'*x[m-n-i:m]
    end
end

function multiplicaNdiagonalVector(A::AbstractMatrix,x::AbstractVector,m::Int,n::Int)
    c=zeros(Float64,m)
    r=2*n+1
    for i=1:n
        c[i]=A[i,n+2-i:r]'*x[1:n+i]
    end
    for i=n+1:m-n
        c[i]=A[i,:]'* x[i-n:i+n]
    end
    for i=n-1:-1:0
        c[m-i]=A[m-i,1:n+1+i]'*x[m-n-i:m]
    end
    return c
end


function resuelvePentadiagonalInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    #Eliminación gaussiana para Aaux triangular inferior bandada de banda n y tamaño m
    x=zeros(Float64, m)
    x[1]=baux[1]/Aaux[1,3]
    x[2]=(baux[2] - Aaux[2,2]*x[1])/Aaux[2,3]
    for i=3:m
        x[i]=(baux[i]-Aaux[i,1]*x[i-2] - Aaux[i,2]*x[i-1])/Aaux[i,3]
    end
    return x
end

function resuelvePentadiagonalInferior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    #Eliminación gaussiana para Aaux triangular inferior bandada de banda n y tamaño m
    x[1]=baux[1]/Aaux[1,3]
    x[2]=(baux[2] - Aaux[2,2]*x[1])/Aaux[2,3]
    for i=3:m
        x[i]=(baux[i]-Aaux[i,1]*x[i-2] - Aaux[i,2]*x[i-1])/Aaux[i,3]
    end
    return x
end


function resuelvePentadiagonalSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    #Eliminación gaussiana para Aaux triangular superior bandada de banda n y tamaño m
    x=zeros(Float64, m)
    x[m]=baux[m]/Aaux[m,3]
    x[m-1]=(baux[m-1]-Aaux[m-1,4]*x[m])/Aaux[m-1,3]
    for i=m-2:-1:1
        x[i]=(baux[i]-Aaux[i,4]*x[i+1]-Aaux[i,5]*x[i+2])/Aaux[i,3]
    end
    return x
end

function resuelvePentadiagonalSuperior(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    #Eliminación gaussiana para Aaux triangular superior bandada de banda n y tamaño m
    x[m]=baux[m]/Aaux[m,3]
    x[m-1]=(baux[m-1]-Aaux[m-1,4]*x[m])/Aaux[m-1,3]
    for i=m-2:-1:1
        x[i]=(baux[i]-Aaux[i,4]*x[i+1]-Aaux[i,5]*x[i+2])/Aaux[i,3]
    end
end

function resuelveCholeskyPentadiagonal(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    L=CholeskyLDLPentadiagonal(Aaux,m)
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,3]
        L[i,3]=1.0
    end
    y=resuelvePentadiagonalInferior(L,baux,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    x=resuelvePentadiagonalSuperior(L,y,m)
    return x
end

function resuelveCholeskyPentadiagonalDada(L::AbstractMatrix,d::AbstractVector,baux::AbstractVector,m::Int)
    y=resuelvePentadiagonalInferior(L,baux,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    x=resuelvePentadiagonalSuperior(L,y,m)
    return x
end

function resuelveCholeskyPentadiagonalDada(L::AbstractMatrix,d::AbstractVector,baux::AbstractVector,m::Int,x::AbstractVector)
    y=resuelvePentadiagonalInferior(L,baux,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    resuelvePentadiagonalSuperior(L,y,m,x)
end

function resuelveLDUPentadiagonalDada(L::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,3]
        L[i,3]=1
    end
    y=resuelvePentadiagonalInferior(L,baux,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    resuelvePentadiagonalSuperior(L,y,m,x)
end

function LDUpentadiagonal(A::AbstractMatrix,m::Int)
    L=zeros(Float64,m,5)
    L[1,3]=A[1,3]
    L[1,4]=A[1,4]/A[1,3]
    L[1,5]=A[1,5]/A[1,3]
    L[2,2]=A[2,2]/A[1,3]
    L[2,3]=A[2,3]-L[2,2]*L[1,3]L[1,4]
    L[2,4]=(A[2,4]-L[2,2]*L[1,3]L[1,5])/L[2,3]
    L[2,5]=A[2,5]/L[2,3]
    for i=3:m-2
        L[i,1]=A[i,1]/L[i-2,3]
        L[i,2]=(A[i,2]-L[i,1]*L[i-2,3]*L[i-2,4])/L[i-1,3]
        L[i,3]=A[i,3]-L[i,1]*L[i-2,3]*L[i-2,5] - L[i,2]*L[i-1,3]*L[i-1,4]
        L[i,4]=(A[i,4]-L[i,2]*L[i-1,3]*L[i-1,5])/L[i,3]
        L[i,5]=A[i,5]/L[i,3]
    end
    L[m-1,1]=A[m-1,1]/L[m-3,3]
    L[m-1,2]=(A[m-1,2]-L[m-1,1]*L[m-3,3]*L[m-3,4])/L[m-2,3]
    L[m-1,3]=A[m-1,3]-L[m-1,1]*L[m-3,3]*L[m-3,5] - L[m-1,2]*L[m-2,3]*L[m-2,4]
    L[m-1,4]=(A[m-1,4]-L[m-1,2]*L[m-2,3]*L[m-2,5])/L[m-1,3]
    L[m,1]=A[m,1]/L[m-2,3]
    L[m,2]=(A[m,2]-L[m,1]*L[m-2,3]*L[m-2,4])/L[m-1,3]
    L[m,3]=A[m,3]-L[m,1]*L[m-2,3]*L[m-2,5] - L[m,2]*L[m-1,3]*L[m-1,4]
    return L
end


function resuelveCholeskyTridiagonal(Aaux::AbstractMatrix,baux::AbstractVector,m::Int)
    L=CholeskyLDLTridiagonal(Aaux,m)
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,2]
        L[i,2]=1.0
    end
    y=resuelveTridiagonalInferior(L,baux,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    x=resuelveTridiagonalSuperior(L,y,m)
    return x
end

function resuelveCholeskyTridiagonal(Aaux::AbstractMatrix,baux::AbstractVector,m::Int,x::AbstractVector)
    L=CholeskyLDLTridiagonal(Aaux,m)
    d=zeros(Float64,m)
    for i=1:m
        d[i]=L[i,2]
        L[i,2]=1.0
    end
    y=resuelveTridiagonalInferior(L,baux,m)
    for i=1:m
        y[i]=y[i]/d[i]
    end
    resuelveTridiagonalSuperior(L,y,m,x)
end