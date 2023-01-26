include("metodosTridiagonales.jl")
function elementoFinitoLineal(z::AbstractVector,x::AbstractVector,y::AbstractVector,w::AbstractVector,lambda::Number)
    m=length(x) #(x,y) son las observaciones (m observaciones)
    n=length(w) #w es el vector con los puntos de la malla
    r=length(z) #z es el punto donde se evalúa el Elemento Finito.
    H=zeros(Float64,n,3)
    b=zeros(Float64,n)
    J=zeros(Int64,r).+(n-1) #En este vector se guardará la entrada i tal que z_k ∈ [x_i,x_{i+1})
    Mi=zeros(Float64,r)
    Md=copy(Mi)
    z=sort(z) #Ordenamos el vector z
    s=sortperm(x) #Ordenamos (x,y) con respecto a x
    x=x[s]
    y.=y[s]
    w=sort(w) #Ordenamos el vector z
    i=1
    j=1
    zj=z[j]
    xi=x[i]
    for k=1:n-1
        wd=w[k+1] #calculamos extremos de intervalos
        wi=w[k]
        hk=wd+wi
        lk=wd-wi
        while zj<=wd && j<=r
            J[j]=k
            ji=(2*zj-hk)/lk
            Mi[j]=(1-ji)/2 #Calculamos las N_i de z_k
            Md[j]=(1+ji)/2 
            j+=1
            if j<=r
                zj=z[j]
            end
        end
        H[k,2]=H[k,2]+lambda/lk #Armamos la matriz tridiagonal para Phi primero con λ/h_i [1 -1]'[1 -1]
        H[k,3]=-lambda/lk
        H[k+1,1]=-lambda/lk
        H[k+1,2]=lambda/lk
        while xi<=wd && i<=m
            ji=(2*xi-hk)/lk
            Ni=(1-ji)/2 #Calculamos las N_i de x_i
            Nd=(1+ji)/2
            H[k,2]=H[k,2]+(Ni*Ni) #Sumamos [N_i N_i+1]' [N_i N_i+1] a la matriz tridiagonal
            H[k,3]=H[k,3]+(Ni*Nd)
            H[k+1,1]=H[k+1,1]+(Ni*Nd)
            H[k+1,2]=H[k+1,2]+(Nd*Nd)
            b[k]=b[k]+y[i]*Ni #Sumamos y_i [N_i N_i+1]' a b
            b[k+1]=b[k+1]+y[i]*Nd
            i+=1
            if i<=m
                xi=x[i]
            end
        end
    end
    phi=resuelveCholeskyTridiagonal(H,b,n) #Resolvemos para encontrar Φ
    fz=zeros(Float64,r)
    for j=1:r
        i=J[j]
        fz[j]=Mi[j]*phi[i] + Md[j]*phi[i+1] #Evaluamos en z_k
    end
    return fz
end