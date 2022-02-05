module GrayCoding

using LinearAlgebra
using Gadfly

# Write your package code here.
"""
Generate Encoding and Decoding matrices for Gray Codes of alphabet.
```julia-repl
julia> G,B,g,b=GrayMatrix(4, 2);
julia> G
    4×4 Matrix{Int64}:
    1  0  0  0
    1  1  0  0
    1  1  1  0
    1  1  1  1
    julia> B
    4×4 Matrix{Int64}:
    1  0  0  0
    1  1  0  0
    0  1  1  0
    0  0  1  1
    julia> g 
    4×16 Matrix{Int64}:
    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1
    0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1
    0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1
    0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1
    julia> b 
    4×16 Matrix{Int64}:
    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1
    0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0
    0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0
    0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0
```
"""
function GrayMatrix(n::Int64=3, q::Int64=2)
    M = q^n
    G = Array(Bidiagonal(ones(1, n)[:], (q - 1) * ones(1, n - 1)[:], :L))
    B = mod.(inv(G), q)
    G = convert.(Int64, G)
    B = convert.(Int64, B)
    x = 0:M-1
    a = digits.(x, base = q, pad = n)  # left-lsb
    # a=digits.(x,base=2,pad=num_bits)|> reverse  # left-msb
    b = hcat(a...)
    b = reverse(b, dims = 1)
    g = Int.(mod.(G * b, q))
    return B, G, b, g
end

"""
Plots a matrix into a 2D with labels. Optional arguments including colors
```julia-repl
julia> using Random;
julia> A= rand(0:9,10,10);
julia> matrixplot(A)
```
"""
function matrixplot(A;kwargs...)
    a,b=size(A)
    X = transpose(repeat(1:b, 1, a))[:]
    Y = repeat(a:-1:1, b)[:]
    Gadfly.plot(x = X, y = Y, Geom.rectbin(),  color = A,alpha=[0.5], Coord.cartesian(fixed = true), Theme(bar_spacing = 0.1cm), Geom.label(position = :centered), label = string.(A)[:], Theme(key_position = :none, grid_line_width = 0pt, panel_stroke = nothing), Guide.xticks(ticks = nothing, label = false), Guide.yticks(ticks = nothing, label = false), Guide.xlabel(nothing), Guide.ylabel(nothing);kwargs...)
end

"""
Plot the DNA codon matrix
```julia-repl
julia> dnamatrix()
```
"""
function dnamatrix()
	U = [(x,y) for x ∈ 0:7, y ∈ 0:7]
	B,G,b,g=GrayMatrix(3,2)
	V = [string(g[:,y+1]...)string(g[:,x+1]...) for x ∈ 0:7, y ∈ 0:7]
	revV = [string(g[:,y+1]...)string(g[:,x+1]...) for x ∈ 0:7, y ∈ 0:7]
	Vx = [string(g[:,y+1]...) for x ∈ 0:7, y ∈ 0:7]
	Vy = [string(g[:,x+1]...) for x ∈ 0:7, y ∈ 0:7]
	VM=[(Vx[i,j][1]Vy[i,j][1]Vx[i,j][2]Vy[i,j][2]Vx[i,j][3]Vy[i,j][3]) for i ∈ 1:8, j ∈ 1:8]
	DNA=parse.(Int,VM,base=2)
	# DM0=[replace(Vx[i,j][1]Vy[i,j][1],"00"=>"C","01"=>"A","10"=>"U","11"=>"G")replace(Vx[i,j][2]Vy[i,j][2],"00"=>"C","01"=>"A","10"=>"U","11"=>"G")replace(Vx[i,j][3]Vy[i,j][3],"00"=>"C","01"=>"A","10"=>"U","11"=>"G") for i ∈ 1:8, j ∈ 1:8]
	DM=[replace(Vx[i,j][1]Vy[i,j][1],"00"=>"C","01"=>"A","10"=>"U","11"=>"G")replace(Vx[i,j][2]Vy[i,j][2],"00"=>"C","01"=>"A","10"=>"U","11"=>"G")replace(Vx[i,j][3]Vy[i,j][3],"00"=>"C","01"=>"A","10"=>"U","11"=>"G") for j ∈ 1:8, i ∈ 1:8]	
    AA=copy(DM) # Amino Acid
	replace!(AA,"CUU"=>"Leucine","CUC"=>"Leucine","CUA"=>"Leucine","CUG"=>"Leucine","UUA"=>"Leucine","UUG"=>"Leucine","UUU"=>"Phenylalanine","UUC"=>"Phenylalanine","AUC"=>"Isoleucine","AUA"=>"Isoleucine","AUU"=>"Isoleucine","AUG"=>"Methionine","GUA"=>"Valine","GUC"=>"Valine","GUU"=>"Valine","GUG"=>"START","UCA"=>"Serine","UCC"=>"Serine","UCU"=>"Serine","UCG"=>"Serine","CCC"=>"Proline","CCA"=>"Proline","CCU"=>"Proline","CCG"=>"Proline","ACU"=>"Threonine","ACA"=>"Threonine","ACC"=>"Threonine","ACG"=>"Threonine","GCC"=>"Alanine","GCU"=>"Alanine","GCA"=>"Alanine","GCG"=>"Alanine","GGU"=>"Glycine","GGA"=>"Glycine","GGC"=>"Glycine","GGG"=>"Glycine","CGA"=>"Arginine","CGC"=>"Arginine","CGU"=>"Arginine","CGG"=>"Arginine","GAU"=>"Aspartic acid","GAC"=>"Aspartic acid","GAA"=>"Glutamic acid","GAG"=>"Glutamic acid","AAU"=>"Asparagine","AAC"=>"Asparagine","UGU"=>"Cysteine","UGC"=>"Cysteine","UGG"=>"Tryptophan","CAA"=>"Glutamine","CAG"=>"Glutamine","UAA"=>"STOP","UAG"=>"STOP","UAU"=>"Tyrosine","UAC"=>"Tyrosine","AAA"=>"Lysine","AAG"=>"Lysine","CAC"=>"Histidine","CAU"=>"Histidine","AGG"=>"Arginine","AGA"=>"Arginine","AGU"=>"Serine","AGC"=>"Serine","UGA"=>"STOP" )
	return DM,VM,AA,Vx,Vy,DNA
end

"""
Decimal to binary conversion
```julia-repl
julia> dec2bin(10,7)
```
"""
function dec2bin(x,n)
	a=digits.(x,base=2,pad=n)  # left-lsb
# a=digits.(x,base=2,pad=num_bits)|> reverse  # left-msb
b=hcat(a...);
	return b
end

"""
Binary to decimal number conversion. Input can be 
-   binary strings, 
-   binary digits or 
-   a vector of binary string or digits
```julia-repl
julia> bin2dec([011,"0111",0111])
```
"""
bin2dec(u) =[parse(Int,string(u[:,ii]...),base=2) for ii in 1:size(u,2)]
bin2dec(x::AbstractVector)  = parse(Int,string(x...),base=2)
bin2dec(x::Int64)  = parse(Int,string(x),base=2)
bin2dec(x::AbstractString) = parse(Int,x,base=2)

"""
Pulse amplitude modulation (PAM) mapping. This is a type of digital modulation mapping used in Communication systems.

"""
function pam_encode(x,M)
	# M --> M-QAM
	n=Int(ceil(log2(M)/2))
	B,G,b,g=find_gray_matrix(n)
	u=digits(x,base=2,pad=n) |> reverse
	v=Int.(mod.(G*u,2))
	w=bin2dec(v)
	y=-sqrt(M)+1+2*w
	return y
end

"""
Generate Gray vectors
"""
function gen_gray(m)
	x=[0 1]
	for i in 2:m
		rev_x=reverse(x)
		x=[x rev_x.+2^(i-1)]
	end
	return x[:]
end

"""
Decimal to binary conversion
"""
function dec2bin(x,n)
	a=digits.(x,base=2,pad=n)  # left-lsb
# a=digits.(x,base=2,pad=num_bits)|> reverse  # left-msb
    b=hcat(a...);
	return b
end

"""
Recursive function to illustrate the reflection+shift property of Gray mapping.
### Arguents
* n - The iteration number `n ≥ 0`
* C - The decimal sequence of the gray mapped bits
* R - The reflected sequence (without shifting)
```julia-repl
julia> C,R = gray_recursion(4)
```
"""
function gray_recursion(n::Int)
    C= n < 1 ? [0] : vcat(gray_recursion(n - 1)[1:end], gray_recursion(n - 1)[end:-1:1] .+ Int(exp2(n-1)) )
	return C
end

"""
Reflected code. 
```julia-repl
julia>reflect_code(3)
[0,1,3,2,2,3,1,0]
```
"""
function reflect_code(n)
    n<1 ? [0] : vcat(gray_recursion(n-1),reverse(gray_recursion(n-1)))
end

"""
Recursive construction of binary Gray code digits.

Gray code ``g[n]`` can be recursively constructed as follows.
Start with `` g[1] = (0,1) = (g_{1},g_{2})`` 

```math
g[n+1] = 0g_{1},...,0g_{N−1},0g_{N},1g_{N},1g_{N−1},...,1g_{1}.
```

### Examples
```julia-repl
julia> gray(3)
3×8 Matrix{Int64}:
 0  0  0  0  1  1  1  1
 0  0  1  1  0  1  1  0
 0  1  1  0  1  1  0  0
```
"""
function gray(n)
    n < 2 ? [0 1] : hcat(vcat(0,gray(n-1)),vcat(1,reverse(gray(n-1))))
end

"""
Find the Givens embedding ``{}^{i}G_{j,k}(A)``
"""
function G(A,i,j,k)
	n = size(A,2)
	G = Matrix(1.0I(n).+0im)
	α = A[k,i]
	β = A[j,i]
	
	Γ0 = [ α' β'
	      -β  α ]
	N = norm(Γ0[1,:],2);
	Γ = Γ0./N

    # Embed the Givens matrix Γ in G
	G[k,k] = Γ[1,1]
	G[k,j] = Γ[1,2]
	G[j,k] = Γ[2,1]
	G[j,j] = Γ[2,2]
	
	return G,G*A
end

"""
Find the matrix which has a Givens matrix embedding ``{}^{i}G_{j,k}(A)``.

For a given matrix ``A``, a generic rotation matrix ``{}^{i}G_{j,k}(A)`` is generated. The matrix ``{}^{i}G_{j,k}(A)`` is such that it helps to selectively nullifies an element of matrix ``V={}^{i}G_{j,k}(A) A``. That is, ``V_{j,i}=0`` where ``V={}^{i}G_{j,k}(A) A``.  The fllowing Givens rotation matrix ``{}^{i}\\Gamma_{j,k}``

```math
\\begin{aligned}
{}^{i}\\Gamma_{j,k} &= \\begin{pmatrix} {}^{i}g_{k,k} & {}^{i}g_{k,j} \\\\ {}^{i}g_{j,k} & {}^{i}g_{j,j}\\end{pmatrix} \\\\
&=\\frac{1}{\\sqrt{\\lvert a_{j,i} \\rvert^{2}+ \\lvert a_{k,i} \\rvert^{2}}}\\begin{pmatrix} a_{k,i}^{*} & a_{j,i}^{*} \\\\ -a_{j,i} & {}^{i}a_{k,i}\\end{pmatrix}.
\\end{aligned}
```

is embedded in an identity matrix ``I(n)`` to produce,
```math
{}^{i}G_{j,k}(A) = \\begin{pmatrix} 
1 & 0 & \\ldots & \\ldots & \\ldots & \\ldots & \\ldots & \\ldots & 0 \\\\ 
0 & 1 & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\vdots \\\\
\\vdots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\vdots \\\\
0 & 0 & \\ddots & {\\color{red} {}^{i}g_{k,k}} & \\ddots &  {\\color{red} {}^{i}g_{k,j}} & \\ddots & \\ddots & \\vdots \\\\
\\vdots & \\ddots & \\ddots & \\ddots & \\ddots \\ddots & & \\ddots & \\ddots & \\vdots \\\\
0 & 0 & \\ddots & {\\color{red}  {}^{i}g_{j,k}} & \\ddots &  {\\color{red} 
 {}^{i}g_{j,j}} & \\ddots & \\ddots & \\vdots \\\\
\\vdots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\ddots & \\vdots \\\\
0 & 0 & \\ldots & \\ldots & \\ldots & \\ldots & \\ldots & 1 & 0  \\\\
0 & 0 & \\ldots & \\ldots & \\ldots & \\ldots & \\ldots & \\ldots & 1 
\\end{pmatrix}.
```
Essentially, ``{}^{i}G_{j,k}(A)`` is a modified identity matrix such that four non trivial elements are taken from the givens rotation matrix ``{}^{i}\\Gamma_{j,k}``.


"""
function GivensG(A,i,j,k)
	n = size(A,2)
	G = Matrix(1.0I(n).+0im)
	α = A[k,i]
	β = A[j,i]
	
	Γ0 = [ α' β'
	      -β  α ]
	N = norm(Γ0[1,:],2);
	Γ = Γ0./N

    # Embed the Givens matrix Γ in G
	G[k,k] = Γ[1,1]
	G[k,j] = Γ[1,2]
	G[j,k] = Γ[2,1]
	G[j,j] = Γ[2,2]
	
	return G,G*A
end

"""
Given unitary matrux ``U(n) \\in S(2^n)``, it uses a a repeated Givens rotations to get a to level matrix as follows.

```math
\\prod_{j=1}^{n-2}{ {}^{1}G_{n-j,n-j-1} U(n)}  = \\begin{pmatrix} 1 & 0 \\\\ 0 & U(n-1)\\end{pmatrix}
```

### Parameters 
* U -- Input. Unitary matrix of size ``2^n``
* V -- Output unitary matrix in two level form `[1 0;0 U']` form where `U'` is a unitary matrix of size ``2^{n-1}``.
* GG -- The ``n-2`` sequence of Given matrices (in augmented form) ``[{}^{1}G_{n}\\lvert{}^{1}G_{n-1}\\lvert\\ldots\\lvert{}^{1}G_{2}]``

### Examples
```julia-repl
julia> level2unitary(U)
```
"""
function level2unitary(U)
	n=size(U,2)
	V=I(n)*U;
	GG=Matrix(I(n))
	for i=1:1
		for j=0:n-1-i
			G0,V=GivensG(V,i,n-j,n-j-1)
			GG=[GG;G0]
		end
	end
	return GG[n+1:end,:],V
end

"""
Decomposition of aribtrary unitary matrix ``U(n) \\in S(2^n)``, as a cascade of two level Givens matrices.

Using the following property,
```math
\\prod_{j=0}^{n-1}{ {}^{1}G_{n-j,n-j-1} U(n)}  = \\begin{pmatrix} 1 & 0 \\\\ 0 & U(n-1)\\end{pmatrix}
```

```math
\\prod_{i=1}^{n-1} \\prod_{j=0}^{n-1-i}{ {}^{i}G_{n-j,n-j-1} U(n)}  = \\begin{pmatrix} 1 & \\ldots & 0 \\\\ \\vdots & \\ddots & \\vdots \\\\ 0 & \\ldots & \\det(U)\\end{pmatrix}.
```
There are ``2^{n-1}(2^{n}-1)`` number of unitary two-level matrices (each of them formed by embedding a Givens rotaton matrix into indentity matrix). Note that ``\\sum_{i=1}^{2^{n}}{N-i}=2^{n-1}(2^{n}-1)``.

### Parameters 
* U -- Input. Unitary matrix of size ``2^n``
* V -- Output unitary matrix in two level form `[1 0;0 U']` form where `U'` is a unitary matrix of size ``2^{n-1}``.
* Gm -- The ``(n-2)(n-1)`` sequence of Given matrices (in augmented form) ``[{}^{1}G_{n}\\lvert{}^{1}G_{n-1}\\lvert\\ldots\\lvert{}^{1}G_{2}\\lvert {}^{2}G_{n}\\lvert{}^{2}G_{n-1}\\lvert\\ldots\\lvert{}^{2}G_{3}\\lvert\\ldots\\lvert{}^{n-2}G_{n}\\lvert{}^{n-2}G_{n-1}\\lvert{}^{n-1}G_{n}]``
* Gs -- The ``(n-2)(n-1)`` sequence of Given matrices (in augmented form) left to right ``[ {}^{n-1}G_{n} \\lvert {}^{n-2}G_{n-1} \\lvert {}^{n-2}G_{n} \\lvert \\ldots \\lvert {}^{2}G_{3} \\lvert \\ldots \\lvert {}^{2}G_{n-1} \\lvert {}^{2}G_{n}  \\lvert{}^{1}G_{2} \\lvert \\ldots \\lvert {}^{1}G_{n-1} \\lvert {}^{1}G_{n} ]``  



### Examples
```julia-repl
julia> using LinearAlgebra
julia> N=4;A=rand(N,N)+im*rand(N,N);S=svd(A);U=S.U
julia> Gn(U)
```
"""
function Gn(A)
	n=size(A,2)
	V=I(n)*A;
	Gm=Matrix(I(n))
	Gs=Matrix(I(n))
	for i=1:n-1
		for j = 0:n-1-i
			G1,V = GivensG(V,i,n-j,n-j-1)
			Gm   = [Gm ; G1]
			Gs   = [G1 ; Gs]
		end
	end
	return V,Gm[n+1:end,:],Gs[1:end-n,:]
	
end

"""
For a given unitary matrix ``U``, it finds the cascaded rotation matrix ``C`` such that ``C\\times U =I``, except for the diagonal element of the ``n``th element which is ``\\det(U)``. The matrix ``C`` is obtained by the cascade of several (i.e., ``2^{n-1}(2^{n}-1)``) two level Givens rotation matrices, namely,

```math
C  = \\prod_{i=1}^{n-1} \\prod_{j=0}^{n-1-i}{ {}^{i}G_{n-j,n-j-1}}
```

```math
\\prod_{i=1}^{n-1} \\prod_{j=0}^{n-1-i}{ {}^{i}G_{n-j,n-j-1} U(n)}  = \\begin{pmatrix} 1 & \\ldots & 0 \\\\ \\vdots & \\ddots & \\vdots \\\\ 0 & \\ldots & \\det(U)\\end{pmatrix}.
```


### Parameters 
* U -- Input. Unitary matrix of size ``n``
* Ic --  Identity matix with the exception that the last diagonal entry (``n``th diagonal element) which is ``\\det(U)``.
* C -- The cascaded rotation matrix ``{}^{n-1}G_{n} \\times {}^{n-2}G_{n-1}  \\times {}^{n-2}G_{n}  \\times \\ldots  \\times {}^{2}G_{3}  \\times \\ldots  \\times {}^{2}G_{n-1}  \\times {}^{2}G_{n}   \\times {}^{1}G_{2}  \\times \\ldots  \\times {}^{1}G_{n-1}  \\times {}^{1}G_{n} ``  



### Examples
```julia-repl
julia> using LinearAlgebra
julia> N=4;A=rand(N,N)+im*rand(N,N);S=svd(A);U=S.U
4×4 Matrix{ComplexF64}:
  -0.4-0.06im   0.23-0.73im  -0.03-0.14im   0.39-0.27im
 -0.61-0.32im   0.07+0.06im   0.32+0.02im  -0.64-0.08im
 -0.38-0.33im  -0.48+0.38im  -0.07-0.4im    0.46+0.07im
  -0.2-0.26im  -0.09-0.14im   -0.7+0.47im  -0.09+0.38im
julia> Gc,Ic=Gcascaded(U)
4×4 Matrix{ComplexF64}:
 1.0+0.0im   0.0+0.0im  -0.0+0.0im   0.0-0.0im
 0.0+0.0im   1.0-0.0im   0.0+0.0im  -0.0-0.0im
 0.0-0.0im  -0.0+0.0im   1.0+0.0im   0.0+0.0im
 0.0-0.0im  -0.0+0.0im   0.0-0.0im   0.4-0.9im
julia> det(A)
0.4166084175706718 - 0.9090860390575042im
```
"""
function Gcascaded(U)
	Gv=Gn(U)[3]
	Gc=Matrix(I(4))
	m=Int(size(Gv,1)/4)
	for q in 0:m-1
		Gc=Gc*Gv[1+q*4:4*(q+1),:]
	end
	Ic=Gc*U
	return Gc,Ic
end

"""
For ``n`` bits, with the corresponding decimal sequence ``x=0,1,2,\\ldots,2^{n-1}``, find the gray ordering sequence using the bitwise `XOR` logic. Namely.
``gπ= x ⊻ \\lfloor x/2 \\rfloor = x \\oplus \\lfloor x/2 \\rfloor`` where ``⊻`` is the bitwise `XOR` operation.

```julia-repl
julia> gray_ordering(3)
[0,1,3,2,6,7,5,4]
```
"""
function gray_ordering(n)
	N=2^n
	x=Int.(vcat(range(start=0,stop=N-1,step=1)))
	y = Int.(floor.(x/2))
	gπ=[x[i] ⊻ y[i] for i=1:N]
	return gπ
end

""" Converts a gray mapped decimal number to simple binary mapped decimal number. This routine follows the inverse mapping of the function `gray_ordering()`
### Parameters
* ``n`` The number of digits in the equivalent binary mapping
* ``g`` The gray mapped decimal number
* ``d`` The simple binary mapped decimal number

```julia-repl
julia> G2D.(0:7,3)
[0,1,3,2,7,6,4,5]

```

"""
function G2D(g::Int,n)
	xx=gray_ordering(n) .+1
	a=collect(1:2^n)
	p=invperm(xx)
	d=p[g+1]
	return d-1
end


"""
Quantum circuit decomposition. For nan arbitrary unitary matix for ``n`` qubits. Arbitrary quantum circuit abstracted by unitary matrix ``U``  decomposed by ``2^{n-1}2^{n}`` unitary two-level matrices, each of which corresponds ``{}^{i}\\Gamma_{j,k}``. The program produce the ``\\Gamma`` matrix and the coefficients ``i,j,k``. The quantum circuit of this decomposition can be visualized as a cascade (from left to right) of this matrix.
"""
function sequenceΓ(k::Int64)
 	n=Int(2^k)
	N=Int(n*(n-1)/2)
	idxI=-100
	idxJ=-100
	m=0
	for i=1:n-1
		for j = n:-1:(i+1)
			idxI=[idxI n-i]
			idxJ=[idxJ j]

		end
	end
	dI,dJ=Int.(idxI[2:end] .+0),Int.(idxJ[2:end] .+0)
 	gV=gray_ordering(k)

	gJ2 = gV[dJ .- 1] .+1

	gI=gV[-dI.+0 .+ n] .+1 

	gJ= gV[dJ] .+1


	ii=reverse(gI)
	jj=reverse(gJ)
	kk=reverse(gJ2)
	Γ = hcat(ii,jj,kk)'
	return ii,jj,kk,Γ
	
end

end
