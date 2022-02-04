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
\\prod_{i=1}^{n-1} \\prod_{j=0}^{n-1-i}{ {}^{i}G_{n-j,n-j-1} U(n)}  = \\begin{pmatrix} 1 & \\ldots & 0 \\\\ \\vdots & \\ddots & \\vdots \\\\ 0 & \\ldots & \\det()\\end{pmatrix}
```


### Parameters 
* U -- Input. Unitary matrix of size ``2^n``
* V -- Output unitary matrix in two level form `[1 0;0 U']` form where `U'` is a unitary matrix of size ``2^{n-1}``.
* Gm -- The ``(n-2)(n-1)`` sequence of Given matrices (in augmented form) ``[{}^{1}G_{n}\\lvert{}^{1}G_{n-1}\\lvert\\ldots\\lvert{}^{1}G_{2}\\lvert {}^{2}G_{n}\\lvert{}^{2}G_{n-1}\\lvert\\ldots\\lvert{}^{2}G_{3}\\lvert\\ldots\\lvert{}^{n-2}G_{n}\\lvert{}^{n-2}G_{n-1}\\lvert{}^{n-1}G_{n}]``

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
	for i=1:n-1
		for j = 0:n-1-i
			G1,V = GivensG(V,i,n-j,n-j-1)
			Gm   = [Gm ; G1]
		end
	end
	return Gm[n+1:end,:],V
	
end

end
