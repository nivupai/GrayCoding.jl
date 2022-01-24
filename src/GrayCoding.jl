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

end
