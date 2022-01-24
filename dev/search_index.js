var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GrayCoding","category":"page"},{"location":"#GrayCoding","page":"Home","title":"GrayCoding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GrayCoding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GrayCoding]","category":"page"},{"location":"#GrayCoding.GrayMatrix","page":"Home","title":"GrayCoding.GrayMatrix","text":"Generate Encoding and Decoding matrices for Gray Codes of alphabet.\n\njulia> G,B,g,b=GrayMatrix(4, 2);\njulia> G\n    4×4 Matrix{Int64}:\n    1  0  0  0\n    1  1  0  0\n    1  1  1  0\n    1  1  1  1\n    julia> B\n    4×4 Matrix{Int64}:\n    1  0  0  0\n    1  1  0  0\n    0  1  1  0\n    0  0  1  1\n    julia> g \n    4×16 Matrix{Int64}:\n    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n    0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1\n    0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1\n    0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1\n    julia> b \n    4×16 Matrix{Int64}:\n    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n    0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0\n    0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0\n    0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0\n\n\n\n\n\n","category":"function"},{"location":"#GrayCoding.bin2dec-Tuple{Any}","page":"Home","title":"GrayCoding.bin2dec","text":"Binary to decimal number conversion. Input can be \n\nbinary strings, \nbinary digits or \na vector of binary string or digits\n\njulia> bin2dec([011,\"0111\",0111])\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.dec2bin-Tuple{Any, Any}","page":"Home","title":"GrayCoding.dec2bin","text":"Decimal to binary conversion\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.dnamatrix-Tuple{}","page":"Home","title":"GrayCoding.dnamatrix","text":"Plot the DNA codon matrix\n\njulia> dnamatrix()\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.gen_gray-Tuple{Any}","page":"Home","title":"GrayCoding.gen_gray","text":"Generate Gray vectors\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.matrixplot-Tuple{Any}","page":"Home","title":"GrayCoding.matrixplot","text":"Plots a matrix into a 2D with labels. Optional arguments including colors\n\njulia> using Random;\njulia> A= rand(0:9,10,10);\njulia> matrixplot(A)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.pam_encode-Tuple{Any, Any}","page":"Home","title":"GrayCoding.pam_encode","text":"Pulse amplitude modulation (PAM) mapping. This is a type of digital modulation mapping used in Communication systems.\n\n\n\n\n\n","category":"method"}]
}