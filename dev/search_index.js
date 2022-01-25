var documenterSearchIndex = {"docs":
[{"location":"tutorials/#Examples","page":"Tutorials","title":"Examples","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"TBD g=Gb and b=Bg, where G is a Jordan matrix, which is ","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"julia> n,q=4,2\njulia> GrayCoding.GrayMatrix(n,q)\n4×4 Matrix{Int64}:\n 1  0  0  0\n 1  1  0  0\n 1  1  1  0\n 1  1  1  1\n4×4 Matrix{Int64}:\n 1  0  0  0\n 1  1  0  0\n 0  1  1  0\n 0  0  1  1\n4×16 Matrix{Int64}:\n 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n 0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1\n 0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1\n 0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1\n4×16 Matrix{Int64}:\n 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n 0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0\n 0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0\n 0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0\n\njulia> G,B,g,b=GrayCoding.GrayMatrix(10,5);\njulia> G\n10×10 Matrix{Int64}:\n 1  0  0  0  0  0  0  0  0  0\n 1  1  0  0  0  0  0  0  0  0\n 1  1  1  0  0  0  0  0  0  0\n 1  1  1  1  0  0  0  0  0  0\n 1  1  1  1  1  0  0  0  0  0\n 1  1  1  1  1  1  0  0  0  0\n 1  1  1  1  1  1  1  0  0  0\n 1  1  1  1  1  1  1  1  0  0\n 1  1  1  1  1  1  1  1  1  0\n 1  1  1  1  1  1  1  1  1  1\n julia>B\n 10×10 Matrix{Int64}:\n 1  0  0  0  0  0  0  0  0  0\n 4  1  0  0  0  0  0  0  0  0\n 0  4  1  0  0  0  0  0  0  0\n 0  0  4  1  0  0  0  0  0  0\n 0  0  0  4  1  0  0  0  0  0\n 0  0  0  0  4  1  0  0  0  0\n 0  0  0  0  0  4  1  0  0  0\n 0  0  0  0  0  0  4  1  0  0\n 0  0  0  0  0  0  0  4  1  0\n 0  0  0  0  0  0  0  0  4  1","category":"page"},{"location":"applications/#Applications","page":"Applications of Gray Codes","title":"Applications","text":"","category":"section"},{"location":"applications/","page":"Applications of Gray Codes","title":"Applications of Gray Codes","text":"Digital Modulation Schemes:\nDNA Codon mapping\nDigital Electronics/Counters\nMusic/Puzzles","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GrayCoding","category":"page"},{"location":"#GrayCoding","page":"Home","title":"GrayCoding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GrayCoding.","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for GrayCoding!","category":"page"},{"location":"#What-is-GrayCocoding.jl?","page":"Home","title":"What is GrayCocoding.jl?","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"GrayCoding is a formal Linear Algebraic framework for q-ry Gray Code. Encoding and Decooding of Gray codes can be treated as a special case of algebraic block coding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Encoding: textbfg=G textbfb\nDecoding: textbfb=B textbfg","category":"page"},{"location":"","page":"Home","title":"Home","text":"tip: Tip\nThis is still under active devlopment.","category":"page"},{"location":"#Resources-for-getting-started","page":"Home","title":"Resources for getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are few ways to get started with GrayCoding:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Read TBD.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GrayCoding]","category":"page"},{"location":"#GrayCoding.GrayMatrix","page":"Home","title":"GrayCoding.GrayMatrix","text":"Generate Encoding and Decoding matrices for Gray Codes of alphabet.\n\njulia> G,B,g,b=GrayMatrix(4, 2);\njulia> G\n    4×4 Matrix{Int64}:\n    1  0  0  0\n    1  1  0  0\n    1  1  1  0\n    1  1  1  1\n    julia> B\n    4×4 Matrix{Int64}:\n    1  0  0  0\n    1  1  0  0\n    0  1  1  0\n    0  0  1  1\n    julia> g \n    4×16 Matrix{Int64}:\n    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n    0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1\n    0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1\n    0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1\n    julia> b \n    4×16 Matrix{Int64}:\n    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n    0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0\n    0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0\n    0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0\n\n\n\n\n\n","category":"function"},{"location":"#GrayCoding.bin2dec-Tuple{Any}","page":"Home","title":"GrayCoding.bin2dec","text":"Binary to decimal number conversion. Input can be \n\nbinary strings, \nbinary digits or \na vector of binary string or digits\n\njulia> bin2dec([011,\"0111\",0111])\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.dec2bin-Tuple{Any, Any}","page":"Home","title":"GrayCoding.dec2bin","text":"Decimal to binary conversion\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.dnamatrix-Tuple{}","page":"Home","title":"GrayCoding.dnamatrix","text":"Plot the DNA codon matrix\n\njulia> dnamatrix()\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.gen_gray-Tuple{Any}","page":"Home","title":"GrayCoding.gen_gray","text":"Generate Gray vectors\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.matrixplot-Tuple{Any}","page":"Home","title":"GrayCoding.matrixplot","text":"Plots a matrix into a 2D with labels. Optional arguments including colors\n\njulia> using Random;\njulia> A= rand(0:9,10,10);\njulia> matrixplot(A)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.pam_encode-Tuple{Any, Any}","page":"Home","title":"GrayCoding.pam_encode","text":"Pulse amplitude modulation (PAM) mapping. This is a type of digital modulation mapping used in Communication systems.\n\n\n\n\n\n","category":"method"},{"location":"algebra/#Algebraic-framework-of-Gray-Codes","page":"Algebra of Gray Codes","title":"Algebraic framework of Gray Codes","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"The classical algorithmic procedure of encoding and decoding are as follows:","category":"page"},{"location":"algebra/#Encoding","page":"Algebra of Gray Codes","title":"Encoding","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"q-ry digits d to gray digits g conversion.","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"g_i = begincases  d_i   textif  modleft(displaystylesum_j=1^i-1g_j2right)=0   q-1-d_i   textif modleft(displaystyle sum_j=1^i-1g_j2right)=1 endcases","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"and g_1 = d_1.","category":"page"},{"location":"algebra/#Decoding","page":"Algebra of Gray Codes","title":"Decoding","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"d_i = begincases  g_i   textif  modleft(displaystylesum_j=1^i-1g_j2right)=0   q-1-g_i   textif modleft(displaystyle sum_j=1^i-1g_j2right)=1  endcases","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"and d_1 = g_1.","category":"page"},{"location":"algebra/#Linear-Algebraic-Formulation-(N.Rethnakar-2020)","page":"Algebra of Gray Codes","title":"Linear Algebraic Formulation (N.Rethnakar 2020)","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"Example of generator matrix G for binary to gray mapping is given by, `` G=\\begin{pmatrix} 1 &  0  &  0 &  0  &  0 &  0   \\\n 1   &  1  & 0 & 0 & 0 & 0   \\\n 0   &  1  & 1 & 0 & 0 & 0   \\\n 0   &  0  & 1 & 1 & 0 & 0   \\\n 0   &  0  & 0 & 1 & 1 & 0   \\\n 0   &  0  & 0 & 0 & 1 & 1 \\end{pmatrix} ``  ","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"The decoding matrix B=G^-1 is given by,","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"``   B=\\begin{pmatrix}   1   &   0  &  0   &   0   &  0  &  0   \\\n1 &  1 &   0  &   0  &  0    &  0   \\\n1  &  1  &   1 & 0  &  0  &   0   \\\n1 &  1 &   1 &  1  &  0   &   0   \\\n1 &   1 &  1 &   1  &   1   & 0   \\\n1 &  1 &   1  &   1  &  1  & 1 \\end{pmatrix} ``","category":"page"}]
}
