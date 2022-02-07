var documenterSearchIndex = {"docs":
[{"location":"tutorials/#Examples","page":"Tutorials","title":"Examples","text":"","category":"section"},{"location":"tutorials/#Recursive-construction","page":"Tutorials","title":"Recursive construction","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"(Image: )","category":"page"},{"location":"tutorials/#Recursive-procedure","page":"Tutorials","title":"Recursive procedure","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"(Image: )","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Reflect Cn-1, shift by q^n-1 and augment (TBD).","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"For n in mathbbN, positive integer and N = 2^n. A Gray code G_n is an tuple G_n = (X_1X_2X_N) which satisfies the following properties:","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"X_1 X_2   X_N are binary sequences (of length n) corresponding to the binary representation of the numbers 0 1 ldots  N  1, arranged in a specific order,\nFor any 0 le j le N-1, adjacent pairs X_jX_j+1 differ in only one position (i.e., only one of the n bits would differ),\nThe start and end sequences (i.e., sequences X_1 and X_N differ in just one position.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Gray code G_n can be recursively constructed as follows. Start with G_1 = (01) and for N=2^n n ge 1, Let G_n = left(X_1ldotsX_N1X_Nright), ","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"G_n+1 = left(0X_1ldots0X_N10X_N1X_N1X_N11X_1right)","category":"page"},{"location":"tutorials/#Illustration","page":"Tutorials","title":"Illustration","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"using Plots\nfunction plotmatrix(A;kwargs...)\n    a,b=size(A)\n    X = transpose(repeat(1:b, 1, a))[:]\n    Y = repeat(a:-1:1, b)[:]\n\tscatter(X,Y, marker_z = A[:], marker=:rect,markersize = 4,  color = :viridis,aspectratio=1,ylims=[0,size(G,1)+1],alpha=1,label=:none,colorkey=:none,axis=:none;kwargs...)\n\njulia> plotmatrix(gray(6));\njulia> plotmatrix(G,size=(800,400),color=:summer)\njulia> plotmatrix(G,size=(800,200),color=:summer,markersize=7,xlims=[1,size(G,2)+0],ylims=[1/2,size(G,1)-0])\nend","category":"page"},{"location":"tutorials/#Binary-Gray-Code-n4","page":"Tutorials","title":"Binary Gray Code n=4","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"(Image: )","category":"page"},{"location":"tutorials/#Binary-Gray-Code-n5","page":"Tutorials","title":"Binary Gray Code n=5","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"(Image: )","category":"page"},{"location":"tutorials/#Binary-Gray-Code-n6","page":"Tutorials","title":"Binary Gray Code n=6","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"(Image: )","category":"page"},{"location":"tutorials/#Linear-Algebraic-method","page":"Tutorials","title":"Linear Algebraic method","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"TBD g=Gb and b=Bg, where G is a Jordan matrix, which is ","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"julia> n,q=4,2\njulia> GrayCoding.GrayMatrix(n,q)\n4×4 Matrix{Int64}:\n 1  0  0  0\n 1  1  0  0\n 1  1  1  0\n 1  1  1  1\n4×4 Matrix{Int64}:\n 1  0  0  0\n 1  1  0  0\n 0  1  1  0\n 0  0  1  1\n4×16 Matrix{Int64}:\n 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n 0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1\n 0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1\n 0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1\n4×16 Matrix{Int64}:\n 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n 0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0\n 0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0\n 0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0\n\njulia> G,B,g,b=GrayCoding.GrayMatrix(10,5);\njulia> G\n10×10 Matrix{Int64}:\n 1  0  0  0  0  0  0  0  0  0\n 1  1  0  0  0  0  0  0  0  0\n 1  1  1  0  0  0  0  0  0  0\n 1  1  1  1  0  0  0  0  0  0\n 1  1  1  1  1  0  0  0  0  0\n 1  1  1  1  1  1  0  0  0  0\n 1  1  1  1  1  1  1  0  0  0\n 1  1  1  1  1  1  1  1  0  0\n 1  1  1  1  1  1  1  1  1  0\n 1  1  1  1  1  1  1  1  1  1\n julia>B\n 10×10 Matrix{Int64}:\n 1  0  0  0  0  0  0  0  0  0\n 4  1  0  0  0  0  0  0  0  0\n 0  4  1  0  0  0  0  0  0  0\n 0  0  4  1  0  0  0  0  0  0\n 0  0  0  4  1  0  0  0  0  0\n 0  0  0  0  4  1  0  0  0  0\n 0  0  0  0  0  4  1  0  0  0\n 0  0  0  0  0  0  4  1  0  0\n 0  0  0  0  0  0  0  4  1  0\n 0  0  0  0  0  0  0  0  4  1","category":"page"},{"location":"applications/#Applications","page":"List of Applications","title":"Applications","text":"","category":"section"},{"location":"applications/","page":"List of Applications","title":"List of Applications","text":"Digital Modulation Schemes:\nDNA Codon mapping\nQuantum Circuits and Gates\nDigital Electronics/Counters\nMusic/Puzzles","category":"page"},{"location":"quantum/#Quantum-Computing","page":"Quantum Algorithms and Circuits","title":"Quantum Computing","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"Broadly stated, the term quantum computation comprieses of the following three elements:","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"A register or a set of registers,\nA unitary matrix U, as an abstract representation of the quantum algorithm,\nMeasurements to extract the information of interest.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"As a mathematial abstraction, a quantum computation is the set mathcalHUMm, where H = mathbbC^2^n is the Hilbert space of an n-qubit register, U in Uleft(2^nright) represents the quantum algorithm and M_m is the set of measurement operators. The hardware circuitry along with equipment to control and manipulate the qubits is called a quantum computer.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"info: More on Quantum Computation Information\nQuantum Computing is a fascinating subject. In this, we cover only the very basic things to connect the Gray Code framework that we've devloped to Quantum circuits. The classic book on this subject by Michael A. Nielsen’s and Isaac L. Chuang titled “Quantum Computation and Quantum Information” is the go to place for more (Michael A. Nielsen and Isaac L. Chuang. Quantum Computation and Quantum Information. Cambridge University Press (2000)). ","category":"page"},{"location":"quantum/#Single-qubit-gates","page":"Quantum Algorithms and Circuits","title":"Single qubit gates","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"These are the simplest set of gates which take one qubit as an input, act upon (change the state) and produce a qubit at the output.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"A  generic 1-qubit quantum gate corresponds to a 2 times 2 unitary matrix, which has the following form:","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"U = e^imath theta beginpmatrix a  -b^*  b  a^* endpmatrix ","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"where ab in mathbbC such that lvert a rvert^2+ lvert b rvert^2 = 1, and alpha in mathbbR results in arbitrary rotation. This matrix is essentially a Givens rotation matrix Wikipedia. ","category":"page"},{"location":"quantum/#Examples","page":"Quantum Algorithms and Circuits","title":"Examples","text":"","category":"section"},{"location":"quantum/#X-gate.","page":"Quantum Algorithms and Circuits","title":"X gate.","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"This gate flips the state of the qubit. In other words, it changes the state of the qubit from avert 0 rangle + bvert 1 rangle to avert 1 rangle+ bvert 0 rangle. \nThe matrix representation forX gate is, sigma_x=beginpmatrix 0  1  1  0 endpmatrix.  Multiplying the vector representing the qubit to the matrix is equivalent to the gate operation. i.e., beginpmatrix b  a endpmatrix =beginpmatrix 0  1  1  0 endpmatrix beginpmatrix a  b endpmatrix ","category":"page"},{"location":"quantum/#Y-gate.","page":"Quantum Algorithms and Circuits","title":"Y gate.","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"This gate perform two flips (a bit flip and a phase flip) on the state of the qubit. In other words, it changes the state of the qubit from avert 0 rangle + bvert 1 rangle, to bvert 0 rangle - avert 1 rangle. \nThe matrix representation forY gate is, sigma_y=beginpmatrix 0  -1  1  0  endpmatrix.  Multiplying the vector representing the qubit to the matrix is equivalent to the gate operation. ","category":"page"},{"location":"quantum/#Z-gate.","page":"Quantum Algorithms and Circuits","title":"Z gate.","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"This gate performs a sign flip on the state of the qubit. In other words, it changes the state of the qubit from avert 0 rangle + bvert 1 rangle, to avert 0 rangle - vert 1 rangle. \nThe matrix representation for X gate is, sigma_z=beginpmatrix 1  0  0  -1 endpmatrix.  Multiplying the vector representing the qubit to the matrix is equivalent to the gate operation. ","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"The matrix representations sigma_x sigma_yand sigma_z are known as Pauli's matrices.","category":"page"},{"location":"quantum/#Hadamard-gate:-H-gate","page":"Quantum Algorithms and Circuits","title":"Hadamard gate: H gate","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"This gate works as follows: vert 0 rangle changes to frac1sqrt2(vert 0 rangle + vert 1 rangle), and the vert 1 rangle changes to frac1sqrt2(vert 0 rangle - vert 1 rangle). \nFor an example, avert 0 rangle + bvert 1 rangle changes to, fracasqrt2(vert 0 rangle + vert 1 rangle) +  fracbsqrt2(vert 0 rangle - vert 1 rangle). \nIt can be simplified to give, fraca+bsqrt2vert 0 rangle + fraca-bsqrt2 vert 1 rangle.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"Mathematically, the following transformation captures the essence of H gate.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"beginpmatrix fraca+bsqrt2  fraca-bsqrt2 endpmatrix = frac1sqrt2 beginpmatrix 1  1  1  -1 endpmatrix beginpmatrix a  b endpmatrix","category":"page"},{"location":"quantum/#Decomposition-of-Quantum-Gates","page":"Quantum Algorithms and Circuits","title":"Decomposition of Quantum Gates","text":"","category":"section"},{"location":"quantum/#Two-level-decomposition","page":"Quantum Algorithms and Circuits","title":"Two level decomposition","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"A 2-level unitary correspond to unitary operation that non-trivially perform on only 2 of the states. Any controlled 1-qubit  gate can be abstracted into a 2-level unitary matrix, e.g. for a single qubit gate U = beginpmatrix a  b  c  d endpmatrix.","category":"page"},{"location":"quantum/#Universal-Decomposition","page":"Quantum Algorithms and Circuits","title":"Universal Decomposition","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"Any arbitrary unitary gate acting on n-qubit can be implemented as a cascade of single-qubit and controlled-NOT (CNOT) gates.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"textbfUpsi rangle = beginpmatrix u_11  u_12  ldots  u_1N  u_11  u_12  ldots  u_1N endpmatrix","category":"page"},{"location":"quantum/#CNOT-gates","page":"Quantum Algorithms and Circuits","title":"CNOT gates","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"CNOT stands for Controlled-NOT, which is one of the key quantum logic gate. It is a two qubit gate. The gate flips the second qubit (called the target qubit) when the first gate (the control gate) is lvert 1 rangle, while the target qubit remain unchanged when the control gate is in state lvert 0 rangle.","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"beginaligned\n U_textCNOT = lvert 0 rangle otimes langle 0 rvert otimes beginpmatrix 1  0  0  1 endpmatrix + lvert 1 rangle otimes langle 1 rvert otimes beginpmatrix 1  0  0  1 endpmatrix\n                 = lvert 0 rangle langle 0 rvert otimes beginpmatrix 1  0  0  1 endpmatrix + lvert 1 rangle  langle 1 rvert otimes beginpmatrix 1  0  0  1 endpmatrix\n                 = beginpmatrix 1  0  0  0  0  1  0  0  0  0  0  1  0  0  1  0endpmatrix \n endaligned","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"In terms of Gray matrix, we can also simply express this as,","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"beginaligned\nbeginpmatrix lvert acutepsi_1 rangle  lvert acutepsi_2 rangle endpmatrix = G_2 beginpmatrix lvert psi_1 rangle  lvert psi_2 rangle endpmatrix \n= beginpmatrix 1  0  1  1 endpmatrix beginpmatrix lvert psi_1 rangle  lvert psi_2 rangle endpmatrix\nendaligned","category":"page"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"","category":"page"},{"location":"quantum/#The-simplest-CNOT-gate-is-the-single-qubit-controlled-CNOT-discussed-above,-which-can-be-explicitly-denoted-as-C{1}\\text{NOT}.-Generalization-of-this-to-multi-quibit-controlled-CNOT,-denoted-by-C{n-1}\\text{NOT}.","page":"Quantum Algorithms and Circuits","title":"The simplest CNOT gate is the single qubit controlled CNOT discussed above, which can be explicitly denoted as C^1textNOT. Generalization of this to multi quibit controlled CNOT, denoted by C^n-1textNOT.","text":"","category":"section"},{"location":"quantum/","page":"Quantum Algorithms and Circuits","title":"Quantum Algorithms and Circuits","text":"In quantum circuit design, applying a rotation for which the binary representations of i − 1 and j − 1 differ in a single bit can be accomplished by a single fully-controlled one-qubit rotation (a particular Givens rotation) and hence costs a small number of gates. All other rotations require a permutation of data before the rotation is applied and thus should be avoided.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GrayCoding","category":"page"},{"location":"#GrayCoding","page":"Home","title":"GrayCoding","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GrayCoding.","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Welcome to the documentation for GrayCoding!","category":"page"},{"location":"#What-is-GrayCocoding.jl?","page":"Home","title":"What is GrayCocoding.jl?","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"GrayCoding is a formal Linear Algebraic framework for q-ry Gray Code. Encoding and Decooding of Gray codes can be treated as a special case of algebraic block coding.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Encoding: textbfg=G textbfb\nDecoding: textbfb=B textbfg","category":"page"},{"location":"","page":"Home","title":"Home","text":"tip: Tip\nThis is still under active devlopment.","category":"page"},{"location":"#Resources-for-getting-started","page":"Home","title":"Resources for getting started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"There are few ways to get started with GrayCoding:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Read TBD.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GrayCoding]","category":"page"},{"location":"#GrayCoding.G-NTuple{4, Any}","page":"Home","title":"GrayCoding.G","text":"Find the Givens embedding ^iG_jk(A)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.G2D-Tuple{Int64, Any}","page":"Home","title":"GrayCoding.G2D","text":"Converts a gray mapped decimal number to simple binary mapped decimal number. This routine follows the inverse mapping of the function gray_ordering()\n\nParameters\n\nn The number of digits in the equivalent binary mapping\ng The gray mapped decimal number\nd The simple binary mapped decimal number\n\njulia> G2D.(0:7,3)\n[0,1,3,2,7,6,4,5]\n\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.Gcascaded-Tuple{Any}","page":"Home","title":"GrayCoding.Gcascaded","text":"For a given unitary matrix U, it finds the cascaded rotation matrix C such that Ctimes U =I, except for the diagonal element of the nth element which is det(U). The matrix C is obtained by the cascade of several (i.e., 2^n-1(2^n-1)) two level Givens rotation matrices, namely,\n\nC  = prod_i=1^n-1 prod_j=0^n-1-i ^iG_n-jn-j-1\n\nprod_i=1^n-1 prod_j=0^n-1-i ^iG_n-jn-j-1 U(n)  = beginpmatrix 1  ldots  0  vdots  ddots  vdots  0  ldots  det(U)endpmatrix\n\nParameters\n\nU – Input. Unitary matrix of size n\nIc –  Identity matix with the exception that the last diagonal entry (nth diagonal element) which is det(U).\nC – The cascaded rotation matrix ^n-1G_n times ^n-2G_n-1  times ^n-2G_n  times ldots  times ^2G_3  times ldots  times ^2G_n-1  times ^2G_n   times ^1G_2  times ldots  times ^1G_n-1  times ^1G_n  \n\nExamples\n\njulia> using LinearAlgebra\njulia> N=4;A=rand(N,N)+im*rand(N,N);S=svd(A);U=S.U\n4×4 Matrix{ComplexF64}:\n  -0.4-0.06im   0.23-0.73im  -0.03-0.14im   0.39-0.27im\n -0.61-0.32im   0.07+0.06im   0.32+0.02im  -0.64-0.08im\n -0.38-0.33im  -0.48+0.38im  -0.07-0.4im    0.46+0.07im\n  -0.2-0.26im  -0.09-0.14im   -0.7+0.47im  -0.09+0.38im\njulia> Gc,Ic=Gcascaded(U)\n4×4 Matrix{ComplexF64}:\n 1.0+0.0im   0.0+0.0im  -0.0+0.0im   0.0-0.0im\n 0.0+0.0im   1.0-0.0im   0.0+0.0im  -0.0-0.0im\n 0.0-0.0im  -0.0+0.0im   1.0+0.0im   0.0+0.0im\n 0.0-0.0im  -0.0+0.0im   0.0-0.0im   0.4-0.9im\njulia> det(A)\n0.4166084175706718 - 0.9090860390575042im\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.GivensG-NTuple{4, Any}","page":"Home","title":"GrayCoding.GivensG","text":"Find the matrix which has a Givens matrix embedding ^iG_jk(A).\n\nFor a given matrix A, a generic rotation matrix ^iG_jk(A) is generated. The matrix ^iG_jk(A) is such that it helps to selectively nullifies an element of matrix V=^iG_jk(A) A. That is, V_ji=0 where V=^iG_jk(A) A.  The fllowing Givens rotation matrix ^iGamma_jk\n\nbeginaligned\n^iGamma_jk = beginpmatrix ^ig_kk  ^ig_kj  ^ig_jk  ^ig_jjendpmatrix \n=frac1sqrtlvert a_ji rvert^2+ lvert a_ki rvert^2beginpmatrix a_ki^*  a_ji^*  -a_ji  ^ia_kiendpmatrix\nendaligned\n\nis embedded in an identity matrix I(n) to produce,\n\n^iG_jk(A) = beginpmatrix \n1  0  ldots  ldots  ldots  ldots  ldots  ldots  0  \n0  1  ddots  ddots  ddots  ddots  ddots  ddots  vdots \nvdots  ddots  ddots  ddots  ddots  ddots  ddots  ddots  vdots \n0  0  ddots  colorred ^ig_kk  ddots   colorred ^ig_kj  ddots  ddots  vdots \nvdots  ddots  ddots  ddots  ddots ddots   ddots  ddots  vdots \n0  0  ddots  colorred  ^ig_jk  ddots   colorred \n ^ig_jj  ddots  ddots  vdots \nvdots  ddots  ddots  ddots  ddots  ddots  ddots  ddots  vdots \n0  0  ldots  ldots  ldots  ldots  ldots  1  0  \n0  0  ldots  ldots  ldots  ldots  ldots  ldots  1 \nendpmatrix\n\nEssentially, ^iG_jk(A) is a modified identity matrix such that four non trivial elements are taken from the givens rotation matrix ^iGamma_jk.\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.Gn-Tuple{Any}","page":"Home","title":"GrayCoding.Gn","text":"Decomposition of aribtrary unitary matrix U(n) in S(2^n), as a cascade of two level Givens matrices.\n\nUsing the following property,\n\nprod_j=0^n-1 ^1G_n-jn-j-1 U(n)  = beginpmatrix 1  0  0  U(n-1)endpmatrix\n\nprod_i=1^n-1 prod_j=0^n-1-i ^iG_n-jn-j-1 U(n)  = beginpmatrix 1  ldots  0  vdots  ddots  vdots  0  ldots  det(U)endpmatrix\n\nThere are 2^n-1(2^n-1) number of unitary two-level matrices (each of them formed by embedding a Givens rotaton matrix into indentity matrix). Note that sum_i=1^2^nN-i=2^n-1(2^n-1).\n\nParameters\n\nU – Input. Unitary matrix of size 2^n\nV – Output unitary matrix in two level form [1 0;0 U'] form where U' is a unitary matrix of size 2^n-1.\nGm – The (n-2)(n-1) sequence of Given matrices (in augmented form) ^1G_nlvert^1G_n-1lvertldotslvert^1G_2lvert ^2G_nlvert^2G_n-1lvertldotslvert^2G_3lvertldotslvert^n-2G_nlvert^n-2G_n-1lvert^n-1G_n\nGs – The (n-2)(n-1) sequence of Given matrices (in augmented form) left to right  ^n-1G_n lvert ^n-2G_n-1 lvert ^n-2G_n lvert ldots lvert ^2G_3 lvert ldots lvert ^2G_n-1 lvert ^2G_n  lvert^1G_2 lvert ldots lvert ^1G_n-1 lvert ^1G_n   \n\nExamples\n\njulia> using LinearAlgebra\njulia> N=4;A=rand(N,N)+im*rand(N,N);S=svd(A);U=S.U\njulia> Gn(U)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.GrayMatrix","page":"Home","title":"GrayCoding.GrayMatrix","text":"Generate Encoding and Decoding matrices for Gray Codes of alphabet.\n\njulia> G,B,g,b=GrayMatrix(4, 2);\njulia> G\n    4×4 Matrix{Int64}:\n    1  0  0  0\n    1  1  0  0\n    1  1  1  0\n    1  1  1  1\n    julia> B\n    4×4 Matrix{Int64}:\n    1  0  0  0\n    1  1  0  0\n    0  1  1  0\n    0  0  1  1\n    julia> g \n    4×16 Matrix{Int64}:\n    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n    0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1\n    0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1\n    0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1\n    julia> b \n    4×16 Matrix{Int64}:\n    0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1\n    0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0\n    0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0\n    0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0\n\n\n\n\n\n","category":"function"},{"location":"#GrayCoding.bin2dec-Tuple{Any}","page":"Home","title":"GrayCoding.bin2dec","text":"Binary to decimal number conversion. Input can be \n\nbinary strings, \nbinary digits or \na vector of binary string or digits\n\njulia> bin2dec([011,\"0111\",0111])\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.dec2bin-Tuple{Any, Any}","page":"Home","title":"GrayCoding.dec2bin","text":"Decimal to binary conversion\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.dnamatrix-Tuple{}","page":"Home","title":"GrayCoding.dnamatrix","text":"Plot the DNA codon matrix\n\njulia> dnamatrix()\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.gen_gray-Tuple{Any}","page":"Home","title":"GrayCoding.gen_gray","text":"Generate Gray vectors\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.gray-Tuple{Any}","page":"Home","title":"GrayCoding.gray","text":"Recursive construction of binary Gray code digits.\n\nGray code gn can be recursively constructed as follows. Start with g1 = (01) = (g_1g_2) \n\ngn+1 = 0g_10g_N10g_N1g_N1g_N11g_1\n\nExamples\n\njulia> gray(3)\n3×8 Matrix{Int64}:\n 0  0  0  0  1  1  1  1\n 0  0  1  1  0  1  1  0\n 0  1  1  0  1  1  0  0\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.gray_ordering-Tuple{Any}","page":"Home","title":"GrayCoding.gray_ordering","text":"For n bits, with the corresponding decimal sequence x=012ldots2^n-1, find the gray ordering sequence using the bitwise XOR logic. Namely. gπ= x  lfloor x2 rfloor = x oplus lfloor x2 rfloor where  is the bitwise XOR operation.\n\njulia> gray_ordering(3)\n[0,1,3,2,6,7,5,4]\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.gray_recursion-Tuple{Int64}","page":"Home","title":"GrayCoding.gray_recursion","text":"Recursive function to illustrate the reflection+shift property of Gray mapping.\n\nArguents\n\nn - The iteration number n ≥ 0\nC - The decimal sequence of the gray mapped bits\nR - The reflected sequence (without shifting)\n\njulia> C,R = gray_recursion(4)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.level2unitary-Tuple{Any}","page":"Home","title":"GrayCoding.level2unitary","text":"Given unitary matrux U(n) in S(2^n), it uses a a repeated Givens rotations to get a to level matrix as follows.\n\nprod_j=1^n-2 ^1G_n-jn-j-1 U(n)  = beginpmatrix 1  0  0  U(n-1)endpmatrix\n\nParameters\n\nU – Input. Unitary matrix of size 2^n\nV – Output unitary matrix in two level form [1 0;0 U'] form where U' is a unitary matrix of size 2^n-1.\nGG – The n-2 sequence of Given matrices (in augmented form) ^1G_nlvert^1G_n-1lvertldotslvert^1G_2\n\nExamples\n\njulia> level2unitary(U)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.matrixplot-Tuple{Any}","page":"Home","title":"GrayCoding.matrixplot","text":"Plots a matrix into a 2D with labels. Optional arguments including colors\n\njulia> using Random;\njulia> A= rand(0:9,10,10);\njulia> matrixplot(A)\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.pam_encode-Tuple{Any, Any}","page":"Home","title":"GrayCoding.pam_encode","text":"Pulse amplitude modulation (PAM) mapping. This is a type of digital modulation mapping used in Communication systems.\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.quantumΓ-NTuple{4, Any}","page":"Home","title":"GrayCoding.quantumΓ","text":"Given a matrix A (usually a unitary matrix) and indices (ijk), find the Gamma Givens rotation matrix and the Givens matrix embedding ^iG_jk(A) matrix.\n\nFor a given matrix A, a generic rotation matrix ^iG_jk(A) is generated. The matrix ^iG_jk(A) is such that it helps to selectively nullifies an element of matrix V=^iG_jk(A) A. That is, V_ji=0 where V=^iG_jk(A) A.  The fllowing Givens rotation matrix ^iGamma_jk\n\nbeginaligned\n^iGamma_jk = beginpmatrix ^ig_kk  ^ig_kj  ^ig_jk  ^ig_jjendpmatrix \n=frac1sqrtlvert a_ji rvert^2+ lvert a_ki rvert^2beginpmatrix a_ki^*  a_ji^*  -a_ji  ^ia_kiendpmatrix\nendaligned\n\nis embedded in an identity matrix I(n) to produce,\n\n^iG_jk(A) = beginpmatrix \n1  0  ldots  ldots  ldots  ldots  ldots  ldots  0  \n0  1  ddots  ddots  ddots  ddots  ddots  ddots  vdots \nvdots  ddots  ddots  ddots  ddots  ddots  ddots  ddots  vdots \n0  0  ddots  colorred ^ig_kk  ddots   colorred ^ig_kj  ddots  ddots  vdots \nvdots  ddots  ddots  ddots  ddots ddots   ddots  ddots  vdots \n0  0  ddots  colorred  ^ig_jk  ddots   colorred \n ^ig_jj  ddots  ddots  vdots \nvdots  ddots  ddots  ddots  ddots  ddots  ddots  ddots  vdots \n0  0  ldots  ldots  ldots  ldots  ldots  1  0  \n0  0  ldots  ldots  ldots  ldots  ldots  ldots  1 \nendpmatrix\n\nEssentially, ^iG_jk(A) is a modified identity matrix such that four non trivial elements are taken from the givens rotation matrix ^iGamma_jk.\n\njulia> using LinearAlgebra\njulia> using GrayCoding\njulia>n=2;N=2^n;A=rand(N,N)+im*rand(N,N);S=svd(A);U=S.U\n4×4 Matrix{ComplexF64}:\n -0.365903-0.405021im   0.442293-0.0769938im  …   0.115307-0.288609im\n -0.285173-0.35669im   -0.671764+0.0698449im     -0.384583+0.295428im\n -0.196831-0.611652im  -0.154487+0.0160399im      0.379159-0.121825im\n -0.177839-0.221435im   0.536228-0.175044im      -0.338822+0.62835im\njulia> i,j,k=1,2,4\njulia> Γ,G,GA=quantumΓ(U,i,j,k);\njulia> round.(quantumΓ(S.U,1,2,4)[1],digits=1)\n2×2 Matrix{ComplexF64}:\n -0.3+0.4im  -0.5+0.7im\n  0.5+0.7im  -0.3-0.4im\njulia> round.(quantumΓ(S.U,1,2,4)[2],digits=1)\n4×4 Matrix{ComplexF64}:\n 1.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im\n 0.0+0.0im  -0.3-0.4im  0.0+0.0im   0.5+0.7im\n 0.0+0.0im   0.0+0.0im  1.0+0.0im   0.0+0.0im\n 0.0+0.0im  -0.5+0.7im  0.0+0.0im  -0.3+0.4im\njulia> round.(quantumΓ(S.U,1,2,4)[3],digits=1)\n4×4 Matrix{ComplexF64}:\n -0.4-0.4im   0.4-0.1im   0.6+0.1im   0.1-0.3im\n  0000000     0.7+0.5im  -0.2-0.4im  -0.3+0.2im\n -0.2-0.6im  -0.2+0.0im  -0.4-0.5im   0.4-0.1im\n  0.5+0.0im   0.2-0.2im  -0.2-0.1im  -0.1-0.8im\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.reflect_code-Tuple{Any}","page":"Home","title":"GrayCoding.reflect_code","text":"Reflected code. \n\njulia>reflect_code(3)\n[0,1,3,2,2,3,1,0]\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.sequenceΓ-Tuple{Int64}","page":"Home","title":"GrayCoding.sequenceΓ","text":"Quantum circuit decomposition. For nan arbitrary unitary matix for n qubits.  Arbitrary quantum circuit abstracted by unitary matrix U  decomposed by 2^n-12^n  unitary two-level matrices, each of which corresponds ^iGamma_jk.  The program produce the Gamma matrix and the coefficients ijk. The quantum  circuit of this decomposition can be visualized as a cascade (from left to right) of this matrix.\n\n\n\n\n\n","category":"method"},{"location":"#GrayCoding.tiffoli_matrix-Tuple{}","page":"Home","title":"GrayCoding.tiffoli_matrix","text":"The unitary matrix U corresponding to of the 3 quibit Tiffoli quantum gate\n\nU=beginpmatrix\n1  0  0  0  0  0  0  0 \n0  1  0  0  0  0  0  0 \n0  0  1  0  0  0  0  0 \n0  0  0  1  0  0  0  0 \n0  0  0  0  1  0  0  0 \n0  0  0  0  0  1  0  0 \n0  0  0  0  0  0  colorred0  colorred1 \n0  0  0  0  0  0  colorred1  colorred0 \nendpmatrix\n\n\n\n\n\n","category":"method"},{"location":"algebra/#Algebraic-framework-of-Gray-Codes","page":"Algebra of Gray Codes","title":"Algebraic framework of Gray Codes","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"The classical algorithmic procedure of encoding and decoding are as follows:","category":"page"},{"location":"algebra/#Encoding","page":"Algebra of Gray Codes","title":"Encoding","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"q-ry digits d to gray digits g conversion.","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"g_i = begincases  d_i   textif  modleft(displaystylesum_j=1^i-1g_j2right)=0   q-1-d_i   textif modleft(displaystyle sum_j=1^i-1g_j2right)=1 endcases","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"and g_1 = d_1.","category":"page"},{"location":"algebra/#Decoding","page":"Algebra of Gray Codes","title":"Decoding","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"d_i = begincases  g_i   textif  modleft(displaystylesum_j=1^i-1g_j2right)=0   q-1-g_i   textif modleft(displaystyle sum_j=1^i-1g_j2right)=1  endcases","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"and d_1 = g_1.","category":"page"},{"location":"algebra/#Linear-Algebraic-Formulation-(N.Rethnakar-2020)","page":"Algebra of Gray Codes","title":"Linear Algebraic Formulation (N.Rethnakar 2020)","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"Example of generator matrix G for binary to gray mapping is given by,","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"G=beginpmatrix 1   colorgray0     colorgray0    colorgray0     colorgray0    colorgray0   \n 1    1   colorgray0   colorgray0   colorgray0   colorgray0   \n colorgray0    1   1  colorgray0   colorgray0   colorgray0   \n colorgray0    colorgray0   1  1  colorgray0   colorgray0   \n colorgray0    colorgray0   colorgray0   1  1  colorgray0   \n colorgray0    colorgray0   colorgray0   colorgray0   1  1 endpmatrix","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"The decoding matrix B=G^-1 is given by,","category":"page"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"  B=beginpmatrix  1    colorgray0     colorgray0     colorgray0     colorgray0     colorgray0   \n1   1   colorgray0   colorgray0     colorgray0    colorgray0   \n1   1   1  colorgray0     colorgray0    colorgray0   \n1   1   1  1    colorgray0    colorgray0   \n1   1   1  1    1   colorgray0   \n1   1   1  1    1   1 endpmatrix","category":"page"},{"location":"algebra/#Generalized-q-ry-Gray-Code","page":"Algebra of Gray Codes","title":"Generalized q-ry Gray Code","text":"","category":"section"},{"location":"algebra/","page":"Algebra of Gray Codes","title":"Algebra of Gray Codes","text":"G=beginpmatrix 1   colorgray0    colorgray0    colorgray0    colorgray0    colorgray0  \n q-1    1   colorgray0   colorgray0   colorgray0   colorgray0  \n colorgray0    q-1   1  colorgray0   colorgray0   colorgray0  \n colorgray0    colorgray0   q-1  1  colorgray0   colorgray0  \n colorgray0    colorgray0   colorgray0   q-1  1  colorgray0  \n colorgray0    colorgray0   colorgray0   colorgray0   q-1  1 endpmatrix equiv beginpmatrix 1   colorgray0    colorgray0    colorgray0    colorgray0    colorgray0  \n -1    1   colorgray0   colorgray0   colorgray0   colorgray0  \n colorgray0    -1   1  colorgray0   colorgray0   colorgray0  \n colorgray0    colorgray0   -1  1  colorgray0   colorgray0  \n colorgray0    colorgray0   colorgray0   -1  1  colorgray0  \n colorgray0    colorgray0   colorgray0   colorgray0   -1  1 endpmatrix_mathbbF_q","category":"page"}]
}
