# Examples

## Recursive construction
![][./assets/gray_recursion2.svg]

TBD
``g=Gb`` and ``b=Bg``, where ``G`` is a Jordan matrix, which is 

```julia-repl
julia> n,q=4,2
julia> GrayCoding.GrayMatrix(n,q)
4×4 Matrix{Int64}:
 1  0  0  0
 1  1  0  0
 1  1  1  0
 1  1  1  1
4×4 Matrix{Int64}:
 1  0  0  0
 1  1  0  0
 0  1  1  0
 0  0  1  1
4×16 Matrix{Int64}:
 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1
 0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1
 0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1
 0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1
4×16 Matrix{Int64}:
 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1
 0  0  0  0  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  1  1  1  1  0  0
 0  1  1  0  0  1  1  0  0  1  1  0  0  1  1  0

julia> G,B,g,b=GrayCoding.GrayMatrix(10,5);
julia> G
10×10 Matrix{Int64}:
 1  0  0  0  0  0  0  0  0  0
 1  1  0  0  0  0  0  0  0  0
 1  1  1  0  0  0  0  0  0  0
 1  1  1  1  0  0  0  0  0  0
 1  1  1  1  1  0  0  0  0  0
 1  1  1  1  1  1  0  0  0  0
 1  1  1  1  1  1  1  0  0  0
 1  1  1  1  1  1  1  1  0  0
 1  1  1  1  1  1  1  1  1  0
 1  1  1  1  1  1  1  1  1  1
 julia>B
 10×10 Matrix{Int64}:
 1  0  0  0  0  0  0  0  0  0
 4  1  0  0  0  0  0  0  0  0
 0  4  1  0  0  0  0  0  0  0
 0  0  4  1  0  0  0  0  0  0
 0  0  0  4  1  0  0  0  0  0
 0  0  0  0  4  1  0  0  0  0
 0  0  0  0  0  4  1  0  0  0
 0  0  0  0  0  0  4  1  0  0
 0  0  0  0  0  0  0  4  1  0
 0  0  0  0  0  0  0  0  4  1
```