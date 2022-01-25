# Algebraic framework of Gray Codes


The classical algorithmic procedure of encoding and decoding are as follows:
## Encoding
q-ry digits ``d`` to gray digits ``g`` conversion.

``
g_{i} = \begin{cases} 
d_{i} , & \text{if}  \mod\left(\displaystyle{\sum_{j=1}^{i-1}{g_{j}}},2\right)=0 \\ q-1-d_{i} , & \text{if} \mod\left(\displaystyle \sum_{j=1}^{i-1}{g_{j}},2\right)=1 \\
\end{cases}
``
and ``g_{1} = d_{1}``.
## Decoding

``
d_{i} = \begin{cases} 
g_{i} , & \text{if}  \mod\left(\displaystyle{\sum_{j=1}^{i-1}{g_{j}}},2\right)=0 \\ 
q-1-g_{i} , & \text{if} \mod\left(\displaystyle \sum_{j=1}^{i-1}{g_{j}},2\right)=1 \\
\end{cases}
``
and ``d_{1} = g_{1}``.

TBD