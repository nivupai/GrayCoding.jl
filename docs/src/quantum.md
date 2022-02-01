# Decomposition of Quantum Gates

## Universal Decomposition
Any arbitrary unitary gate acting on n-qubit can be implemented as a cascade of single-qubit and controlled-NOT (CNOT) gates.

```math
\textbf{U}|\psi \rangle = \begin{pmatrix} u_{11} & u_{12} & \ldots & u_{1,N} \\ u_{11} & u_{12} & \ldots & u_{1,N} \end{pmatrix}
```

## CNOT gates
 CNOT stands for Controlled-NOT, which is one of the key quantum logic gate. It is a two qubit gate. The gate flips the second qubit (called the target qubit) when the first gate (the control gate) is ``\lvert 1 \rangle``, while the target qubit remain unchanged when the control gate is in state ``\lvert 0 \rangle``.

 ```math
 U_{\text{CNOT}} &= \lvert 0 \rangle \otimes \langle 0 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \lvert 1 \rangle \otimes \langle 1 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}
 ```
and
 ```math
 \begin{aligned}
                 &= \lvert 0 \rangle \langle 0 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \lvert 1 \rangle  \langle 1 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\\
                 &= \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{pmatrix} 
 \end{aligned}
 ```