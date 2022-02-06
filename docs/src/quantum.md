# Quantum Computing
Broadly stated, the term quantum computation comprieses of the following three elements:
* A register or a set of registers,
* A unitary matrix ``U``, as an abstract representation of the quantum algorithm,
* Measurements to extract the information of interest.

As a mathematial abstraction, a quantum computation is the set ``\{\mathcal{H},U,\{Mm\}\}``, where ``H = \mathbb{C}^{2^{n}}`` is the Hilbert space of an ``n-``qubit register, ``U \in U\left(2^{n}\right)`` represents the quantum algorithm and ``\{M_{m}\}`` is the set of measurement operators. The hardware circuitry along with equipment to control and manipulate the qubits is called a quantum computer.
# Decomposition of Quantum Gates

## Universal Decomposition
Any arbitrary unitary gate acting on n-qubit can be implemented as a cascade of single-qubit and controlled-NOT (CNOT) gates.

```math
\textbf{U}|\psi \rangle = \begin{pmatrix} u_{11} & u_{12} & \ldots & u_{1,N} \\ u_{11} & u_{12} & \ldots & u_{1,N} \end{pmatrix}
```

## CNOT gates
 CNOT stands for Controlled-NOT, which is one of the key quantum logic gate. It is a two qubit gate. The gate flips the second qubit (called the target qubit) when the first gate (the control gate) is ``\lvert 1 \rangle``, while the target qubit remain unchanged when the control gate is in state ``\lvert 0 \rangle``.


```math
\begin{aligned}
 U_{\text{CNOT}} &= \lvert 0 \rangle \otimes \langle 0 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \lvert 1 \rangle \otimes \langle 1 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\\
                 &= \lvert 0 \rangle \langle 0 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix} + \lvert 1 \rangle  \langle 1 \rvert \otimes \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}\\
                 &= \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0\end{pmatrix} 
 \end{aligned}
```

In terms of Gray matrix, we can also simply express this as,
```math
\begin{aligned}
\begin{pmatrix} \lvert \acute{\psi}_{1} \rangle \\ \lvert \acute{\psi}_{2} \rangle \end{pmatrix} &= G_{2} \begin{pmatrix} \lvert \psi_{1} \rangle \\ \lvert \psi_{2} \rangle \end{pmatrix} \\
&= \begin{pmatrix} 1 & 0 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} \lvert \psi_{1} \rangle \\ \lvert \psi_{2} \rangle \end{pmatrix}
\end{aligned}
```

---
The simplest CNOT gate is the single qubit controlled CNOT discussed above, which can be explicitly denoted as ``C^{1}\text{NOT}``. Generalization of this to multi quibit controlled CNOT, denoted by ``C^{n-1}\text{NOT}``.
---

In quantum circuit design, applying a rotation for which the binary represen- tations of i − 1 and j − 1 differ in a single bit can be accomplished by a single fully-controlled one-qubit rotation (a particular Givens rotation) and hence costs a small number of gates. All other rotations require a permutation of data before the rotation is applied and thus should be avoided.