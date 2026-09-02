---
jupytext:
    formats: md:myst
    text_representation:
        extension: .md
        format_name: myst
kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Basis transformations

## Bases

Until now, we have been working with one specific coordinate system in
which vectors are identified with their coordinate vectors with
respect to the *standard basis vectors*

$$
\ve_1=\begin{pmatrix}1\\ 0\\ 0\\ \vdots\\ 0\end{pmatrix},\quad
\ve_2=\begin{pmatrix}0\\ 1\\ 0\\ \vdots\\ 0\end{pmatrix},\quad
\ldots,\quad
\ve_n=\begin{pmatrix}0\\ 0\\ 0\\ \vdots\\ 1\end{pmatrix}.
$$

In other words, a vector $\vv=\begin{pmatrix}v_1\\ v_2\\ \vdots\\
v_n\end{pmatrix}$ means the same as the corresponding linear
combination of the standard basis vectors:

$$
\vv = v_1 \ve_1 + \cdots + v_n \ve_n.
$$

In many situations, it is useful to express vectors in a different
coordinate system.  For this we need the general concept of a *basis*.

:::{prf:definition} Basis
A *basis* of $\RR^n$ consists of $n$ vectors $(\vb_1,\ldots,\vb_n)$
such that every vector $\vv\in\RR^n$ can be written in exactly one way
as a *linear combination*

$$
\vv=c_1\vb_1+\cdots+c_n\vb_n
$$

of $\vb_1,\ldots,\vb_n$ with scalars $c_1,\ldots,c_n$.  These scalars
are called the *coordinates* of $\vv$ relative to
$(\vb_1,\ldots,\vb_n)$.  The vector $\begin{pmatrix}c_1\\ \vdots\\
c_n\end{pmatrix}$ is called the *coordinate vector* of $\vv$ relative
to $(\vb_1,\ldots,\vb_n)$.

:::

Alternatively, a basis consists of $n$ vectors $\vb_1,\ldots,\vb_n$
with the following two properties:

- they *span* the space $\RR^n$, i.e. every vector is a linear
  combination of $\vb_1,\ldots,\vb_n$;

- they are *linearly independent*, i.e. the only way to write the zero
  vector as a linear combination of $\vb_1,\ldots,\vb_n$ is as
  $0\vb_1+\cdots+0\vb_n$.

:::{note}
A priori, we could have defined a basis as a collection of *some*
number (say $m$) of vectors that span the space and are linearly
independent.  However, it can be shown that in fact any two bases contain
the same number of vectors.  Since the standard basis of $\RR^n$
consists of $n$ vectors, the same therefore holds for *any* basis.
:::

## Basis transformation matrices

Suppose we have two bases of $\RR^n$, say

$$
B=(\vb_1, \ldots, \vb_n)\quad\text{and}\quad
C=(\vc_1, \ldots, \vc_n).
$$

:::{prf:definition} Basis transformation matrix

The *basis transformation matrix* (or *change-of-basis matrix*) from
$B$ to $C$ is the matrix $P=P_{B\to C}$ such that the $j$-th column of
$P$ contains the coordinates of $\vb_j$ relative to the basis
$(\vc_1,\ldots,\vc_n)$.
:::

The definition of $P$ means (note the difference with the formula
{eq}`matrix-vector` for matrix-vector multiplication)

$$
\vb_j = \sum_{i=1}^n P_{i,j}\vc_i.
$$

Given a vector $\vx$, we write
$\vx^{(B)}=(x_1^{(B)},\ldots,x_n^{(B)})$ for the coordinate vector of
$\vx$ with respect to the basis $B=(\vb_1,\ldots,\vb_n)$, and we
define $\vx^{(C)}$ similarly.  Then we calculate

$$
\begin{aligned}
\vx &= \sum_{j=1}^n x_j^{(B)} \vb_j\\
&= \sum_{j=1}^n x_j^{(B)}\sum_{i=1}^n P_{i,j}\vc_i\\
&= \sum_{i=1}^n\left(\sum_{j=1}^n P_{i,j} x_j^{(B)}\right) \vc_i\\
&= \sum_{i=1}^n\bigl(P \vx^{(B)}\bigr)_i \vc_i,
\end{aligned}
$$

which shows that

$$
\vx^{(C)} = P_{B\to C}\vx^{(B)}.
$$

Similarly, given a linear map $A$, we write $A^{(B)}$ for the matrix
of $A$ with respect to the basis $B=(\vb_1,\ldots,\vb_n)$, and
likewise for $A^{(C)}$.  Then $A^{(B)}$ and $A^{(C)}$ are related by

$$
A^{(C)} = PA^{(B)}P^{-1}.
$$

Namely, consider what happens when we multiply this matrix to a
coordinate vector of the form $\vx^{\vw}$.  First, applying $P^{-1}$
to $\vx^{\vw}$ gives $\vx^{\vv}$.  Next, applying $A^{(B)}$ gives
$(A\vx)^{(B)}$.  Finally, applying $P$ gives $(A\vx)^{(C)}$, which
is what we wanted to show.

```{code-cell} ipython3
:tags: [remove-input]
from jupyterquiz import display_quiz
display_quiz("_static/quizzes/basis_trans.json")
```

## Invariance under basis transformations

In [](matrix-operations) we attached several quantities to a square
matrix: the trace, the determinant and the characteristic polynomial.
In {ref}`trace-det-basis-invariant`, you will check that the trace and
the determinant are invariant under basis transformation: if $A$ and
$P$ are square matrices of the same size with $P$ invertible, then we
have

$$
\tr(PAP^{-1})=\tr A
\quad\text{and}\quad
\det(PAP^{-1})=\det A.
$$

Because the characteristic polynomial is also defined as a
determinant, it too is invariant under basis transformation.

***

## Exercises

:::{exercise}
Which of the following collections of vectors are bases for $\RR^2$?

1. $\begin{pmatrix}1\\ 2\end{pmatrix}$, $\begin{pmatrix}3\\ 4\end{pmatrix}$,
   $\begin{pmatrix}5\\ 6\end{pmatrix}$

2. $\begin{pmatrix}0\\ 0\end{pmatrix}$, $\begin{pmatrix}1\\ 1\end{pmatrix}$

3. $\begin{pmatrix}1\\ 2\end{pmatrix}$, $\begin{pmatrix}-3\\ 4\end{pmatrix}$
:::

:::{exercise}
Given the two bases

$$
B=(\vb_1,\vb_2)=\biggl(\begin{pmatrix}1\\ 2\end{pmatrix},
\begin{pmatrix}2\\ 1\end{pmatrix}\biggr)
\quad\text{and}\quad
C=(\vc_1,\vc_2)=\biggl(\begin{pmatrix}0\\ 1\end{pmatrix},
\begin{pmatrix}2\\ 1\end{pmatrix}\biggr)
$$

determine the basis transformation matrix $P_{B\to C}$.
:::

:::{exercise}
:label: trace-det-basis-invariant

Suppose $A$ and $B$ are two $n\times n$-matrices and $P$ is an
invertible $n\times n$-matrix such that $B=PAP^{-1}$.  Show that

$$
\tr B = \tr A\quad\text{and}\quad \det B = \det A.
$$
:::

:::{exercise}
Consider two $n\times n$-matrices $A$ and $B$.

1.  Suppose there exists an invertible $n\times n$-matrix $P$
    satisfying $B=PAP^{-1}$.  Show that $A$ and $B$ have the same
    characteristic polynomial.

2.  Give an example where $A$ and $B$ have the same characteristic
    polynomial but a matrix $P$ as above does not exist.
:::
