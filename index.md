***
## Introduction <br>
***

This page was created to exchange functions and computations obtained in $\texttt{Sage}$ for the research regarding the classes $\mathcal{C}$ and $\mathcal{D}$ of bent functions.

The following lines need to be added to your $\texttt{Sage}$ notebook to support cryptographical functions that are available (for more info see [Sage-BooleanFunctions](https://doc.sagemath.org/html/en/reference/cryptography/sage/crypto/boolean_function.html) and [Sage-SBox](https://doc.sagemath.org/html/en/reference/cryptography/sage/crypto/sbox.html)).




```python
from sage.crypto.sbox import SBox
from sage.crypto.boolean_function import BooleanFunction
```

The Maiorana-McFarland class $\mathcal{M}$ is the set of $m$-variable ($m=2n$) Boolean
functions of the form

$$f(x,y)=x \cdot \pi(y)+ g(y), \mbox{ for all } x, y\in\mathbb{F}_2^n,$$

where $\pi$ is a permutation on $\mathbb{F}_2^n$, and $g$ is an arbitrary Boolean function on
$\mathbb{F}_2^n$. From this class Carlet derived the  $\mathcal{C}$ class of bent functions that contains all functions of the form

$$f(x,y) = x\cdot\pi(y)+ \mathbf{1}_{L^{\perp}}(x)$$

where $L$ is any linear subspace of $\mathbb{F}_2^n$, $\mathbf{1}_{L^{\perp}}$ is the indicator function of the space $L^{\perp}$, and
$\pi$ is any permutation on $\mathbb{F}_2^n$ such that:


$(C)$ $\phi(a + L)$ is a flat (affine subspace),  for all $a\in \mathbb{F}_2^n$, where $\phi:=\pi^{-1}$.

The permutation $\phi$ and the subspace $L$  are then said to satisfy the $(C)$ property, or for short _$(\phi,L)$ has property $(C)$_.

Another class introduced by Carlet, called $\mathcal{D}$, is defined similarly as
$$f(x, y)=x\cdot\pi(y) + \mathbf{1}_{E_1}(x) \mathbf{1}_{E_2}(y)$$

where $\pi$ is a permutation on $\mathbb{F}_2^{n}$ and $E_1,E_2$ two linear subspaces of $\mathbb{F}_2^n$ such that $\pi(E_2)=E_1^{\perp}$.

<br>
We will define three $\texttt{Sage}$ functions: 
- the indicator function $\mathbf{1}_A(x)$ which will be denoted with $\texttt{ind(x,A)}$, 
- the trace function $Tr_m^n(x)$ denoted with $\texttt{tr(x,m,n)}$ and 
- the anihilator of a linear space $E$ of $\mathbb{F}_{2^m}$ denoted with $\texttt{dualSpace(E,m)}$.


```python
#Fm=GF(2^m);
def ind(x,L):
    if x in L:
        return Fm.one();
    else:
        return Fm.zero();


def tr(x,m,n):
    return sum([x^(2^(i*m)) for i in [0..(n/m)-1]])


def dualSpace(E,m):
    Fm=GF(2^m);
    t=[];
    for x in Fm:
        b=0;
        for y in E:
            if (x*y).trace()==0:
                b=b+1;
        if b==len(E):
            t.append(x)
    return t;

def shift_space(W,u):
    new=[];
    for v in W:
        new.append(vector(GF(2),v)+vector(GF(2),u));
    return new;

```

**Example:** Let $m=9$ and $s=3$. Then $m/s$ is odd. Furthermore, for $d=284$ we have that $d(2^s+1) \mod (2^9-1)=1$. Let $U=\{1,\alpha,\alpha^2\}$ where $\alpha$ is a primitive element of $\mathbb{F}_{2^3}$ such that $\alpha^3+\alpha+1=0$. 

For $L=\langle U\rangle$ we have that $(\pi^{-1},L)$ satisfies the $(C)$ property. Hence, the function $f:\mathbb{F}_{2^m}\times\mathbb{F}_{2^m}\to\mathbb{F}_2$ defined with $$f(x,y) = x\cdot\pi(y)+ \mathbf{1}_{L^{\perp}}(x)$$ is bent.


```python
Fm=GF(2^9);
sFm=sorted(Fm);

p=Fm.primitive_element();
a=p^((2^9-1)/(2^3-1));
u=sFm[13];

Lperp=[x for x in Fm if ((x).trace()+1)*((x*a).trace()+1)==1];

f=[(x*y^284).trace()+ind(x,Lperp) for x in sFm for y in sFm];
f=BooleanFunction(f);
f.is_bent()
```

**Example:** Let $m=9$ and $s=3$. Then $m/s$ is odd. Furthermore, for $d=284$ we have that $d(2^s+1) \mod (2^9-1)=1$. Let $\alpha$ be a primitive element of $\mathbb{F}_{2^3}$ such that $\alpha^3+\alpha+1=0$. 

For $E_2=\langle \alpha,\alpha^2 \rangle$ we have that $\pi(E_2)=E_2^{\perp}=E_1$. Hence, the function $f:\mathbb{F}_{2^m}\times\mathbb{F}_{2^m}\to\mathbb{F}_2$ defined with $$f(x,y) = x\cdot\pi(y)+ \mathbf{1}_{E_1}(x)\mathbf{1}_{E_2}(y)$$ is bent.


```python
Fm=GF(2^9);
sFm=sorted(Fm);

p=Fm.primitive_element();
a=p^((2^9-1)/(2^3-1));

E2=[Fm.zero(),a^2,a+a^2,a];
E1=dualSpace(E2,9);

f=[(x*y^284).trace()+ind(x,E1)*ind(y,E2) for x in sFm for y in sFm];
f=BooleanFunction(f);
f.is_bent()
```

Let $\pi$ be a permutation on $\mathbb{F}_{2^m}$, $L\subset\mathbb{F}_{2^m}$ be a linear subspace of $\mathbb{F}_{2^m}$ such that $(\pi^{-1},L)$ satisfies the $(C)$ property, and let $E_1,E_2\neq \{0\}$ be two linear subspaces of $\mathbb{F}_{2^m}$ such that $\pi(E_2)=E_1^{\perp}$. Then the class of bent  functions $f:\mathbb{F}_{2^m}\times\mathbb{F}_{2^m}\to\mathbb{F}_{2}$ containing all functions of the form  
$$f(x,y)=Tr_1^m(x\pi(y))+a_0\mathbf{1}_{L^{\perp}}(x)+a_1\mathbf{1}_{E_1}(x)\mathbf{1}_{E_2}(y), \,\;  a_i \in \mathbb{F}_{2}$$
is called $\mathcal{CD}$ and is a superclass of $\mathcal{C}$ and $\mathcal{D}$.

**Example:** Let $m=9$ and $s=3$. Then $m/s$ is odd. Furthermore, for $d=284$ we have that $d(2^s+1) \mod (2^9-1)=1$.  Let $U=\{1,\alpha,\alpha^2\}$ where $\alpha$ is a primitive element of $\mathbb{F}_{2^3}$ such that $\alpha^3+\alpha+1=0$. 

For $L=\langle U\rangle$ we have that $(\pi^{-1},L)$ satisfies the $(C)$ property. For $E_2=\langle \alpha,\alpha^2 \rangle$ we have that $\pi(E_2)=E_2=E_1^{\perp}$. Hence, the function $f:\mathbb{F}_{2^m}\times\mathbb{F}_{2^m}\to\mathbb{F}_2$ defined with $$f(x,y) = x\cdot\pi(y)+ \mathbf{1}_{L^{\perp}}(x)+\mathbf{1}_{E_1}(x)\mathbf{1}_{E_2}(y)$$ is bent.


