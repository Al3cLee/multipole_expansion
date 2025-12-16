#import "preamble.typ": *

= Multipole Expansion

We briefly review multipole expansion before diving into the code.

#notation[
The Einstein summation convention is used, except for the Taylor expansion.
Tensors formed by juxtaposing vector components $x^(i), x^(j), x^(k)$, etc. are sometimes written as
the vector's name with many indices like $x^(i j k)$, as is the standard shorthand.
]

#definition("Green's function for static charge")[
The Green's function $phi_(a)(va(x)) := e_(a) slash abs(va(x) - va(x)_(a))$ is the solution to
the point-charge Poisson equation

$ laplacian phi_(a) = delta^((3)) (va(x) - va(x)_(a)). $
]<dfn:static_greens_function>

#motivation[
If $abs(va(x)_(a)) << abs(va(x))$, physically we are observing at $va(x)$ which
is "far away from the source", and hence can expand $phi(va(x))$ with respect to $va(x)_(a)$:

$ phi_(a)(va(x)) = 
  sum_(k=0)^(oo) 1/k! x_(a)^(i_(1)) ... x_(a)^(i_(k))
  pdv(,x_(a)^(i_(1)),...,x_(a)^(i_(k)),total:k)
  e_(a)/abs(va(x)-va(x)_(a)) quad "at" va(x)_(a)=0 $ <eqn:expansion_phi>
]
Using substitution of variables

$ pdv(,x_(a)^(i)) = (-1) times pdv(,(x^(i)-x_(a)^(i))) $

we can write the derivative part in @eqn:expansion_phi as

$ partial_(x_(a)^(i_(q))) ... partial_(x_(a)^(i_(k))) 1/abs(va(x)-va(x)_(a))
=
(-1)^(k) partial_(y^(i_(1))) ... partial_(y^(i_(k))) 1/abs(va(y)), $

where the derivatives are taken at $va(x)_(a)=0$ and $va(y)=va(x)$, respectively.

== Derivative Tensors

#definition("Derivative tensor")[
We define the derivative tensor of order $k$ as

$ D_(i_(1),...,i_(k)) (va(x)) := (-1)^(k) partial_(y^(i_(1))) ... partial_(y^(i_(k))) 1/abs(va(y)) quad "at" va(y)=va(x). $
]<dfn:derivative_tensor>

Because partial derivatives commute, whenever a pair of indices $i_(p),i_(q)$ are contracted (i.e. traced over),
we can commute them to the rightmost positions and use

$ delta^(p,q) partial_(i_(p)) partial_(i_(q)) =: laplacian. $

#result[
The derivative tensor for the Green's function $phi_(a)(va(x))$ is traceless:
$ forall (p,q), tr_((p,q)) D := delta^(i_(p),i_(q))D =0. $
]

#remark[
This traceless property comes from $laplacian (1 slash r) =0$.
Each source moment $x_(a)^(i_(1))...x_(a)^(i_(k))$ has coefficient $D_(i_(1),...,i_(k))$ means:
although we are expanding with respect to $va(x)_(a)$,
the expansion is carried out around something rotation-invariant: the point charge at the origin.
]

#motivation[
Let us now combine @dfn:derivative_tensor and @dfn:static_greens_function to see how the
expansion of $phi$ can be expressed conveniently with $D$.
For this purpose, we need to estimate the order of $D$ in $abs(va(x))$.
]
Inside $D$, the term with the lowest power of $abs(va(x))$ is the one where
every partial derivative operator acted on the $1 slash abs(va(x))^(...)$ term.
From direct calculation, this term is
$ & (-1)^(k) (-1/2 -0)(-1/2 -1)...(-1/2 - (k-1) )  times 2^(k) times y^(i_(1))... y^(i_(k)) times (r^(2))^(-1 slash 2 -k) \
  =& (2k-1)!! times ( y^(i_(1))...y^(i_(k)) ) / r^( 2k +1 ) quad "with" r := abs(va(x)), va(y)=va(x) . $
Starting from @eqn:expansion_phi, we insert $D$ and manually multiply it with $r^(2k + 1)$
such that its explicit $r$ dependence is in non-negative powers. We have

$ phi_(a)(va(x))
  &= e_(a) sum_(k=0)^(oo) 1/k! x_(a)^(i_(1)...i_(k)) D_(i_(1)...i_(k))\
  &= e_(a) sum_(k=0)^(oo)1/k! x_(a)^(i_(1)...i_(k))/r^(k) (D_(i_(1)...i_(k)) times r^(2k+1))/r^(k+1), $

with $x_(a)^(i_(1),...,i_(k)) := x_(a)^(i_(1)) ... x_(a)^(i_(k))$ is
the source moment tensor.

#definition("Scaled derivative tensor")[
The $k$th-order scaled derivative tensor in the Cartesian basis is defined as
$ D'_(i_(1),...,i_(k)) := D_(i_(1),...,i_(k)) times r^(2k+1) $
where the factor $r^(2k+1)$ makes sure $D'$ has non-negative power in $r$.
]

#result[
Using the scaled derivative tensor $Q$, the expansion of $phi_(a)$ can be expressed as
$ phi_(a)(va(x))
  =
  e_(a) sum_(k=0)^(oo) 1/(r^(k+1) k!) n_(a)^(i_(1),...,i_(k))(va(x)) D'_(i_(1),...,i_(k)), $<eqn:phi_as_Dprime>
where $n_(a) := x_(a) slash r^(k)$ is the _normalized_ source moment tensor.
]

#motivation[
We have worked out an expansion, but can we see
from @eqn:phi_as_Dprime that the low-order terms really dominate the series?
Physically this makes sense: seen from very far apart, any charge distribution
becomes a point charge. Let us verify this intuition.
]

The object $n_(a)^(i_(1)...i_(k))$ requires no attention, since
it is already normalized by $r$, which we assume to be much larger than $r_(a):=abs(va(x)_(a))$.
With increasing $k$ it will only become smaller.

The tensor $D'$ has a leading term with no $r$ dependence:
$ (2k-1)!! times x^(i_(1)) ... x^(i_(k)). $
This term is manifestly not traceless: tracing over any pair of indices will produce
a prefactor $r^(2)$.
The other terms in $D'$, therefore, are responsible for
subtracting away such traces. For $k=2m$, the highest-power trace term is
produced by taking $m$ traces over $m$ pairs of indices, and is of power $2m = k$.
Therefore,
$ D' = "leading symmetric term" - ( "trace terms" ) $
where the leading symmetric term scales as $r^(k)$ and
any trace term scaling as no more than $r^(k)$.
This warrants no mathematical proof, but this mental model helps us understand that
the scaling of $D'$ is contained by $r^(k)$.

#result[
$D' ~ r^(k)$ is suppressed by the prefactor $1 slash r^(k+1)$,
and since the tensor $n_(a)$ decreases with larger $k$, @eqn:phi_as_Dprime is indeed dominated by low-order terms.
]

== Multipole Moment Tensor

#motivation[
In @eqn:phi_as_Dprime, the dependence of $phi$ on $va(x)$ is not very clear,
because $D'$ is rather involved. We want an expression of $phi$ whose
dependence on $va(x)$ is clear, and upon fixing $va(x)$, would give a
physical quantity that depends on $va(x)_(a)$.

More specifically,
we want to move the complexity of $D'$ to $n_(a)$.
A tempting approach is to replace the tensor
$ n_(a)^(i_(1),...,i_(k)) := x_(a)^(i_(1),...,i_(k))/r^(k) $
with the actual normalized position tensor
$ n^(i_(1),...,i_(k)) := x^(i_(1),...,i_(k))/r^(k), $
and at the same time, replace $D'$ with a tensor that depends on $va(x)_(a)$.
]

The tensor $D'_(i_(1), dots.c, i_(k))$ has leading term proportional to $x^(i_(1) dots.c i_(k))$,
so for this term we can indeed swap the $k$-th tensor moment $x_(a)$ with $x$.
Its sub-leading terms are more complicated.
The term with $m$ Kronecker deltas is
$ (-1)^(k) (-1)^(k-m) (2k-2m-1)!! times
 r^(2m) delta_(i_(1),j_(1))dots.c delta_(i_(m),j_(m)) 
 times y^(ell_(1)) dots.c y^(ell_(k-2m)), $<eqn:sub-leading-in-Dprime>
where the leading $(-1)^(k)$ comes from the definition of $D'$.
Each $delta_(i,j)$ is due to $partial_(j)$ acting on the denominator,
followed by $partial_(i)$ that acts on the numerator,
so each $delta$ means _one less derivative operator acting on the denominator_.
This manifests in the coefficient $(-1)^(k-m) (2k-2m-1)!!$.
The indices $(ell_(1), dots.c ell_(k-2m) )$
are the ones that survived the $m$ contractions.
This sub-leading term is contracted with
$ n_(a) := x_(a)^(i_(1),dots.c,i_(k))/r^(k) $
to form the contribution to $phi_(a)$.

Recall that we want to extract
$y^(i_(1) dots.c i_(k))$
from this term in $D'$.
We assign the delta symbols and prefactors in @eqn:sub-leading-in-Dprime
to $n_(a)$ and focus on the remaining part
$ r^(2m) y^(ell_(1)) dots.c y^(ell_(k-2m)) $
Our task of constructing $x^(i_(1), dots.c, i_(k))$ from this term
is faced with a difficulty: @eqn:sub-leading-in-Dprime has
$k-2m<k$ indices on $y$.
We overcome this by manually restoring $2m$ indices, using
$ y^(ell_(1),dots.c,ell_(k-2m)) r^(2m) = y^(ell_(1),dots.c,ell_(k)) times
  delta_(i_(1),j_(1))times dots.c times delta_(i_(m),j_(m)). $<eqn:index-restore>
Now, these $y$ components can be combined with $1/r^(k)$ to form $n$,
and we can collect the remaining terms.

From the perspective of $n_(a)$,
it received the prefactors and $delta$ symbols, which
would cause it to contract $m$ pairs and produce $r_(a)^(2m)$
with $r_(a):=abs(va(x)_(a))$.
Then it received some $delta$ symbols again,
this time generated from the index-restore trick @eqn:index-restore.
Finally, it donated its $1 slash r^(k)$ so that
the extracted $k$-th order $y$ tensor becomes $n$.

In all, the $n_(a)$ tensor has become
$ (-1)^(k) (-1)^(k-m) (2k-2m-1)!! times
 r_(a)^(2m) delta_(i_(1),j_(1))dots.c delta_(i_(m),j_(m)) 
 times x_(a)^(ell_(1)) dots.c x_(a)^(ell_(k-2m)), $
which is _exactly the same_ as the general term in $D'$,
except for the new subscripts $a$.
We have completed the task to "move"
the complexity from $D'$ to another object: the swap
$ D' mapsto n,\ n_(a) mapsto "the above object"  $
leaves $phi_(a)(va(x))$ unchanged.

#definition("Multipole Tensor")[
The $k$-th order multipole tensor $Q$ is defined as
$ Q := sum_(m=0)^(floor.l k slash 2 floor.r) (-1)^(m) (2k-2m-1)!! times
 r_(a)^(2m) delta_(i_(1),j_(1))dots.c delta_(i_(m),j_(m))
 times x_(a)^(ell_(1)) dots.c x_(a)^(ell_(k-2m)), $
where the indices $(ell_(1),dots.c,ell_(k-2m))$
are the ones except those inside $delta$ symbols,
and
$ floor.l k slash 2 floor.r equiv min{m | m in ZZ, m >= k slash 2} $
is the lower-rounded integer of $k slash 2$.
]

#result[
Using the multipole moment tensor $Q$, the expansion of $phi_(a)$ can be expressed as
$ phi_(a)(va(x))
  =
  e_(a) sum_(k=0)^(oo) 1/(r^(k+1) k!) n^(i_(1),...,i_(k))(va(x)) Q_(i_(1),...,i_(k)), $<eqn:phi-as-Q>
where $n := x slash r^(k)$ is the _normalized_ field (as in "far-field") moment tensor.
]

#remark[
This expression is physically very clear. At order $k$,
for a given spatial moment component $n$,
the tensor $Q$ completely specifies the contribution
to $phi_(a)$ by this moment $n$.
]

#remark[
The tensors $D, D'$ and $Q$ are similar: they are all traceless, and
they share the same algebraic structure.
Most importantly, they are all contracted with a totally symmetric
position moment to give the contribution to $phi_(a)$.
]
