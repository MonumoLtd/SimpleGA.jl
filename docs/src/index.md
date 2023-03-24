```@meta
CurrentModule = GADraft
```

# GeometricAlgebra.jl

_A Julia package implementing Geometric Algebra_. 

GeometricAlgebra provides implementations of various low-dimensional geometric algebras that are designed to be lightweight and fast. Each has its own unique representation and namespace, and it is possible to have multiple geometric algebras insantiated simultaneously. This allows the programmer to operate in whichever algebra is most efficient for a given task.

Two further algebras are also supported, GA(4,4) and GA(32,32). These are less optimised, and are provided more for investigation.

All algebras are defined over the reals and follow Julia's internal promotion rules. 


## Installation
```
pkg> add GeometricAlgebra
```

## First Example
```
using GeometricAlgebra
newbasis = basis(2,0) #Creates a 2D basis
a = e1 + e2
b = 2.0*e1 + 3*e2 
a*b #Evaluates the geometric product of a and b.
```

# API

## Bases and the Even / Odd trick
This library employs an optimisation achieved by separating even and odd multivectors into separate types. This reduces the size of the representation and halves the time to compute a geometric product. An algebra-specific map is used to convert odd elements to even, so that both share the same representation. 

But there is an important consequence: **you cannot combine even and odd elements in any of the optimised algebras**. 

In practice it is highly unusual that you should want to combine even and odd elements, but there are some occasions where it does make sense. For these you should work in an algebra one dimension higher that contains the desired algebra in its even subalgebra. The standard example is special relativity. You might want to add a scalar to a vector in G(3,0) to simulate a 4-vector, but the implementation here will not allow it. So instead you should work one dimension higher, in the spacetime algebra G(1,3), which is the appropriate setting for relativistic calculations. 

The standard interface to an algebra is the ```basis``` function. This can be called at any point in code to call into existence a copy of the basis vectors for the desired algebra. The coefficients in these basis vectors are all ```Int8```, which get promoted to whatever the proagrmmer desires by multiplying by reals. 

There are various interfaces to ```basis```
```
basis(p) # Returns a basis for G(p,0), if it exists. If not it returns the nearest algebra that contains G(p,0).
basis(p,q) # Returns a basis for G(p,q) if it exists. If not it returns the nearest basis that contains G(p,q).
basis(p,q,r) # Returns a basis for G(p,q,r), where r generators are null if the algebra is implemented. If not it returns the nearest basis that contains G(p,q,r)
basis("name") # Returns the basis of the named algebra.
```

As well as accessing the basis via the basis function, each basis can be called directly by name, so for example ```bas20``` is the name of the basis for G(2,0).

The following algebras are currently implemented, together with their name for ```basis("name")``` and the name of the basis:
* G(2,0), name: "GA20", basis: ```bas20```
* G(3,0), name: "GA30", basis: ```bas30```
* G(4,0), name: "GA40", basis: ```bas40```
* G(1,3), name: "STA", basis: ```basSTA```
* G(3,0,1), name: "PGA", basis: ```basPGA```
* G(4,1), name: "CGA", basis: ```basCGA```
* G(2,4), name: "GA24", basis: ```bas24```
* G(3,3), name: "GA33, basis: ```bas33```

Note that the constructors for each type are not exposed directly. They are all of the form ```GeometricAlgebra.GA20.Even()```, for example. You should fine working directly with the basis vectors is the simplest way to create multivector objects.

In addition, two further algebras are provided with a separate implementation. For these there is no distinction between even and odd elements; everything is a single ```Multivector``` type. Some optimisations have been applied to these algebras, but they are not nearly as fast as the ones listed above.

* G(4,4), name: "GA44", basis: ```bas44```
* G(32,32), name: "GA64", basis: ```bas64```

Each of these algebras is described in more detail below, including listing the additional functions that are exposed that are algebra specific.



## Arithmetic Functions

Most of the functionality is provided by overloading the Base Julia functions, so +,-,* all behave as expected, as does division by a real. The library also builds on LinearAlebra and overloads two functions from there:

```
tr(A) # The scalar part of A
dot(A,B) # The scalar part of the product AB, usually written <AB>.
A' = adjoint(A) # The reverse of A. 
```
Both of these functions return Reals, so take us out of whatever type is being employed for the arguments. These are key functions for extracting values at the end of geometric algebra calculation.

### Projection
Projection onto a single grade is performed with the project operator

```
Project(A,n) # Returns the grade n component of A.
```
This operation returns a multivector object. 


### Exponentiation
An exponential operator is included for Even multivectors. Two version exist

```
Base.exp(A) # Exponentiates the full even multivector A.
bivector_exp(A) # Exponentiates the bivector part of A.
```
Clearly if the argument is a bivector both ```exp(A)``` and and ```bivector_exp(A)``` will written the same answer (to some precision). The reason for having a separate ```bivector_exp()``` command is that it can take advantage of various optimisations to ensure rotors are calculated as efficiently as possible. A separate command was thought to be more efficient than implementing a runtime check on the argumed of ```exp(A)``` to determine if the argument was a bivector.

### Further Functions
```dot``` is a simple utility function for doting scalars into a basis set. The code is simple:
```
dot(xs,ys) = reduce(+,map((x,y)->x*y,xs,ys))
```

One advantage of how this library uses Julia's namespace and module structure is that it is straightforward to have multiple algbras in flight at one time. In a given application you may want to move elements between algebras, which can be achieved with the ```dot``` function mainly. There were too many cases to try and provide a single conversion function, and too many different ways of converting (projective split, conformal split ...). One function that is provided is a universal ```embed``` function that lifts elements into the G(4,4) algebra. This is primarily for testing purposes, but might form the basis for writing more specific applications.
```
embed(A) # Takes a multivector from one of the smaller algebras and embeds it in G(4,4)
```

# The Core Algebras

## G(2,0)
G(2,0) is the algebra of the Euclidean plane. Even elements are stored as complex numbers. Odd elements (vectors) are mapped to even elements via
```
a -> e[1]*a
```

As the even elements are complex numbers, some additional operations are permitted for these:
```
A/B # Divides the even multivector A by B. Complex numbers are a division algebra, and the order is unambiguous so this is fine to allow in the this algebra.
Base.log(), Base.Real and Base.Imag all act as expected.
```





Documentation for [GADraft](https://github.com/tpgillam/GADraft.jl).

```@index
```

```@autodocs
Modules = [GADraft]
```
