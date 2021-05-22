## GaussianIntegers

This repository is for all the code, I have written concerning the Gaussian integers. The code contains algorithms for the following tasks:

* Simple arithmetic including +, -, * and Euclidean division with the Gaussian integers
* Taking conjugates, determining associates, computing powers, modular exponentiation
* The Euclidean algorithm and the extended Euclidean algorithm
* Solving modular linear equations and the Chinese remainder theorem
* Prime factorization
* Computing the quartic/biquadratic residue symbol

# Prerequisites and setup

Compiling the code requires GMP (The GNU Multiple Precision Arithmetic Library) installed properly. For information on download and setup, consult their [website](https://gmplib.org/).

The files "GMPinterface.h" and "GMPinterface.cpp" handle all the interactions with the main files "gaussian_integers.h" and "gaussian_integers.cpp" with the GMP library. Note that the multiple precision integer class `mpz_class` is simply abbreviated to `mpz` and similarly for `mpf_class`. The GMP interface files also contain some basic number theoretic algorithms including the Jacobi symbol, a simple factorization algorithm, Euler's phi function and an algorithm for computing modular square roots (modulo a prime). 

# The `GI` class

The main object is the `GI` class which has three constructors. `GI()` is just 0, `GI(mpz a)` is just a (a Gaussian integer with imaginary part 0) and `GI(mpz a, mpz b)` is a + bi. Ordinary arithmetic and writing to the console is supported. For example,

```c++
GI alpha {1, 2};
GI beta {-3, 5};
cout << alpha + beta << endl;
cout << alpha - beta << endl;
cout << alpha * beta << endl;
```
will output
```c++
-2+7i
4-3i
-13-i
```

The division operator / is also supported. This applies Euclidean division and outputs the quotient, so the output of
```c++
cout << beta / alpha;
```
becomes
```c++
1+2i
```
because -3 + 5i = (1 + 2i)(1 + 2i) + i. The modulo operator % also works as intended:
```c++
cout << beta % alpha; // outputs i
```

# Some examples of usage

Here are further examples of using the library. What is the greatest common divisor of, say, 347 + 89i and 117 - 547i?

```c++
cout << GIgcd({347, 89}, {117, -547}); // outputs 1 + i
```
The extended Euclidean algorithm takes two Gaussian integers and outputs a vector with the greatest common divisor and the two coefficients in Bézout's lemma. For example

```c++
vector<GI> v = GIexgcd({347, 89}, {117, -547});
for (GI a : v)
  cout << a << endl;
```
outputs:

```c++
1+i
164-41i
-22-106i
```
and we can verify that (164 - 41i)(347 + 89i) + (-22 - 106i)(117 - 547i) = 1 + i. Let us now solve a modular linear equation. For example, let us solve (347 + 89i)x = 117 - 547i mod 2003:

```c++
cout << GImodlinearsolve({347, 89}, {117, -547}, {2003}); // outputs 522+436i
```
And we can verify that (347 + 89i)(522 + 436i) = 117 - 547i modulo 2003. On the other hand, there is no solution to (347 + 89i)x = 117 - 547i mod 1 + 2i since

```c++
GImodlinearsolve({347, 89}, {117, -547}, {2003}); // outputs "No solution to equation (347+89i)rho = 117-547i (mod 1+2i)"
```
The function `GIprimefactor` factors a Gaussian integer into a product og primes and a single unit (1, -1, i or -i). The output is a vector of these Gaussian integers. To write the factorization in a nice way as a product, use the function `printGIproduct`. As an example,

```c++
printGIproduct(GIprimefactor({1284, -418764}));
```
will output

```c++
(-i)(1+i)^5(3)(-7+2i)(-32+37i)(-65-24i)
```
as the factorization of 1284 - 418764i. As a final example, we may compute the quartic residue symbol (alpha/beta)_4, where beta is not divisible by 1 + i (so beta must be an odd Gaussian integer). For example

```c++
cout << QuarticResSymbol({347, 89}, {5, 6}); // outputs -1
```
This means that 347 + 89i is not a quartic residue modulo the prime 5 + 6i (use the function `isGIprime` with input {5, 6} to verify this). It is however a quadratic residue modulo 5 + 6i.
