#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../GMPinterface.h"

using namespace std;

struct GI;

void printGIvector(const vector<GI>& v);

class GI
    // Gaussian integer
{
    // real and imaginary part, norm
    mpz Re;
    mpz Im;
    mpz N = Re*Re + Im*Im;
public:

    // constructors
    GI() : Re{ 0 }, Im{ 0 } {}
    GI(mpz a) : Re{ a }, Im{ 0 } {}
    GI(mpz a, mpz b) : Re{ a }, Im{ b } {}

    // get member values
    const mpz real() const { return Re; }
    const mpz imag() const { return Im; }
    const mpz norm() const { return N; }

    // set member values
    void set_Re(const mpz& real) { Re = real; N = Re*Re + Im*Im; }
    void set_Im(const mpz& imag) { Im = imag; N = Re*Re + Im*Im; }

    // operator overloads
    GI operator + (const GI& alpha) const { return GI(Re + alpha.Re, Im + alpha.Im); }
    GI operator - (const GI& alpha) const { return GI(Re - alpha.Re, Im - alpha.Im); }
    GI operator * (const GI& alpha) const { return GI(Re * alpha.Re - Im * alpha.Im, Re*alpha.Im + Im * alpha.Re); }
    GI operator / (const GI& alpha) const;      // calculates the quotient from the Euclidean division
    GI operator % (const GI& alpha) const;
    GI& operator += (const GI& alpha) { *this = *this + alpha; return *this; }
    GI& operator -= (const GI& alpha) { *this = *this - alpha; return *this; }
    GI& operator *= (const GI& alpha) { *this = *this * alpha; return *this; }
    GI& operator /= (const GI& alpha) { *this = *this / alpha; return *this; }
    GI& operator %= (const GI& alpha) { *this = *this % alpha; return *this; }
    bool operator < (const GI& alpha) const { return N < alpha.norm(); }     // for sorting by norm
};

// operator functions
GI conj(const GI& alpha);
GI quo(const GI& alpha, const GI& beta);
GI rem(const GI& alpha, const GI& beta);
vector<GI> associates(const GI& alpha);
bool operator == (const GI& alpha, const GI& beta);
bool operator != (const GI& alpha, const GI& beta);
ostream& operator << (ostream& os, const GI& alpha);

// utilities
string to_string(const GI& alpha);
struct GIprimepower { GI alpha; mpz n; };
void printGIproduct(const vector<GI>& v);
GI GIproduct(const vector<GI>& v);
GI GIpow(const GI& alpha, mpz n);
GI GImodularExponentiation(GI alpha, mpz n, GI beta);
GI GIprimary(GI alpha);

// primes and algorithms
bool isGIprime(const GI& alpha);
GI GIgcd(const GI& alpha, const GI& beta);
GI GIgcd(const vector<GI>& GIs);
bool GIpairwisecoprime(const vector<GI>& GIs);
vector<GI> GIexgcd(const GI& alpha, const GI& beta);
GI GIlcm(const GI& alpha, const GI& beta);
GI GIlcm(const vector<GI>& GIs);
GI GImodlinearsolve(const GI& alpha, const GI& beta, const GI& gamma);
GI GIChRem(const vector<GI> alpha, const vector<GI> moduli);
vector<GI> GIprimefactor(GI alpha);
GI QuarticRes(GI alpha, GI pi);
GI QuarticResSymbol(GI alpha, GI beta);

