#include "gaussian_integers.h"

void printGIvector(const vector<GI>& v)
// prints a vector of GIs in '{}'
{
    cout << "{";
    for (int i = 0; i < v.size(); ++i) {
        if (i != v.size() - 1)
            cout << v[i] << ", ";
        else
            cout << v[i];
    }
    cout << "}";
}

// operator and member functions
// -----------------------------------------------------------------------------------

GI conj(const GI& alpha)
{
    return GI(alpha.real(), -alpha.imag());
}

GI quo(const GI& alpha, const GI& beta)
{
    GI prod = alpha * conj(beta);
    mpf p = mpf(prod.real())/beta.norm();
    mpf q = mpf(prod.imag())/beta.norm();
    GI gamma {round(p), round(q)};
    return gamma;
}

GI rem(const GI& alpha, const GI& beta)
{
    return alpha % beta;
}

GI GI::operator / (const GI& alpha) const
{
    return quo(*this, alpha);
}

GI GI::operator % (const GI& alpha) const
{
    GI gamma = quo(*this, alpha);       // determine the quotient gamma
    GI rho = *this - alpha*gamma;       // calculate remainder
    return rho;
}

vector<GI> associates(const GI& alpha)
{
    vector<GI> assoc = {alpha, GI(-1,0)*alpha, GI(0,1)*alpha, GI(0,-1)*alpha};
    return assoc;
}

bool operator == (const GI& alpha, const GI& beta)
{
    return alpha.real() == beta.real() && alpha.imag() == beta.imag();
}

bool operator != (const GI& alpha, const GI& beta)
{
    return !(alpha == beta);
}

ostream& operator << (ostream& os, const GI& alpha)
{
    string s = to_string(alpha);
    return os << s;
}

// utilities
// -----------------------------------------------------------------------------------

string to_string(const GI& alpha)
// converts a Gaussian integer to a string of the form Re + Im*i
{
    string s;
    if (alpha.real() != 0) {
        s += alpha.real().get_str();    // we add the real part to the string
        if (alpha.imag() != 0) {        // three cases for the imaginary part
            if (alpha.imag() == 1)
                s += "+i";
            else if (alpha.imag() == -1)
                s += "-i";
            else if (alpha.imag() < 0) {
                mpz temp = abs(alpha.imag());
                s += "-" + temp.get_str() + "i";
            }
            else
                s += "+" + alpha.imag().get_str() + "i";
        }
    }
    else {
        if (alpha.imag() == 0)
            s = "0";
        if (alpha.imag() == 1)
            s = "i";
        if (alpha.imag() == -1)
            s = "-i";
        if (alpha.imag() != 1 && alpha.imag() != -1 && alpha.imag() != 0)
            s = alpha.imag().get_str() + "i";
    }
    return s;
}

void printGIproduct(const vector<GI>& v)
// write out a product of GIs with parentheses and powers
{
    if (v.size() == 0) {
        cout << "Error, empty vector" << endl;
        return;
    }
    vector<GI> temp = v;
    sort(temp.begin(), temp.end());

    int n = 1;
    cout << "(" << temp[0] << ")";
    for (int i = 1; i < temp.size(); ++i) {
        if (temp[i] != temp[i - 1]) {
            if (n > 1)
                cout << "^" << n;
            n = 1;
            cout << "(" << temp[i] << ")";
        }
        else
            ++n;
    }
    if (n > 1)
        cout << "^" << n;
    cout << endl;
}

GI GIproduct(const vector<GI>& v)
// iteratively compute a product of Gaussian integers
{
    GI res = GI(1);
    for (GI alpha : v)
        res = res * alpha;
    return res;
}

GI GIpow(const GI& alpha, mpz n)
// calculates alpha^n
{
    if (alpha == GI(0) && n == 0) return GI(0);     // we define 0^0 = 0
    string bin = decToBin(n);                       // get the power n in binary
    GI res = GI(1);
    for (int i = bin.size() - 1; i >= 0; --i) {     // use repeated squarings
        res = res*res;
        if (bin[i] == '1') res = res*alpha;
    }
    return res;
}

GI GImodularExponentiation(GI alpha, mpz n, GI beta)
// calculates alpha^n mod beta, where n is a positive integer
{
    string bin = decToBin(n);               // get the power n in binary
    GI res = GI(1);
    for (int i = bin.size() - 1; i >= 0; --i) {
        res = res*res % beta;
        if (bin[i] == '1') res = res*alpha % beta;
    }
    return res;
}

GI GIprimary(GI alpha)
// returns the unique primary associate of an odd Gaussian integer alpha
{
    if (alpha % GI(1,1) == GI(0))
        throw runtime_error("input must be an odd Gaussian integer");
    if (alpha.real() % 2 == 0)
        alpha = GI(0,1)*alpha;  // make the real part odd
    if (alpha.real() % 4 == 3 || alpha.real() % 4 == -1) {
        if (alpha.imag() % 4 == 2 || alpha.imag() % 4 == -2) return alpha;
        else return GI(-1)*alpha;
    }
    if (alpha.real() % 4 == 1 || alpha.real() % 4 == -3) {
        if (alpha.imag() % 4 == 2 || alpha.imag() % 4 == -2) return GI(-1)*alpha;
        else return alpha;
    }
}

// primes and algorithms
// -----------------------------------------------------------------------------------

bool isGIprime(const GI& alpha)
// check if alpha is a Gaussian prime using the characterization
// note that the function uses a probabilistic primality test, so it may in rare cases
// return false even if alpha is a Gaussian prime
{
    if (is_prop_prime(alpha.norm()) == 2) return true;
    if (alpha.real() != 0 && alpha.imag() == 0)
        if (abs(alpha.real()) % 4 == 3 && is_prop_prime(alpha.real()) == 2) return true;
    if (alpha.real() == 0 && alpha.imag() != 0)
        if (abs(alpha.imag()) % 4 == 3 && is_prop_prime(alpha.imag()) == 2) return true;
    return false;
}

GI GIgcd(const GI& alpha, const GI& beta)
// Euclidean algorithm for Gaussian integers
{
    if (alpha % beta == GI(0))
        return beta;
    return GIgcd(beta, alpha % beta);
}

GI GIgcd(const vector<GI>& GIs)
// the gcd of a vector of Gaussian integers
{
    if (GIs.size() == 0)
        return GI(0);
    if (GIs.size() == 1)
        return GIs[0];

    GI gcd = GIgcd(GIs[0], GIs[1]);
    for (int i = 2; i < GIs.size(); ++i) {
        gcd = GIgcd(gcd, GIs[i]);
    }
    return gcd;
}

bool GIpairwisecoprime(const vector<GI>& GIs)
// returns true, if the GIs in the vector are pairwise coprime
{
    for (int i = 0; i < GIs.size(); ++i) {
        for (int j = i + 1; j < GIs.size(); ++j) {
            if (GIgcd(GIs[i], GIs[j]).norm() != 1) return false;
        }
    }
    return true;
}

vector<GI> GIexgcd(const GI& alpha, const GI& beta)
// Extended Euclidean algorithm for Gaussian integers
// Returns {gcd, s, t} where s and t are coefficients satisfying alpha*s + beta*t = gcd
{
    if (beta == GI(0)) {
        vector<GI> res = {alpha, GI(1), GI(0)};
        return res;
    }
    else {
        vector<GI> temp = GIexgcd(beta, alpha % beta);
        vector<GI> res = {temp[0], temp[2], temp[1] - quo(alpha, beta)*temp[2]};
        return res;
    }
}

GI GIlcm(const GI& alpha, const GI& beta)
// calculate the lcm of two Gaussian integers using the Euclidean algorithm
{
    return (alpha * beta) / GIgcd(alpha, beta);
}

GI GIlcm(const vector<GI>& GIs)
// calculate the lcm of a vector of Gaussian integers
{
    if (GIs.size() == 0)
        return GI(0);
    if (GIs.size() == 1)
        return GIs[0];

    GI lcm = GIlcm(GIs[0], GIs[1]);
    for (int i = 2; i < GIs.size(); ++i) {
        lcm = GIlcm(lcm, GIs[i]);
    }
    return lcm;
}

GI GImodlinearsolve(const GI& alpha, const GI& beta, const GI& gamma)
// returns a solution rho (if one exists) to the equation alpha*rho = beta
// (mod gamma)
{
    vector<GI> s = GIexgcd(alpha, gamma);
    if (beta % s[0] == GI(0)) {
        GI rho = (s[1] * beta)/s[0];
        rho = rho % gamma;
        return rho;
    }
    cout << "No solution to equation: (" << alpha << ")rho" << " = " << beta << " (mod " << gamma << ")" << endl;
    return GI(0);
}

GI GIChRem(const vector<GI> alpha, const vector<GI> moduli)
// The Chinese Remainder Theorem for Z[i], finds rho such that
// alpha[0] = rho mod moduli[0], ..., alpha[n] = rho mod moduli[n] where
// n = alpha.size() = moduli.size()
{
    if (alpha.size() != moduli.size()) {
        cout << "Error, input lengths are not identical" << endl;
        return GI(0);
    }
    if (alpha.size() == 0)
        return GI(0);

    GI rho = alpha[0];      // the solution rho
    GI gamma = moduli[0];   // the product of the moduli
    if (GIpairwisecoprime(moduli)) {
        for (int i = 1; i < alpha.size(); ++i) {
            rho = rho + gamma * GImodlinearsolve(gamma, alpha[i] - rho, moduli[i]);
            gamma = gamma * moduli[i];
        }
        return rho;
    }
    cout << "Error, moduli must be pairwise coprime" << endl;
    return GI(0);
}

vector<GI> GIprimefactor(GI alpha)
// factors alpha as a product of irreducibles and a single unit
{
    vector<mpz> p = primefactor(alpha.norm());
    vector<GI> factors;
    for (int i = 0; i < p.size(); ++i) {
        if (p[i] == 2) {                // if p[i] = 2, 1+i is a factor
            factors.push_back(GI(1,1));
            alpha = alpha / GI(1,1);
        }
        else if (p[i] % 4 == 3) {       // if p[i] % 4 == 3, remove two factors of p[i]
            factors.push_back(p[i]);
            alpha = alpha / p[i];
            ++i;
        }
        else if (p[i] % 4 == 1) {
            mpz n = 2;
            mpz k;
            while(modularExponentiation(n, (p[i]-1)/2, p[i]) == 1) ++n;
            k = modularExponentiation(n, (p[i]-1)/4, p[i]);     // find k, such that k^2 = -1 mod p[i] (and reduce mod p[i])
            GI factor = GIgcd(p[i], {k, 1});
            if (alpha % factor != GI(0)) factor = conj(factor); // if the factor does not divide, the conjugate is the desired factor
            factors.push_back(factor);
            alpha = alpha / factor;
        }
    }
    factors.push_back(alpha);       // when the previous loop is done, alpha is a unit
    return factors;
}

GI QuarticRes(GI alpha, GI pi)
// uses the definition to compute the residue symbol (alpha, pi)_4 using modular exponentiation
// pi is assumed to be prime
{
    if (pi % GI(1,1) == GI(0))
        throw runtime_error("second parameter must be an odd Gaussian prime");
    return GImodularExponentiation(alpha, (pi.norm() - 1)/4, pi);
}

GI QuarticResSymbol(GI alpha, GI beta)
// calculates the quartic residue symbol (alpha, beta)_4 where beta is an odd Gaussian integer
{
    if (beta % GI(1, 1) == GI(0))
        throw runtime_error("second parameter must be an odd Gaussian integer");
    GI res = GI(1);

    while(true) {
        beta = GIprimary(beta);
        alpha %= beta;

        if (alpha == GI(0)) {
            if (beta.norm() != 1) return GI(0);     // alpha and beta have a common factor
            else return res;
        }
        while (alpha % GI(1,1) == GI(0)) {          // remove factors of 1 + i and apply
            alpha /= GI(1,1);                       // the supplementary law for 1 + i
            res *= GIpow(GI(0,1), ((beta.real() - beta.imag() - beta.imag()*beta.imag() - 1)/4) % 4 + 4);
        }
        GI u = alpha / GIprimary(alpha);              // the inverse of the unit for the primary associate
        if (u == GI(-1)) res *= GIpow(GI(0,1), (1 - beta.real()) % 4 + 4);
        if (u == GI(0,1)) res *= GIpow(GI(0,1), ((1 - beta.real())/2) % 4 + 4);
        if (u == GI(0,-1)) res *= GIpow(GI(0,1), (3*(1 - beta.real())/2) % 4 + 4);
        alpha = GIprimary(alpha);
        if ((alpha.real() % 4 == 3 || alpha.real() % 4 == -1) && (beta.real() % 4 == 3 || beta.real() % 4 == -1)) {    // quartic reciprocity
            res *= GI(-1);
        }
        swap(alpha, beta);
    }
}
