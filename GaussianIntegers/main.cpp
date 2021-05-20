#include <iostream>
#include "gaussian_integers.h"

using namespace std;

int main()
{
    GI alpha {1, 2};
    GI beta {-3, 5};
    GI gamma {1284, -418764};
    cout << alpha + beta << endl;
    cout << alpha - beta << endl;
    cout << alpha * beta << endl;
    cout << GIgcd(alpha, beta) << endl;
    cout << gamma << " = "; printGIproduct(GIprimefactor(gamma));
    cout << QuarticResSymbol(gamma, alpha);
    return 0;
}
