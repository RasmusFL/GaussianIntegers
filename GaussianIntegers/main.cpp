#include <iostream>
#include "gaussian_integers.h"

using namespace std;

int main()
{
    GI alpha {1, 2};
    GI beta {-3, 5};
    cout << alpha + beta << endl;
    cout << alpha - beta << endl;
    cout << alpha * beta << endl;
    cout << beta / alpha << endl;
    cout << beta % alpha << endl;
    cout << GIgcd({347, 89}, {117, -547}) << endl;
    vector<GI> v = GIexgcd({347, 89}, {117, -547});
    for (GI a : v)
        cout << a << endl;
    cout << GImodlinearsolve({347, 89}, {117, -547}, {2003}) << endl;
    GImodlinearsolve({347, 89}, {117, -547}, {1, 2});
    printGIproduct(GIprimefactor({1284, -418764}));
    cout << QuarticResSymbol({347, 89}, {5, 6}) << endl;
    cout << isGIprime({5, 6});
    return 0;
}
