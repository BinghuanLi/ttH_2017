#include "../interference/utils_functions.h"
 
double deltaPhi(double phi1, double phi2){
    double result = phi1 - phi2;
    while (result > M_PI) result -= 2*M_PI;
    while (result <= -M_PI) result += 2*M_PI;
    return result;
}


double deltaEta(double eta1, double eta2){
    return (eta1-eta2);
};


double deltaR(double dphi, double deta){
    return sqrt(pow(dphi,2)+pow(deta,2));
};


