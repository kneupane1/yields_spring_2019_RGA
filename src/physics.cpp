
#include "physics.hpp"

namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  p_mu.SetVectM(p_mu_3, MASS_P);
  return (p_mu + q_mu).Mag();
}

// overload with 4 vectors
double xb_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  double Q2 = Q2_calc(e_mu, e_mu_prime);
  TLorentzVector q = e_mu - e_mu_prime;
  TLorentzVector target(0, 0, 0, MASS_P);
  return (Q2 / (2 * (q.Dot(target))));
}
double vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

double deltat(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r) {
  double relatavistic_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return electron_vertex_time - vertex_time(sc_t, sc_r, relatavistic_beta);
}

std::shared_ptr<TLorentzVector> fourVec(double px, double py, double pz, double mass) {
  auto x = std::make_shared<TLorentzVector>();
  x->SetXYZM(px, py, pz, mass);
  return x;
}
float invTan(const float &y, const float &x) {
  if (x > 0 && y > 0)
    return atan(y / x);  // 1st Quad.
  else if (x < 0 && y > 0)
    return atan(y / x) + PI;  // 2nd Quad
  else if (x < 0 && y < 0)
    return atan(y / x) + PI;  // 3rd Quad
  else if (x > 0 && y < 0)
    return atan(y / x) + 2 * PI;  // 4th Quad
  else if (x == 0 && y > 0)
    return PI / 2;
  else if (x == 0 && y < 0)
    return 3 * PI / 2;
  return NAN;
}

float phi_boosted(const std::shared_ptr<TLorentzVector> &vec) { return invTan(vec->Py(), vec->Px()); }

}  // namespace physics
