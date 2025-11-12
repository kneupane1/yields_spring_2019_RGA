
#include "cuts.hpp"
#include <iostream>
#include "TFile.h"
#include "reaction.hpp"

Cuts::Cuts(const std::shared_ptr<Branches12> &data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : _data(data), _dt(dt) {}
Cuts::~Cuts() {}

bool Cuts::IsPip(int i) {
  if (_data->gpart() <= i) return false;
  bool _pip = true;
  // _pip &= (_data->pid(i) == 0);
  // _pip &= (_data->sc_cnd_layer(i) ==3 );
  // // _pip &= ((_data->p(i) < 0.5) || (_data->sc_extras_dedx(i) < (-3.43 * (_data->p(i)) + 5.8)));
  // _pip &= (_data->charge(i) == POSITIVE);

  _pip &= (_data->charge(i) == POSITIVE);
  _pip &= (_data->pid(i) == PIP);
  // _pip &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.4);
  // _pip &= !(abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.2);
  //   _pip &= (_data->p(i) > 0.2);
  _pip &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000);

  // // min/max mom cuts
  // if (abs(_data->status(i)) < 4000) {
  //   _pip &= (_data->p(i) > 0.3);
  //   _pip &= (_data->p(i) < 5.0);
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _pip &= (_data->p(i) > 0.0);
  //   _pip &= (_data->p(i) < 3.0);
  // }
  ////_pip &= (abs(_data->chi2pid(i)) < 0.5);

  //_pip &= DC_fiducial_cut_theta_phi(i);
  //_pip &= Hadron_Delta_vz_cut(i);
  //_pip&= Hadron_Chi2pid_cut(i);*/

  return _pip;
}
bool Cuts::IsProton(int i) {
  if (_data->gpart() <= i) return false;
  bool _proton = true;
  // _proton &= (_data->pid(i) == 0);
  // _proton &= (_data->charge(i) == POSITIVE);
  // _proton &= (_data->sc_cnd_layer(i) == 3);
  // _proton &= (_data->p(i) > 0.5);
  // _proton &= (_data->sc_extras_dedx(i) > (-3.43 * (_data->p(i)) + 5.8));

  _proton &= (_data->charge(i) == POSITIVE);
  _proton &= (_data->pid(i) == PROTON);
  // _proton &= (abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.4);
  // _proton &= !(abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
  _proton &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000);
  // min/max mom cuts
  // if (abs(_data->status(i)) < 4000) {
  //   _proton &= (_data->p(i) > 0.4);
  //   _proton &= (_data->p(i) < 5.0);
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _proton &= (_data->p(i) > 0.0);
  //   _proton &= (_data->p(i) < 3.0);
  // }
  //   _proton &= (_data->p(i) > 0.2);
  // //_proton &= (abs(_data->chi2pid(i)) < 0.5);*/
  return _proton;
}
bool Cuts::IsPim(int i) {
  if (_data->gpart() <= i) return false;
  bool _pim = true;
  //   _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (_data->pid(i) == PIM);
  // _pim &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.5);
  _pim &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000);
  // min/max mom cuts
  // if (abs(_data->status(i)) < 4000){
  //   _pim &= (_data->p(i) > 0.4);
  //   _pim &= (_data->p(i) < 5.0);
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _pim &= (_data->p(i) > 0.0);
  //   _pim &= (_data->p(i) < 3.0);
  // }

  // //_pim &= (abs(_data->chi2pid(i)) < 0.5);

  //_pim &= DC_fiducial_cut_theta_phi(i);
  //_pim &= Hadron_Delta_vz_cut(i);
  //_pim &= Hadron_Chi2pid_cut(i);

  return _pim;
}

// /////////////////////// Pid_Cuts ///////////////////////
bool Pid_Cuts::ElectronCuts() {
  bool cut = true;
  cut &= (_data->gpart() > 0);
  if (!cut) return false;

  cut &= (_data->gpart() < 20);
  //
  cut &= (_data->charge(0) == NEGATIVE);
  cut &= (_data->pid(0) == ELECTRON);
  // cut &= (_data->p(0) > 1.50);
  cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
  // cut &= (abs(_data->chi2pid(0)) < 3);  ////////////// check it.......
  // cut &= CC_nphe_cut();
  // cut &= EC_outer_vs_EC_inner_cut();
  // cut &= EC_sampling_fraction_cut();
  // cut &= EC_hit_position_fiducial_cut_homogeneous();
  // cut &= DC_fiducial_cut_XY();
  // cut &= DC_z_vertex_cut();
  return cut;
}
bool Pid_Cuts::HadronsCuts(int i) {
  bool cut = true;
  // if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000) cut &= DC_fiducial_cut_theta_phi(i);
  // cut &= Hadron_Delta_vz_cut(i);
  // cut &= Hadron_Chi2pid_cut(i);
  return cut;
}

///////////////////// Pid_Cuts ///////////////////////
