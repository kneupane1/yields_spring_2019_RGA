
#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
class Reaction {
 protected:
  std::shared_ptr<Branches12> _data;

  // double _beam_energy = 10.6;
  double _beam_energy = 6.394;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _other;
  std::unique_ptr<TLorentzVector> _neutron;
  std::unique_ptr<TLorentzVector> _boosted_gamma;
  std::unique_ptr<TLorentzVector> _boosted_prot;
  std::unique_ptr<TLorentzVector> _boosted_pip;
  std::unique_ptr<TLorentzVector> _boosted_pim;

  std::unique_ptr<TLorentzVector> _rotated_prot;
  std::unique_ptr<TLorentzVector> _rotated_pip;
  std::unique_ptr<TLorentzVector> _rotated_pim;
  std::unique_ptr<TLorentzVector> _rotated_pim_measured;

  TVector3 _prot_Vect3;
  TVector3 _pip_Vect3;
  TVector3 _pim_Vect3;

  // std::unique_ptr<TLorentzVector> _missingPim;
  std::unique_ptr<TLorentzVector> _boosted_pim_measured;

  bool _is_boosted = false;

  bool _mc = false;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  bool _is_FD = false;
  bool _is_CD = false;
  bool _is_lower_band = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  short _sector = -1;

  float _MM_mProt = NAN;
  float _MM2_mPim = NAN;
  float _MM2_exclusive = NAN;
  float _MM_exclusive = NAN;
  float _excl_Energy = NAN;
  float _MM2_mPip = NAN;
  float _MM2_mProt = NAN;

  float _W = NAN;
  float _Q2 = NAN;
  float _elec_mom = NAN;
  float _elec_E = NAN;
  float _theta_e = NAN;

  float _W_after = NAN;

  float _P_elec = NAN;

  float _theta_star = NAN;
  float _phi_star = NAN;

  float _rec_pim_mom = NAN;
  float _rec_pim_theta = NAN;
  float _rec_pim_phi = NAN;

  float _elec_phi = NAN;
  float _x_mu_phi = NAN;
  float _diff_elec_x_mu_theta = NAN;
  float _diff_elec_x_mu_phi = NAN;

  float _beam_phi = NAN;
  float _diff_beam_x_mu_theta = NAN;
  float _diff_beam_x_mu_phi = NAN;

  float _prot_status = NAN;
  float _pip_status = NAN;
  float _pim_status = NAN;

  int _sectorPim = -1;
  int _sectorPip = -1;
  int _sectorProt = -1;

  float _inv_Ppip = NAN;
  float _inv_Ppim = NAN;
  float _inv_pip_pim = NAN;

  void SetElec();

  /// finished momentum corrections earlier

  double _elec_mom_corrected = NAN;

  double _cx = NAN;
  double _cy = NAN;
  double _cz = NAN;

  float _chi2pid_Ele = NAN;
  float _chi2pid_Prot = NAN;
  float _chi2pid_Pip = NAN;
  float _chi2pid_Pim = NAN;

  //
  static const int CD_SEC = 3;
  static const int FD_SEC = 6;
  double _pip_mom = NAN;

 public:
  Reaction() {};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();
  inline float weight() {
    // return _data->mc_weight();
    return 1.0;
  }
  // Check lists when you swich from mc to exp or vice-versa
  // 1. inline weight function above
  // 2. gamma, _w, _q2 and dpp function in electron four vector set up at reaction.cpp because of momentum corrections
  // for elec included only for exp data
  // 3. turn on the SetMomCorrElec() function on clas12_yields.hpp
  // 4. clas12_yields: auto data = std::make_shared<Branches12>(_chain, true);  turn off true for data
  // 5. from if (data->mc_npart() < 1) to all particle set up im mc events.
  // 6. all mc bank related (generated) output parameters will not work in exp data

  // momentum correction
  void SetMomCorrElec();
  double dpp(float px, float py, float pz, int sec_mom_corr, int ivec);
  double dppC(float px, float py, float pz, int sec, int ivec);
  double Corr_elec_mom();

  inline bool mc() { return _mc; }
  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  float rec_pim_px();
  float rec_pim_py();
  float rec_pim_pz();
  float rec_pim_E();
  float rec_pim_P();
  float rec_pim_mm2();

  float Diff_elec_x_mu_theta();
  float Diff_elec_x_mu_phi();

  float Diff_beam_x_mu_theta();
  float Diff_beam_x_mu_phi();

  float pim_px();
  float pim_py();
  float pim_pz();
  float pim_E();
  float pim_P();

  float pip_px();
  float pip_py();
  float pip_pz();
  float pip_E();

  float prot_px();
  float prot_py();
  float prot_pz();
  float prot_E();

  float elec_px();
  float elec_py();
  float elec_pz();
  inline float elec_mom() { return _elec_mom; }
  inline float elec_En() { return _elec_E; }
  inline float Theta_Elec() { return _theta_e; }

  float beam_px();
  float beam_py();
  float beam_pz();
  float beam_E();

  float target_px();
  float target_py();
  float target_pz();
  float target_E();

  float pim_momentum_corrected();
  float pim_theta_corrected();
  float pim_Phi_corrected();

  float pip_momentum_corrected();
  float pip_theta_corrected();
  float pip_Phi_corrected();

  float prot_momentum_corrected();
  float prot_theta_corrected();
  float prot_Phi_corrected();

  // missingPim
  float pim_momentum();
  float pim_theta_lab();
  float pim_Phi_lab();
  float pim_momentum_measured();
  float pim_theta_lab_measured();
  float pim_Phi_lab_measured();

  float pim_theta_cm();
  float pim_Phi_cm();
  float pim_momentum_cm();
  float pim_theta_cm_measured();
  float pim_Phi_cm_measured();
  float pim_momentum_cm_measured();

  // missingPip
  float pip_momentum();
  float pip_theta_lab();
  float pip_Phi_lab();
  float pip_momentum_measured();
  float pip_theta_lab_measured();
  float pip_Phi_lab_measured();

  // missingProt
  float prot_momentum();
  float prot_theta_lab();
  float prot_Phi_lab();
  float prot_momentum_measured();
  float prot_theta_lab_measured();
  float prot_Phi_lab_measured();

  void boost();

  inline float Theta_star() { return _theta_star; }
  inline float Phi_star() { return _phi_star; }

  void CalcMissMass();
  float MM2_mPim();
  float MM2_exclusive();
  float Energy_excl();
  float MM2_mPip();
  float MM2_mProt();
  float MM_mPim();
  float MM_mPip();
  float MM_mProt();
  float MM_exclusive();

  float inv_Ppip();
  float inv_Ppim();
  float inv_Pippim();

  void invMassPpip();
  void invMassPpim();
  void invMasspippim();

  virtual std::string CsvHeader();
  virtual std::string ReacToCsv();

  inline float chi2pid_Elec() { return _chi2pid_Ele; }
  inline float chi2pid_Prot() { return _chi2pid_Prot; }
  inline float chi2pid_Pip() { return _chi2pid_Pip; }
  inline float chi2pid_Pim() { return _chi2pid_Pim; }

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  inline float W_after() { return _W_after; }

  float_t scalar_triple_product();

  inline short sec() { return _data->dc_sec(0); }
  inline short pimSec() { return _sectorPim; }
  inline short pipSec() { return _sectorPip; }
  inline short protSec() { return _sectorProt; }

  inline int det() { return abs(_data->status(0) / 1000); }

  inline bool Inclusive() {
    bool _channelIncl = true;
    _channelIncl &= (_hasE);
    return _channelIncl;
  }

  inline bool TwoPion_missingPim() {
    bool _channelTwoPi = true;
    _channelTwoPi &= ((_numProt == 1 && _numPip == 1) && (_hasE && _hasP && _hasPip));
    return _channelTwoPi;
  }

  inline bool TwoPion_exclusive() {
    bool _channelTwoPi_excl = true;

    _channelTwoPi_excl &= ((_numProt == 1 && _numPip == 1 && _numPim == 1) &&
                           (_hasE && _hasP && _hasPip && _hasPim /*&& !_hasNeutron && !_hasOther*/));
    return _channelTwoPi_excl;
  }
  inline bool TwoPion_missingPip() {
    bool _channelTwoPi_mpip = true;

    _channelTwoPi_mpip &=
        ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && _hasPim /*&&!_hasPip && !_hasNeutron && !_hasOther*/));
    return _channelTwoPi_mpip;
  }
  inline bool TwoPion_missingProt() {
    bool _channelTwoPi_mprot = true;
    _channelTwoPi_mprot &=
        ((_numPip == 1 && _numPim == 1) && (_hasE && _hasPip && _hasPim /*&&!_hasP  && !_hasOther*/));
    return _channelTwoPi_mprot;
  }

  const TLorentzVector &e_mu() { return *_beam; }
  const TLorentzVector &e_mu_prime() { return *_elec; }
  const TLorentzVector &gamma() { return *_gamma; }
};

class MCReaction : public Reaction {
 private:
  float _weight_mc = NAN;
  float _W_mc = NAN;
  float _Q2_mc = NAN;

  std::unique_ptr<TLorentzVector> _elec_mc;
  std::unique_ptr<TLorentzVector> _gamma_mc;
  std::unique_ptr<TLorentzVector> _prot_mc;
  std::unique_ptr<TLorentzVector> _pip_mc;
  std::unique_ptr<TLorentzVector> _pim_mc;
  std::unique_ptr<TLorentzVector> _other_mc;

  float _MM2_exclusive_mc = NAN;
  float _excl_Energy_mc = NAN;

  float _rec_x_mu_mom_mc = NAN;
  float _rec_x_mu_theta_mc = NAN;
  float _x_mu_phi_mc = NAN;

  float _beam_phi_mc = NAN;
  float _elec_phi_mc = NAN;
  float _diff_elec_x_mu_theta_mc = NAN;
  float _diff_elec_x_mu_phi_mc = NAN;

  float _diff_beam_x_mu_theta_mc = NAN;
  float _diff_beam_x_mu_phi_mc = NAN;

  float _elec_mom_mc = NAN;
  float _elec_E_mc = NAN;
  float _theta_e_mc = NAN;

 public:
  void SetMCProton(int i);
  void SetMCPip(int i);
  void SetMCPim(int i);
  void SetMCOther(int i);

  MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  void SetMCElec();
  inline float weight() { return _data->mc_weight(); }

  inline float elec_mom_mc() { return _elec_mom_mc; }
  inline float elec_En_mc() { return _elec_E_mc; }
  inline float Theta_Elec_mc() { return _theta_e_mc; }
  inline float W_mc() { return _W_mc; }
  inline float Q2_mc() { return _Q2_mc; }

  float pim_mom_mc_gen();
  float pip_mom_mc_gen();
  float prot_mom_mc_gen();

  float pim_theta_mc_gen();
  float pip_theta_mc_gen();
  float prot_theta_mc_gen();

  float pim_phi_mc_gen();
  float pip_phi_mc_gen();
  float prot_phi_mc_gen();

  void CalcMissMass_mc();

  std::string CsvHeader();
  std::string ReacToCsv();
};

#endif
