#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  short electron_sector;
  short pim_sec;
  short pip_sec;
  short prot_sec;
  float w;
  float q2;
  float w_mc;
  float q2_mc;
  float w_had;
  float w_diff;
  float w_had_corr;
  float w_diff_corr;
  float w_after;

  float elec_mom;
  float elec_energy;
  float elec_theta;
  float elec_mom_mc;
  float elec_energy_mc;
  float elec_theta_mc;

  float corr_elec_mom;
  float mom_part;
  float beta_part;
  float charge_part;
  float pid_part;
  int status_part;

  float scalar_product;
  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;
  float mm2_mPim;
  float mm2_mPim_corr;
  float weight_mPim;
  float pim_mom_exclusive;

  float pim_theta_exclusive;
  float pim_phi_exclusive;
  float mm2_exclusive;
  float mm2_exclusive_at_zero;
  float weight_exclusive;
  float mm_exclusive_at_zero;

  float pip_mom_mPip;
  float pip_theta_mPip;
  float pip_phi_mPip;
  float mm2_mPip;
  float mm2_mPip_corr;
  float weight_mPip;
  float pip_mom_exclusive;
  float pip_theta_exclusive;
  float pip_phi_exclusive;
  float energy_x_mu;

  float prot_mom_mProt;
  float prot_theta_mProt;
  float prot_phi_mProt;
  float mm2_mProt;
  float mm2_mProt_corr;
  float weight_mProt;

  float prot_mom_exclusive;
  float prot_theta_exclusive;
  float prot_phi_exclusive;
  float prot_dcr1theta_exclusive;
  float pip_dcr1theta_exclusive;
  float pim_dcr1theta_exclusive;

  float x_mu_mom_exclusive;
  float x_mu_theta_exclusive;
  float x_mu_phi_exclusive;

  float gen_pim_mom;
  float gen_pim_theta;
  float gen_pim_phi;

  float gen_pip_mom;
  float gen_pip_theta;
  float gen_pip_phi;

  float gen_prot_mom;
  float gen_prot_theta;
  float gen_prot_phi;

  int status_Pim;
  int status_Pip;
  int status_Prot;

  float inv_ppip;
  float inv_ppim;
  float inv_pip_pim;

  float chi2pid_e;
  float chi2pid_p;
  float chi2pid_pip;
  float chi2pid_pim;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    return "sec_pim,sec_pip,sec_prot,w,q2,mm_exclusive_at_zero,mm2_exclusive_at_zero,energy_x_mu,mm_mProt,mm_mPip,mm_"
           "mPim,"
           "inv_pPip,inv_pPim,inv_pipPim,"
           "status_Pim,"
           "status_Pip,status_Prot,weight";

    // mccase
    // return "w,q2,w_had,w_mc,q2_mc,mm2_exclusive_at_zero,energy_x_mu,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    ////.......................................
    os << std::setprecision(1);

    // // os << data.electron_sector << ",";
    os << data.pim_sec << ",";
    os << data.pip_sec << ",";
    os << data.prot_sec << ",";

    os << std::setprecision(7);

    os << data.w << ",";
    os << data.q2 << ",";

    // // // // measured
    // os << std::setprecision(10);

    // os << data.prot_mom_exclusive << ",";
    // os << std::setprecision(7);

    // os << data.prot_theta_exclusive << ",";
    // os << data.prot_phi_exclusive << ",";

    // os << std::setprecision(10);

    // os << data.pip_mom_exclusive << ",";
    // os << std::setprecision(7);

    // os << data.pip_theta_exclusive << ",";
    // os << data.pip_phi_exclusive << ",";

    // os << data.pim_mom_exclusive << ",";
    // // os << std::setprecision(7);

    // os << data.pim_theta_exclusive << ",";
    // os << data.pim_phi_exclusive << ",";
    os << data.mm_exclusive_at_zero << ",";
    os << data.mm2_exclusive_at_zero << ",";
    os << data.energy_x_mu << ",";
    os << data.mm2_mProt << ",";
    // os << data.mm2_mProt_corr << ",";

    os << data.mm2_mPip << ",";
    // os << data.mm2_mPip_corr << ",";
    os << data.mm2_mPim << ",";

    // os << data.mm2_mPim << ",";
    // os << data.mm2_mPim_corr << ",";

    // os << std::setprecision(7);
    os << data.inv_ppip << ",";
    os << data.inv_ppim << ",";
    os << data.inv_pip_pim << ",";

    os << std::setprecision(3);

    os << data.status_Pim << ",";
    os << data.status_Pip << ",";
    os << data.status_Prot << ",";

    // os << data.chi2pid_e<< ",";
    // os << data.chi2pid_p << ",";
    // os << data.chi2pid_pip << ",";
    // os << data.chi2pid_pim << ",";

    os << std::setprecision(1);

    os << data.weight_exclusive << ",";

    return os;
  }
};

#endif
