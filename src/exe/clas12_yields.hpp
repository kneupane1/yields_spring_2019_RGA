
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  float beam_energy = 10.2;
  if (std::is_same<CutType, rga_Cuts>::value) {
    beam_energy = 10.2;
  } else if (std::is_same<CutType, uconn_Cuts>::value) {
    beam_energy = 10.2;
    // } else if (std::is_same<CutType, rgf_Cuts>::value) {
    //         beam_energy = rgf_E0;
    // }
    // else if (std::is_same<CutType, rgk_Cuts>::value) {
    //         beam_energy = rgk_E0;
  }

  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  // for sim data use it
  // auto data = std::make_shared<Branches12>(_chain, true);
  // for exp data use it
  auto data = std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // for (size_t current_event = 0; current_event < 350; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);

    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    int statusPim = -9999;
    int statusPip = -9999;
    int statusProt = -9999;

    // if (data->mc_npart() < 1) continue;

    // // If we pass electron cuts the event is processed
    // total++;

    // // Make a reaction class from the data given
    // auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

    // for (int part = 1; part < data->mc_npart(); part++) {
    //   // Check particle ID's and fill the reaction class

    //   if (data->mc_pid(part) == PIP) {
    //     mc_event->SetMCPip(part);
    //   } else if (data->mc_pid(part) == PROTON) {
    //     mc_event->SetMCProton(part);
    //   } else if (data->mc_pid(part) == PIM) {
    //     mc_event->SetMCPim(part);
    //     // } else {
    //     //   mc_event->SetMCOther(part);
    //   }
    // }

    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<uconn_Cuts>(data);
    // auto cuts = std::make_shared<rga_Cuts>(data);
    if (!cuts->ElectronCuts()) continue;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy);
    // event->SetMomCorrElec();

    // // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);

      // Check particle ID's and fill the reaction class
      if (cuts->IsProton(part)) {
        event->SetProton(part);
        statusProt = abs(data->status(part));
        // if (statusProt < 4000 && statusProt > 2000) sectorProt = data->dc_sec(part);
        // std::cout << "_prot px : " << data->px(part) << "_prot py : " << data->py(part) << "_prot pz : " <<
        // data->pz(part)
        //           << "_prot E : " << MASS_P << std::endl;

      } else if (cuts->IsPip(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetPip(part);
          statusPip = abs(data->status(part));
          // if (statusPip<4000 && statusPip> 2000) sectorPip = data->dc_sec(part);
        }
      } else if (cuts->IsPim(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetPim(part);
          statusPim = abs(data->status(part));
          // if (statusPim < 4000 && statusPim > 2000) sectorPim = data->dc_sec(part);
        }
      } else {
        event->SetOther(part);
      }
    }

    // if (event->TwoPion_missingPim() || event->TwoPion_missingPip() || event->TwoPion_missingProt() ||
    //     event->TwoPion_exclusive()) {
    // if (event->TwoPion_missingPim()) {
      // if (event->TwoPion_missingPip()) {
      // if (event->TwoPion_missingProt()) {
      if (event->TwoPion_exclusive()) {
      // if (event->Inclusive()) {
      // if (event->W() > 1.3 && event->W() < 2.5 && event->Q2() > 1.5 && event->Q2() < 10.5){
      // &&
      // abs(event->MM2_exclusive()) < 0.03 && abs(event->Energy_excl()) < 0.3){
      // &&(event->pim_Phi_lab() > 330 || event->pim_Phi_lab() < 30)) {
      //   //&&
      //   // abs(event->MM2_exclusive()) < 0.03) {
      //   // total++;
      csv_data output;

      // // //// using exclusive topology ...................................

      // output.electron_sector = event->sec();
      output.pim_sec = event->pimSec();
      output.pip_sec = event->pipSec();
      output.prot_sec = event->protSec();
      // output.w = event->W();
      // output.q2 = event->Q2();
      // output.w_had = event->w_hadron();
      // // output.w_diff = event->w_difference();
      // output.w_had_corr = event->w_hadron_corr();
      // // output.w_diff_corr = event->w_difference_corr();

      // output.elec_mom = event->elec_mom();
      // output.elec_energy = event->elec_En();
      // output.elec_theta = event->Theta_Elec();
      // output.corr_elec_mom = event->Corr_elec_mom();
      // output.scalar_product = event->scalar_triple_product();

      // // //   // // for generated case
      //   output.w_mc = mc_event->W_mc();
      //   output.q2_mc = mc_event->Q2_mc();

      //   output.elec_mom_mc = mc_event->elec_mom_mc();
      //   output.elec_energy_mc = mc_event->elec_En_mc();
      //   output.elec_theta_mc = mc_event->Theta_Elec_mc();
      //   output.weight_exclusive = mc_event->weight();

      //   // // for energy loss corrections : gen
        // output.gen_prot_mom = (mc_event->prot_mom_mc_gen());
        // output.gen_prot_theta = (mc_event->prot_theta_mc_gen());
        // output.gen_prot_phi = (mc_event->prot_phi_mc_gen());

        // output.gen_pip_mom = (mc_event->pip_mom_mc_gen());
        // output.gen_pip_theta = (mc_event->pip_theta_mc_gen());
        // output.gen_pip_phi = (mc_event->pip_phi_mc_gen());

        // output.gen_pim_mom = (mc_event->pim_mom_mc_gen());
        // output.gen_pim_theta = (mc_event->pim_theta_mc_gen());
        // output.gen_pim_phi = (mc_event->pim_phi_mc_gen());

        // // // // missing
        output.prot_mom_mProt = event->prot_momentum();
        output.prot_theta_mProt = event->prot_theta_lab();
        output.prot_phi_mProt = event->prot_Phi_lab();

        output.pip_mom_mPip = event->pip_momentum();
        output.pip_theta_mPip = event->pip_theta_lab();
        output.pip_phi_mPip = event->pip_Phi_lab();

        output.pim_mom_mPim = event->pim_momentum();
        output.pim_theta_mPim = event->pim_theta_lab();
        output.pim_phi_mPim = event->pim_Phi_lab();

        // // // recon mes
        // output.prot_mom_exclusive = event->prot_momentum_corrected();
        // output.prot_theta_exclusive = event->prot_theta_corrected();
        // output.prot_phi_exclusive = event->prot_Phi_corrected();
        output.prot_mom_exclusive = event->prot_momentum_measured();
        output.prot_theta_exclusive = event->prot_theta_lab_measured();
        output.prot_phi_exclusive = event->prot_Phi_lab_measured();
        output.prot_dcr1theta_exclusive = event->thetaDCr1Prot();

        output.prot_mom_corr = event->prot_momentum_corrected();
        // // output.prot_theta_corr = event->prot_theta_corrected();
        // // output.prot_phi_corr = event->prot_Phi_corrected();

        output.pip_mom_exclusive = event->pip_momentum_measured();
        output.pip_theta_exclusive = event->pip_theta_lab_measured();
        output.pip_phi_exclusive = event->pip_Phi_lab_measured();
        output.pip_dcr1theta_exclusive = event->thetaDCr1Pip();

        output.pip_mom_corr = event->pip_momentum_corrected();
        // // output.pip_theta_corr = event->pip_theta_corrected();
        // // output.pip_phi_corr = event->pip_Phi_corrected();

        output.pim_mom_exclusive = event->pim_momentum_measured();
        output.pim_theta_exclusive = event->pim_theta_lab_measured();
        output.pim_phi_exclusive = event->pim_Phi_lab_measured();
        output.pim_dcr1theta_exclusive = event->thetaDCr1Pim();

        output.pim_mom_corr = event->pim_momentum_corrected();
        // // output.pim_theta_corr = event->pim_theta_corrected();
        // // output.pim_phi_corr = event->pim_Phi_corrected();

        output.mm2_mProt = event->MM2_mProt();
        output.mm2_mProt_corr = event->MM2_mProt_corr();
        output.mm2_mPip = event->MM2_mPip();
        output.mm2_mPip_corr = event->MM2_mPip_corr();
        output.mm2_mPim = event->MM2();
        output.mm2_mPim_corr = event->MM2_mPim_corr();

        output.mm2_exclusive_at_zero = event->MM2_exclusive();
        output.energy_x_mu = event->Energy_excl();

        output.status_Pim = statusPim;
        output.status_Pip = statusPip;
        output.status_Prot = statusProt;
        // output.inv_ppip = event->inv_Ppip();
        // output.inv_ppim = event->inv_Ppim();
        // output.inv_pip_pim = event->inv_Pippim();

        output.weight_exclusive = event->weight();

        // /// ..........................................

        // //   // // for mom correction pim
        // //   // output.pim_mom_mPim = event->pim_momentum();
        // //   // output.pim_mom_exclusive = event->pim_momentum_measured();
        // //   // output.pim_mom_corr = event->pim_momentum_corrected();

        // //   // output.pim_theta_mPim = event->pim_theta_lab();
        // //   // output.pim_theta_exclusive = event->pim_theta_lab_measured();
        // //   // output.pim_theta_corr = event->pim_theta_corrected();

        // //   // output.pim_phi_mPim = event->pim_Phi_lab();
        // //   // output.pim_phi_exclusive = event->pim_Phi_lab_measured();
        // //   // output.pim_phi_corr = event->pim_Phi_corrected();

        // //   // // for mom correction pip
        // //   // output.pip_mom_mPip = event->pip_momentum();
        // //   // output.pip_mom_exclusive = event->pip_momentum_measured();
        // //   // output.pip_mom_corr = event->pip_momentum_corrected();

        // //   // output.pip_theta_mPip = event->pip_theta_lab();
        // //   // output.pip_theta_exclusive = event->pip_theta_lab_measured();
        // //   // output.pip_theta_corr = event->pip_theta_corrected();

        // //   // output.pip_phi_mPip = event->pip_Phi_lab();
        // //   // output.pip_phi_exclusive = event->pip_Phi_lab_measured();
        // //   // output.pip_phi_corr = event->pip_Phi_corrected();

        // //   // // for mom correction prot
        // //   // output.prot_mom_mProt = event->prot_momentum();
        // //   // output.prot_mom_exclusive = event->prot_momentum_measured();
        // //   // output.prot_mom_corr = event->prot_momentum_corrected();

        // //   // output.prot_theta_mProt = event->prot_theta_lab();
        // //   // output.prot_theta_exclusive = event->prot_theta_lab_measured();
        // //   // output.prot_theta_corr = event->prot_theta_corrected();

        // //   // output.prot_phi_mProt = event->prot_Phi_lab();
        // //   // output.prot_phi_exclusive = event->prot_Phi_lab_measured();
        // //   // output.prot_phi_corr = event->prot_Phi_corrected();

        // //   // output.mm2_exclusive_at_zero = event->MM2_exclusive();
        // //   // output.energy_x_mu = event->Energy_excl();
        // //   // output.weight_exclusive = event->weight();

        // // // // // // // // // mPim .......................................

        // output.pim_mom_mPim = event->pim_momentum();
        // output.pim_theta_mPim = event->pim_theta_lab();
        // output.pim_phi_mPim = event->pim_Phi_lab();
        // output.mm2_mPim = event->MM2();
        // // output.mm2_mPim_corr = event->MM2_mPim_corr();
        // output.weight_mPim = event->weight();

        // // // // // for rec pim
        // // // // // output.elec_mom = event->elec_mom();
        // // // // // output.corr_elec_mom = event->Corr_elec_mom();

        // // output.prot_mom_mProt = event->prot_momentum();
        // // output.prot_theta_mProt = event->prot_theta_lab();
        // // output.prot_phi_mProt = event->prot_Phi_lab();
        // // output.mm2_mProt = event->MM2();

        // // // // for mes prot
        // output.scalar_product = event->scalar_triple_product();
        // output.prot_mom_exclusive = event->prot_momentum_measured();
        // output.prot_theta_exclusive = event->prot_theta_lab_measured();
        // output.prot_phi_exclusive = event->prot_Phi_lab_measured();
        // output.mm2_exclusive = event->MM2();
        // output.mm2_exclusive_at_zero = event->MM2_exclusive();
        // output.energy_x_mu = event->Energy_excl();
        // // output.mm2_mPip = event->MM2_mPip();
        // // output.mm2_mProt = event->MM2_mProt();

        // // // output.diff_ex_theta = event->Diff_elec_x_mu_theta();
        // // // output.diff_ex_phi = event->Diff_elec_x_mu_phi();
        // // // output.diff_bx_theta = event->Diff_beam_x_mu_theta();
        // // // output.diff_bx_phi = event->Diff_beam_x_mu_phi();

        // // output.status_Pim = statusPim;
        // // output.status_Pip = statusPip;
        // // output.status_Prot = statusProt;

        // output.weight_exclusive = event->weight();

        // // for generated case
        // output.w_mc = mc_event->W_mc();
        // output.q2_mc = mc_event->Q2_mc();

        // output.x_mu_mom_exclusive = mc_event->x_mu_momentum_mc();
        // output.x_mu_theta_exclusive = mc_event->x_mu_theta_lab_mc();
        // output.x_mu_phi_exclusive = mc_event->x_mu_Phi_lab_mc();
        // output.mm2_exclusive_at_zero = mc_event->MM2_exclusive_mc();
        // output.energy_x_mu = mc_event->Energy_excl_mc();

        // output.diff_ex_theta = mc_event->Diff_elec_x_mu_theta_mc();
        // output.diff_ex_phi = mc_event->Diff_elec_x_mu_phi_mc();
        // output.diff_bx_theta = mc_event->Diff_beam_x_mu_theta_mc();
        // output.diff_bx_phi = mc_event->Diff_beam_x_mu_phi_mc();

        // output.weight_exclusive = mc_event->weight();

        // mPip  .......................................
        /*
                     output.pip_mom_mPip = event->pip_momentum();
                     output.pip_theta_mPip = event->pip_theta_lab();
                     output.pip_phi_mPip = event->pip_Phi_lab();
                     output.mm2_mPip = event->MM2_mPip();
                     output.mm2_mPip_corr = event->MM2_mPip_corr();
                     output.weight_mPip = event->weight();
        */
        /*
              output.scalar_product = event->scalar_triple_product();
              output.pip_mom_exclusive = event->pip_momentum_measured();
              output.pip_theta_exclusive = event->pip_theta_lab_measured();
              output.pip_phi_exclusive = event->pip_Phi_lab_measured();
              output.mm2_exclusive = event->MM2_mPip();
              output.weight_exclusive = event->weight();
        */
        // mProt .......................................
        /*
                output.prot_mom_mProt = event->prot_momentum();
                output.prot_theta_mProt = event->prot_theta_lab();
                output.prot_phi_mProt = event->prot_Phi_lab();
                output.mm2_mProt = event->MM2_mProt();
                output.mm2_mProt_corr = event->MM2_mProt_corr();
                output.weight_mProt = event->weight();
        */
        /*
              output.scalar_product = event->scalar_triple_product();
              output.prot_mom_exclusive = event->prot_momentum_measured();
              output.prot_theta_exclusive = event->prot_theta_lab_measured();
              output.prot_phi_exclusive = event->prot_Phi_lab_measured();
              output.mm2_exclusive = event->MM2_mProt();
              output.weight_exclusive = event->weight();

        */
        // std::cout << "beam px  " << event->beam_px() << std::endl;
        // std::cout << "beam py  " << event->beam_py() << std::endl;
        // std::cout << "beam pz  " << event->beam_pz() << std::endl;
        // std::cout << "beam E  " << event->beam_E() << std::endl;

        // std::cout << "elec px  " << event->elec_px() << std::endl;
        // std::cout << "elec py  " << event->elec_py() << std::endl;
        // std::cout << "elec pz  " << event->elec_pz() << std::endl;
        // std::cout << "elec E  " << event->elec_E() << std::endl;

        // std::cout << "target px  " << event->target_px() << std::endl;
        // std::cout << "target py  " << event->target_py() << std::endl;
        // std::cout << "target pz  " << event->target_pz() << std::endl;
        // std::cout << "target E  " << event->target_E() << std::endl;

        // std::cout << "pip px  " << event->pip_px() << std::endl;
        // std::cout << "pip py  " << event->pip_py() << std::endl;
        // std::cout << "pip pz  " << event->pip_pz() << std::endl;
        // std::cout << "pip E  " << event->pip_E() << std::endl;

        // std::cout << "prot px  " << event->prot_px() << std::endl;
        // std::cout << "prot py  " << event->prot_py() << std::endl;
        // std::cout << "prot pz  " << event->prot_pz() << std::endl;
        // std::cout << "prot E  " << event->prot_E() << std::endl;

        // std::cout << "W  " << event->W() << std::endl;
        // std::cout << "MM2_root_calc " << event->MM2() << std::endl;
        // std::cout << "MM2_our_hand_calc  " << event->rec_pim_mm2() << std::endl;
        // std::cout << "root_rec_pim P " << event->pim_momentum() << std::endl;
        // std::cout << "root_rec_pim E " << event->MM() << std::endl;
        // std::cout << "hand_rec_pim P " << event->rec_pim_P() << std::endl;
        // std::cout << "mes_pim P " << event->pim_P() << std::endl;

        // std::cout << "hand_rec_pim px " << event->rec_pim_px() << std::endl;
        // std::cout << "hand_rec_pim py " << event->rec_pim_py() << std::endl;
        // std::cout << "hand_rec_pim pz " << event->rec_pim_pz() << std::endl;
        // std::cout << "hand_rec_pim E " << event->rec_pim_E() << std::endl;

        // std::cout << "mes_pim px " << event->pim_px() << std::endl;
        // std::cout << "mes_pim py " << event->pim_py() << std::endl;
        // std::cout << "mes_pim pz " << event->pim_pz() << std::endl;
        // std::cout << "mes_pim E " << event->pim_E() << std::endl;

        _sync->write(output);
    // }
      }
  }
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // Return the total number of events
  return num_of_events;
}
#endif

/*
   ep -> e x+

   W for all events
   W for 2 particles
   W for 2 Part 2nd positive
   hist Phi_e - Phi_pos ~ 90
   W for cut around 90
   pos_mom vs pos_theta
   theta_p_pos_calc_from_electron - theta_pos_measured

   calc theta from magnitude of pos momentum

 */
