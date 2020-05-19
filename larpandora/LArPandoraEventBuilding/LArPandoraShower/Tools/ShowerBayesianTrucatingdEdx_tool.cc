//############################################################################
//### Name:        ShowerBayesianTrucatingdEdx                             ###
//### Author:      Dominic Batker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Recursively adds values from the dEdx vectors and stops ###
//###              when the probability of getting that dEdx value is too  ###
//###              low. This is done for both the a electron prior and     ###
//###              photon prior. The posterior is calculated and the prior ###
//###              with the highest posterior is taken. Currently can      ###
//###              only be used with the sliding calo dEdx                 ###
//############################################################################

#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PFParticle.h"

#include "TFile.h"

#include <vector>
#include <string>

namespace ShowerRecoTools {


  class ShowerBayesianTrucatingdEdx: public IShowerTool {

    public:

      ShowerBayesianTrucatingdEdx(const fhicl::ParameterSet& pset);

      ~ShowerBayesianTrucatingdEdx();

      //Generic Direction Finder
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      double CalculatePosterior(std::string priorname, std::vector<double>& values, int& minprob, float& mean, float& likelihood);
      double CalculatePosterior(std::string priorname, std::vector<double>& values);


      bool isProbabilityGood(float& old_prob, float& new_prob){
        return (old_prob-new_prob) < fProbSeedCut;
      }

      bool isPosteriorProbabilityGood(double& prob, double& old_posteior){
        return (old_posteior - prob) < fPostiorCut;
      }

      bool CheckPoint(std::string priorname, double& value);

      std::vector<double> GetLikelihooddEdxVec(double& electronprob, double& photonprob,
          std::string prior,
          std::vector<double>& dEdxVec
          );

      std::vector<double> MakeSeed(std::vector<double>& dEdxVec);

      void ForceSeedToFit(std::vector<double>& SeedTrack, std::string& prior, float& mean,double& posterior);

      void RecurivelyAddHit(std::vector<double>& SeedTrack,
          std::vector<double>& dEdxVec,
          std::string& prior,
          int& SkippedHitsNum,
          float& old_mean,
          double& old_posteior
          );

      TH1F* electronpriorHist;
      TH1F* photonpriorHist;

      //fcl params
      int fVerbose;
      std::string fdEdxInputLabel;
      int fNumSeedHits;
      float fProbSeedCut;
      float fProbPointCut;
      float fPostiorCut;
      int fnSkipHits;
      std::string fShowerdEdxOutputLabel;
      bool fDefineBestPlane;
      std::string fShowerBestPlaneOutputLabel;
  };


  ShowerBayesianTrucatingdEdx::ShowerBayesianTrucatingdEdx(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fVerbose(pset.get<int>("Verbose")),
    fdEdxInputLabel(pset.get<std::string>("dEdxInputLabel")),
    fNumSeedHits(pset.get<int>("NumSeedHits")),
    fProbSeedCut(pset.get<float>("ProbDiff")),
    fProbPointCut(pset.get<float>("ProbDiffSeed")),
    fPostiorCut(pset.get<float>("PostiorCut")),
    fnSkipHits(pset.get<int>("nSkipHits")),
    fShowerdEdxOutputLabel(pset.get<std::string>("ShowerdEdxOutputLabel")),
    fDefineBestPlane(pset.get<bool>("DefineBestPlane")),
    fShowerBestPlaneOutputLabel(pset.get<std::string>("ShowerBestPlaneOutputLabel"))
  {

    //Get the prior file name
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    auto PriorPath = pset.get<std::string>("PriorFname");
    if (!sp.find_file(PriorPath, fname)) {
      throw cet::exception("ShowerBayesianTrucatingdEdx") << "Could not find the prior file";
    }
    std::string electron_histoname = pset.get<std::string>("PriorElectronHistoName");
    std::string photon_histoname = pset.get<std::string>("PriorPhotonHistoName");

    TFile fin(fname.c_str(), "READ");
    if (!fin.IsOpen()) {
      throw cet::exception("ShowerBayesianTrucatingdEdx") << "Could read the prior file. Stopping";
    }

    //Get the histograms.
    electronpriorHist = dynamic_cast<TH1F*>(fin.Get(electron_histoname.c_str()));
    if (!electronpriorHist) {
      throw cet::exception("ShowerBayesianTrucatingdEdx") << "Could not read the electron hist";
    }
    photonpriorHist = dynamic_cast<TH1F*>(fin.Get(photon_histoname.c_str()));
    if (!photonpriorHist) {
      throw cet::exception("ShowerBayesianTrucatingdEdx") << "Could not read the photon hist ";
    }

    if(electronpriorHist->GetNbinsX() != photonpriorHist->GetNbinsX()){
      throw cet::exception("ShowerBayesianTrucatingdEdx") << "Histrogram bins do not match";
    }


    //Normalise the histograms.
    electronpriorHist->Scale(1/electronpriorHist->Integral());
    photonpriorHist->Scale(1/photonpriorHist->Integral());

  }

  ShowerBayesianTrucatingdEdx::~ShowerBayesianTrucatingdEdx()
  {
  }

  int ShowerBayesianTrucatingdEdx::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event, reco::shower::ShowerElementHolder& ShowerEleHolder){

    //The idea , to some peoples distaste, is to attempt to improve the dEdx value by assuming
    //The particle is either a) electron b) a e-e+ pair.
    //We will take the start of track and work down until a few hits destory our postier probability.

    //Note: tried takeing the postior with the highest sum of the probabilitys on all three
    //      planes and on the 2 planes with the most hits. Made things worse.

    //Get the vectors of the dEdx Elements
    if(!ShowerEleHolder.CheckElement(fdEdxInputLabel)){
      fVerbose(pset.get<int>("Verbose")),
      mf::LogError("ShowerSlidingStandardCalodEdx") << "Start position not set, returning "<< std::endl;
      return 1;
    }



    std::map<int,std::vector<double > > dEdx_plane_final;
    std::map<int,std::vector<double > > dEdx_vec_planes;
    ShowerEleHolder.GetElement(fdEdxInputLabel,dEdx_vec_planes);

    //Do this for each plane;
    for(auto const& dEdx_vec_plane: dEdx_vec_planes){


      //Set up out final value if we don't have any points.
      if(dEdx_vec_plane.second.size() < 1){
        dEdx_plane_final[dEdx_vec_plane.first] = {};
        continue;
      }

      std::vector<double> dEdx_vec = dEdx_vec_plane.second;

      double electronprob_eprior = 0;
      double photonprob_eprior   = 0;

      double electronprob_pprior = 0;
      double photonprob_pprior   = 0;

      std::vector<double> dEdx_electronprior = GetLikelihooddEdxVec(electronprob_eprior,photonprob_eprior,"electron",dEdx_vec);
      std::vector<double> dEdx_photonprior   = GetLikelihooddEdxVec(electronprob_pprior,photonprob_pprior,"photon",dEdx_vec);


      //Use the vector which maximises both priors.
      if(electronprob_eprior < photonprob_pprior){
        dEdx_plane_final[dEdx_vec_plane.first] = dEdx_photonprior;
      }
      else{
        dEdx_plane_final[dEdx_vec_plane.first] = dEdx_electronprior;
      }

    }//Plane Loop

    //Calculate the median of the of dEdx.
    std::vector<double> dEdx_final;
    std::vector<double> dEdx_finalErr;

    int max_hits   = -999;
    int best_plane = -999;

    bool check = false;
    for(auto const& dEdx_plane: dEdx_plane_final){

      //Redefine the best plane
      if((int) (dEdx_plane.second).size() > max_hits){
        best_plane = dEdx_plane.first;
        max_hits   = (dEdx_plane.second).size();
      }

      if((dEdx_plane.second).size() == 0){
        dEdx_final.push_back(-999);
        dEdx_finalErr.push_back(-999);
        continue;
      }

      dEdx_final.push_back(TMath::Median((dEdx_plane.second).size(), &(dEdx_plane.second)[0]));
      dEdx_finalErr.push_back(-999);
      check = true;
    }

    //Check at least one plane has the information
    if(!check)
      return 1;

    if(fDefineBestPlane){
      ShowerEleHolder.SetElement(best_plane,fShowerBestPlaneOutputLabel);
    }

    ShowerEleHolder.SetElement(dEdx_final,dEdx_finalErr,fShowerdEdxOutputLabel);

    return 0;
  }

  double ShowerBayesianTrucatingdEdx::CalculatePosterior(std::string priorname, std::vector<double>& values){
    int minprob_iter = -999;
    float mean       = -999;
    float likelihood = -999;
    return CalculatePosterior(priorname,values,minprob_iter,mean,likelihood);
  }

  double ShowerBayesianTrucatingdEdx::CalculatePosterior(std::string priorname, std::vector<double>& values, int& minprob_iter, float& mean, float& likelihood){

    //Posterior prob;
    float posterior  = 1;
    float meanprob  = 0;
    float likelihood_other = 1;
    likelihood = 1;

    //Minimum probability temp
    float minprob_temp = 9999;
    minprob_iter = 0;

    TH1F* prior_hist = NULL;
    TH1F* other_hist = NULL;

    if(priorname=="electron"){prior_hist = electronpriorHist; other_hist = photonpriorHist;}
    if(priorname=="photon")  {prior_hist = photonpriorHist; other_hist = electronpriorHist;}

    TAxis *xaxis = prior_hist->GetXaxis();

    //Loop over the hits and calculate the probability
    for(int i=0; i<(int)values.size(); ++i){

      float value = values[i];

      Int_t bin = xaxis->FindBin(value);

      float prob = -9999;
      float other_prob =-9999;

      if(bin != xaxis->GetNbins() || bin == 0){
        //Calculate the likelihood
        prob = prior_hist->GetBinContent(bin);
        other_prob = other_hist->GetBinContent(bin);
      }
      else{
        prob = 0;
        other_prob = 0;
      }

      if(prob < minprob_temp){
        minprob_temp = prob;
        minprob_iter = i;
      }

      if(prob == 0 && other_prob == 0){continue;}

      //Calculate the posterior the mean probability and liklihood
      meanprob   += prior_hist->GetBinContent(bin);
      likelihood *= prob;
      likelihood_other *= other_prob;
    }

    posterior = likelihood/(likelihood+likelihood_other);

    meanprob /= values.size();
    mean = meanprob;
    return posterior;
  }

  bool ShowerBayesianTrucatingdEdx::CheckPoint(std::string priorname, double& value){

    TH1F* prior_hist = NULL;

    if(priorname=="electron"){prior_hist = electronpriorHist;}
    if(priorname=="photon")  {prior_hist = photonpriorHist;}

    TAxis *xaxis = prior_hist->GetXaxis();

    Int_t bin = xaxis->FindBin(value);

    float prob = -9999;

    if(bin != xaxis->GetNbins()+1 || bin == 0){
      //Calculate the likelihood
      prob = prior_hist->GetBinContent(bin);
    }
    else{
      prob = 0;
    }

    //Return the probability of getting that point.
    return prob > fProbPointCut;
  }


  std::vector<double> ShowerBayesianTrucatingdEdx::GetLikelihooddEdxVec(double& electronprob, double& photonprob,std::string prior,std::vector<double>& dEdxVec){

    //have a pool
    std::vector<double> dEdxVec_temp = dEdxVec;

    //Get The seed track.
    std::vector<double> SeedTrack = MakeSeed(dEdxVec_temp);
    //    if(SeedTrack.size() < 1){
    //  return SeedTrack;
    // }

    //Force the seed the be a good likelihood.
    float mean = 999;
    double posterior = 999;
    ForceSeedToFit(SeedTrack,prior,mean,posterior);

    //Recursively add dEdx
    int SkippedHitsNum = 0;
    RecurivelyAddHit(SeedTrack,dEdxVec_temp,prior,SkippedHitsNum,mean,posterior);

    //Calculate the likelihood of the vector  with the photon and electron priors.
    electronprob = CalculatePosterior("electron",SeedTrack);
    photonprob   = CalculatePosterior("photon",SeedTrack);

    return SeedTrack;

  }

  std::vector<double> ShowerBayesianTrucatingdEdx::MakeSeed(std::vector<double>& dEdxVec){

    std::vector<double> seed_vector;

    //Add the first hits to the seed
    int MaxHit = fNumSeedHits;
    if(fNumSeedHits > (int) dEdxVec.size()){MaxHit = (int) dEdxVec.size();}

    //    if(MaxHit == 0){
    //  mf::LogError("ShowerBayesianTrucatingdEdx") << "Size of the vector is 0 cannot perform the dEdx cutting "<< std::endl;
    //}

    for(int hit_iter=0; hit_iter<MaxHit; ++hit_iter){
      seed_vector.push_back(dEdxVec[0]);
      dEdxVec.erase(dEdxVec.begin());
    }

    return seed_vector;
  }

  void ShowerBayesianTrucatingdEdx::ForceSeedToFit(std::vector<double>& SeedTrack, std::string& prior, float& mean, double& posterior){

    int minprob_iter = 999;
    float likelihood = -999;
    float prob = CalculatePosterior(prior,SeedTrack,minprob_iter,mean,likelihood);
    while((mean < fProbSeedCut || prob <= 0) && SeedTrack.size() > 1){

      //Remove the the worse point.
      // std::cout << "removing hit with dEdx: " << SeedTrack.at(minprob_iter) << std::endl;
      SeedTrack.erase(SeedTrack.begin() + minprob_iter);
      minprob_iter = 999;

      //Recalculate
      prob = CalculatePosterior(prior,SeedTrack,minprob_iter,mean,likelihood);
    }
    posterior = prob;
    // std::cout << "seed has been fit with size: " << SeedTrack.size() << std::endl;
    return;
  }

  void ShowerBayesianTrucatingdEdx::RecurivelyAddHit(std::vector<double>& SeedTrack, std::vector<double>& dEdxVec, std::string& prior, int& SkippedHitsNum, float& old_mean, double& old_posteior){

    //If we have no more hits to add then lets finish.
    if(dEdxVec.size() < 1){return;}


    bool ok = CheckPoint(prior,dEdxVec[0]);
    // int minprob_iter=999;
    // float mean = -999;
    // float likelihood =999;
    // if(ok){std::cout << "passed first cut" << std::endl;}
    // else{std::cout << "failed first cut" << std::endl;}
    // double prob = CalculatePosterior(prior,SeedTrack,minprob_iter,mean,likelihood);
    // ok *= isProbabilityGood(mean,old_mean);
    // if(ok){std::cout << "passed second cut" << std::endl;}
    // else{std::cout << "failed second cut" << std::endl;}
    // ok *= isPosteriorProbabilityGood(prob,old_posteior);
    // if(ok){std::cout << "passed this cut" << std::endl;}
    // else{std::cout << "failed third cut" << std::endl;}


    //If we failed lets try the next hits
    if(!ok){
      // std::cout << "failed the pass point: " << dEdxVec[0] << " trying another hit" << SkippedHitsNum << std::endl;
      //if(SeedTrack.size() > 1){SeedTrack.pop_back();}
      ++SkippedHitsNum;
      if(SkippedHitsNum > fnSkipHits){return;}
    }
    else{
      //Add the next point in question.
      // std::cout << "adding value: " << dEdxVec[0] << std::endl;
      //Reset the skip number
      SkippedHitsNum = 0;
      SeedTrack.push_back(dEdxVec[0]);
    }

    //We have analysed this hit now lets remove it.
    dEdxVec.erase(dEdxVec.begin());

    RecurivelyAddHit(SeedTrack,dEdxVec,prior,SkippedHitsNum,old_mean,old_posteior);

    return;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerBayesianTrucatingdEdx)
