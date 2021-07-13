//############################################################################
//### Name:        ShowerTrackPCADirection                                 ###
//### Author:      Dominic Barker                                          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using PCA         ###
//###              methods using the initial track hits                    ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Root Includes
#include "TPrincipal.h"

namespace ShowerRecoTools {

  class ShowerTrackPCADirection : IShowerTool {

  public:
    ShowerTrackPCADirection(const fhicl::ParameterSet& pset);

    //Generic Direction Finder
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    TVector3 ShowerPCAVector(const detinfo::DetectorClocksData& clockData,
                             const detinfo::DetectorPropertiesData& detProp,
                             std::vector<art::Ptr<recob::SpacePoint>>& spacePoints_pfp,
                             const art::FindManyP<recob::Hit>& fmh,
                             TVector3& ShowerCentre);

    //fcl
    art::InputTag fPFParticleLabel;
    art::InputTag fHitModuleLabel;
    int fVerbose;
    bool fChargeWeighted;       //Should we charge weight the PCA.
    unsigned int fMinPCAPoints; //Number of spacepoints needed to do the analysis.

    std::string fShowerStartPositionInputLabel;
    std::string fInitialTrackSpacePointsInputLabel;
    std::string fShowerDirectionOutputLabel;
  };

  ShowerTrackPCADirection::ShowerTrackPCADirection(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fChargeWeighted(pset.get<bool>("ChargeWeighted"))
    , fMinPCAPoints(pset.get<unsigned int>("MinPCAPoints"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel"))
    , fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
  {}

  int
  ShowerTrackPCADirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                            art::Event& Event,
                                            reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrackPCA") << "Start Position not set. Stopping" << std::endl;
      ;
      return 1;
    }
    if (!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)) {
      if (fVerbose)
        mf::LogError("ShowerTrackPCA") << "TrackSpacePoints not set, returning " << std::endl;
      return 1;
    }

    //Get the spacepoints handle and the hit assoication
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    const art::FindManyP<recob::Hit>& fmh =
      ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    std::vector<art::Ptr<recob::SpacePoint>> trackSpacePoints;
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel, trackSpacePoints);

    //We cannot progress with no spacepoints.
    if (trackSpacePoints.size() < fMinPCAPoints) {
      if (fVerbose)
        mf::LogError("ShowerTrackPCA") << "Not enough spacepoints for PCA, returning " << std::endl;
      return 1;
    }

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    //Find the PCA vector
    TVector3 trackCentre;
    TVector3 Eigenvector = ShowerPCAVector(clockData, detProp, trackSpacePoints, fmh, trackCentre);

    //Get the General direction as the vector between the start position and the centre
    TVector3 StartPositionVec = {-999, -999, -999};
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, StartPositionVec);
    TVector3 GeneralDir = (trackCentre - StartPositionVec).Unit();

    //Dot product
    double DotProduct = Eigenvector.Dot(GeneralDir);

    //If the dotproduct is negative the Direction needs Flipping
    if (DotProduct < 0) {
      Eigenvector[0] = -Eigenvector[0];
      Eigenvector[1] = -Eigenvector[1];
      Eigenvector[2] = -Eigenvector[2];
    }

    TVector3 EigenvectorErr = {-999, -999, -999};

    ShowerEleHolder.SetElement(Eigenvector, EigenvectorErr, fShowerDirectionOutputLabel);

    return 0;
  }

  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  TVector3
  ShowerTrackPCADirection::ShowerPCAVector(const detinfo::DetectorClocksData& clockData,
                                           const detinfo::DetectorPropertiesData& detProp,
                                           std::vector<art::Ptr<recob::SpacePoint>>& sps,
                                           const art::FindManyP<recob::Hit>& fmh,
                                           TVector3& ShowerCentre)
  {

    //Initialise the the PCA.
    TPrincipal* pca = new TPrincipal(3, "");

    float TotalCharge = 0;

    //Get the Shower Centre
    ShowerCentre =
      IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(clockData, detProp, sps, fmh, TotalCharge);

    //Normalise the spacepoints, charge weight and add to the PCA.
    for (auto& sp : sps) {

      TVector3 sp_position = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp);

      float wht = 1;

      //Normalise the spacepoint position.
      sp_position = sp_position - ShowerCentre;

      if (fChargeWeighted) {

        //Get the charge.
        float Charge = IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(sp, fmh);

        //Get the time of the spacepoint
        float Time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(sp, fmh);

        //Correct for the lifetime at the moment.
        Charge *= std::exp((sampling_rate(clockData) * Time) / (detProp.ElectronLifetime() * 1e3));

        //Charge Weight
        wht *= std::sqrt(Charge / TotalCharge);
      }

      double sp_coord[3];
      sp_coord[0] = sp_position.X() * wht;
      sp_coord[1] = sp_position.Y() * wht;
      sp_coord[2] = sp_position.Z() * wht;

      //Add to the PCA
      pca->AddRow(sp_coord);
    }

    //Evaluate the PCA
    pca->MakePrincipals();

    //Get the Eigenvectors.
    const TMatrixD* Eigenvectors = pca->GetEigenVectors();

    TVector3 Eigenvector = {(*Eigenvectors)[0][0], (*Eigenvectors)[1][0], (*Eigenvectors)[2][0]};

    delete pca;

    return Eigenvector;
  }
}

DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerTrackPCADirection)
