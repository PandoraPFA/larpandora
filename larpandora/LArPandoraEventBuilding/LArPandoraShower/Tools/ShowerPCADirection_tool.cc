//############################################################################
//### Name:        ShowerPCADirection                                      ###
//### Author:      Dominic Barker and Ed Tyley                             ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using PCA         ###
//###              methods. Derived from LArPandoraModularShowers Method.  ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Shower.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//C++ Includes
#include <Eigen/Dense>

namespace ShowerRecoTools {

  class ShowerPCADirection : public IShowerTool {

  public:
    ShowerPCADirection(const fhicl::ParameterSet& pset);

    //Calculate the direction of the shower.
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    void InitialiseProducers() override;

    //Function to add the assoctions
    int AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr,
                        art::Event& Event,
                        reco::shower::ShowerElementHolder& ShowerEleHolder) override;

    // Define standard art tool interface
    recob::PCAxis CalculateShowerPCA(
      const detinfo::DetectorClocksData& clockData,
      const detinfo::DetectorPropertiesData& detProp,
      const std::vector<art::Ptr<recob::SpacePoint>>& spacePoints_pfp,
      const art::FindManyP<recob::Hit>& fmh,
      TVector3& ShowerCentre);

    TVector3 GetPCAxisVector(recob::PCAxis& PCAxis);

    //fcl
    art::InputTag fPFParticleLabel;
    int fVerbose;
    unsigned int
      fNSegments; //Used in the RMS gradient. How many segments should we split the shower into.
    bool fUseStartPosition; //If we use the start position the drection of the
    //PCA vector is decided as (Shower Centre - Shower Start Position).
    bool fChargeWeighted; //Should the PCA axis be charge weighted.

    std::string fShowerStartPositionInputLabel;
    std::string fShowerDirectionOutputLabel;
    std::string fShowerCentreOutputLabel;
    std::string fShowerPCAOutputLabel;
  };

  ShowerPCADirection::ShowerPCADirection(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fNSegments(pset.get<unsigned int>("NSegments"))
    , fUseStartPosition(pset.get<bool>("UseStartPosition"))
    , fChargeWeighted(pset.get<bool>("ChargeWeighted"))
    , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
    , fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel"))
    , fShowerCentreOutputLabel(pset.get<std::string>("ShowerCentreOutputLabel"))
    , fShowerPCAOutputLabel(pset.get<std::string>("ShowerPCAOutputLabel"))
  {}

  void
  ShowerPCADirection::InitialiseProducers()
  {
    InitialiseProduct<std::vector<recob::PCAxis>>(fShowerPCAOutputLabel);
    InitialiseProduct<art::Assns<recob::Shower, recob::PCAxis>>("ShowerPCAxisAssn");
    InitialiseProduct<art::Assns<recob::PFParticle, recob::PCAxis>>("PFParticlePCAxisAssn");
  }

  int
  ShowerPCADirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                       art::Event& Event,
                                       reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    // Get the assocated pfParicle vertex PFParticles
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    const art::FindManyP<recob::SpacePoint>& fmspp =
      ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

    //Get the spacepoints handle and the hit assoication
    auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);

    const art::FindManyP<recob::Hit>& fmh =
      ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

    //Spacepoints
    std::vector<art::Ptr<recob::SpacePoint>> spacePoints_pfp = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if (spacePoints_pfp.size() < 3) {
      mf::LogWarning("ShowerPCADirection")
        << spacePoints_pfp.size() << " spacepoints in shower, not calculating direction"
        << std::endl;
      return 1;
    }

    auto const clockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
    auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

    //Find the PCA vector
    TVector3 ShowerCentre;
    recob::PCAxis PCA = CalculateShowerPCA(clockData, detProp, spacePoints_pfp, fmh, ShowerCentre);
    TVector3 PCADirection = GetPCAxisVector(PCA);

    //Save the shower the center for downstream tools
    TVector3 ShowerCentreErr = {-999, -999, -999};
    ShowerEleHolder.SetElement(ShowerCentre, ShowerCentreErr, fShowerCentreOutputLabel);
    ShowerEleHolder.SetElement(PCA, fShowerPCAOutputLabel);

    //Check if we are pointing the correct direction or not, First try the start position
    if (fUseStartPosition) {
      if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
        if (fVerbose)
          mf::LogError("ShowerPCADirection")
            << "fUseStartPosition is set but ShowerStartPosition is not set. Bailing" << std::endl;
        return 1;
      }
      //Get the General direction as the vector between the start position and the centre
      TVector3 StartPositionVec = {-999, -999, -999};
      ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, StartPositionVec);

      // Calculate the general direction of the shower
      TVector3 GeneralDir = (ShowerCentre - StartPositionVec).Unit();

      //Calculate the dot product between eigenvector and general direction
      double DotProduct = PCADirection.Dot(GeneralDir);

      //If the dotproduct is negative the Direction needs Flipping
      if (DotProduct < 0) {
        PCADirection[0] = -PCADirection[0];
        PCADirection[1] = -PCADirection[1];
        PCADirection[2] = -PCADirection[2];
      }

      //To do
      TVector3 PCADirectionErr = {-999, -999, -999};
      ShowerEleHolder.SetElement(PCADirection, PCADirectionErr, fShowerDirectionOutputLabel);
      return 0;
    }

    //Otherwise Check against the RMS of the shower. Method adapated from EMShower Thanks Mike.
    double RMSGradient = IShowerTool::GetLArPandoraShowerAlg().RMSShowerGradient(
      spacePoints_pfp, ShowerCentre, PCADirection, fNSegments);

    // If the alg fails to calculate the gradient it will return 0. In this case do nothing
    // If the gradient is negative, flip the direction of the shower
    if (RMSGradient < -std::numeric_limits<double>::epsilon()) {
      PCADirection[0] = -PCADirection[0];
      PCADirection[1] = -PCADirection[1];
      PCADirection[2] = -PCADirection[2];
    }

    //To do
    TVector3 PCADirectionErr = {-999, -999, -999};

    ShowerEleHolder.SetElement(PCADirection, PCADirectionErr, fShowerDirectionOutputLabel);
    return 0;
  }

  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  recob::PCAxis
  ShowerPCADirection::CalculateShowerPCA(const detinfo::DetectorClocksData& clockData,
                                         const detinfo::DetectorPropertiesData& detProp,
                                         const std::vector<art::Ptr<recob::SpacePoint>>& sps,
                                         const art::FindManyP<recob::Hit>& fmh,
                                         TVector3& ShowerCentre)
  {

    float TotalCharge = 0;
    float sumWeights = 0;
    float xx = 0;
    float yy = 0;
    float zz = 0;
    float xy = 0;
    float xz = 0;
    float yz = 0;

    //Get the Shower Centre
    if (fChargeWeighted) {
      ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(
        clockData, detProp, sps, fmh, TotalCharge);
    }
    else {
      ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(sps);
    }

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

      xx += sp_position.X() * sp_position.X() * wht;
      yy += sp_position.Y() * sp_position.Y() * wht;
      zz += sp_position.Z() * sp_position.Z() * wht;
      xy += sp_position.X() * sp_position.Y() * wht;
      xz += sp_position.X() * sp_position.Z() * wht;
      yz += sp_position.Y() * sp_position.Z() * wht;
      sumWeights += wht;
    }

    // Using Eigen package
    Eigen::Matrix3f matrix;

    // Construct covariance matrix
    matrix << xx, xy, xz, xy, yy, yz, xz, yz, zz;

    // Normalise from the sum of weights
    matrix /= sumWeights;

    // Run the PCA
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMatrix(matrix);

    Eigen::Vector3f eigenValuesVector = eigenMatrix.eigenvalues();
    Eigen::Matrix3f eigenVectorsMatrix = eigenMatrix.eigenvectors();

    // Put in the required form for a recob::PCAxis
    const bool svdOk = true; //TODO: Should probably think about this a bit more
    const int nHits = sps.size();
    // For some reason eigen sorts the eigenvalues from smallest to largest, reverse it
    const double eigenValues[3] = {
      eigenValuesVector(2), eigenValuesVector(1), eigenValuesVector(0)};
    std::vector<std::vector<double>> eigenVectors = {
      {eigenVectorsMatrix(0, 2), eigenVectorsMatrix(1, 2), eigenVectorsMatrix(2, 2)},
      {eigenVectorsMatrix(0, 1), eigenVectorsMatrix(1, 1), eigenVectorsMatrix(2, 1)},
      {eigenVectorsMatrix(0, 0), eigenVectorsMatrix(1, 0), eigenVectorsMatrix(2, 0)}};
    const double avePos[3] = {ShowerCentre[0], ShowerCentre[1], ShowerCentre[2]};

    return recob::PCAxis(svdOk, nHits, eigenValues, eigenVectors, avePos);
  }

  TVector3
  ShowerPCADirection::GetPCAxisVector(recob::PCAxis& PCAxis)
  {

    //Get the Eigenvectors.
    std::vector<double> EigenVector = PCAxis.getEigenVectors()[0];

    return TVector3(EigenVector[0], EigenVector[1], EigenVector[2]);
  }

  int
  ShowerPCADirection::AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr,
                                      art::Event& Event,
                                      reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    //First check the element has been set
    if (!ShowerEleHolder.CheckElement(fShowerPCAOutputLabel)) {
      if (fVerbose) mf::LogError("ShowerPCADirection: Add Assns") << "PCA not set." << std::endl;
      return 1;
    }

    int ptrSize = GetVectorPtrSize(fShowerPCAOutputLabel);

    const art::Ptr<recob::PCAxis> pcaPtr =
      GetProducedElementPtr<recob::PCAxis>(fShowerPCAOutputLabel, ShowerEleHolder, ptrSize - 1);
    const art::Ptr<recob::Shower> showerPtr =
      GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::PCAxis>>(showerPtr, pcaPtr, "ShowerPCAxisAssn");
    AddSingle<art::Assns<recob::PFParticle, recob::PCAxis>>(pfpPtr, pcaPtr, "PFParticlePCAxisAssn");

    return 0;
  }
}
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCADirection)
