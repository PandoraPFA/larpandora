//############################################################################
//### Name:        ShowerPCADirection                                      ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the shower direction using PCA         ###
//###              methods. Derived from LArPandoraModularShowers Method.            ###
//############################################################################

//Warning! Currently as pandora gives each hit a spacepoint, rather than
//         matching up some energy depositions are double counted.
//         This could lead to a bais in the PCA analysis.

#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

//Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/FindManyP.h"

//LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/Shower.h"

//C++ Includes
#include <iostream>
#include <Eigen/Dense>

//Root Includes
#include "TVector3.h"
#include "TMath.h"

namespace ShowerRecoTools {

  class ShowerPCADirection: public IShowerTool {

    public:

      ShowerPCADirection(const fhicl::ParameterSet& pset);

      ~ShowerPCADirection();

      //Calculate the direction of the shower.
      int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
          art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder
          ) override;

    private:

      void InitialiseProducers() override;

      //Function to add the assoctions
      int AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
          reco::shower::ShowerElementHolder& ShowerEleHolder) override;

      // Define standard art tool interface
      recob::PCAxis CalculateShowerPCA(std::vector<art::Ptr<recob::SpacePoint> >& spacePoints_pfp,
          art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre);

      TVector3 GetPCAxisVector(recob::PCAxis& PCAxis);

      double  RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps,
          TVector3& ShowerCentre,
          TVector3& Direction);

      // Function to calculate the RMS spread of perpendicular distances from PCA
      double CalculateRMS(const std::vector<double>& perpVec);

      //Services
      detinfo::DetectorProperties const* fDetProp;

      //fcl
      art::InputTag fPFParticleModuleLabel;
      float fNSegments;        //Used in the RMS gradient. How many segments should we split the shower into.
      bool fUseStartPosition;  //If we use the start position the drection of the
      //PCA vector is decided as (Shower Centre - Shower Start Position).
      bool fChargeWeighted;    //Should the PCA axis be charge weighted.

      std::string fShowerStartPositionInputLabel;
      std::string fShowerDirectionOutputLabel;
      std::string fShowerCentreOutputLabel;
      std::string fShowerPCAOutputLabel;

  };

  ShowerPCADirection::ShowerPCADirection(const fhicl::ParameterSet& pset) :
    IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools")),
    fDetProp(lar::providerFrom<detinfo::DetectorPropertiesService>()),
    fPFParticleModuleLabel(pset.get<art::InputTag>("PFParticleModuleLabel")),
    fNSegments(pset.get<float>("NSegments")),
    fUseStartPosition(pset.get<bool>("UseStartPosition")),
    fChargeWeighted(pset.get<bool>("ChargeWeighted")),
    fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel")),
    fShowerDirectionOutputLabel(pset.get<std::string>("ShowerDirectionOutputLabel")),
    fShowerCentreOutputLabel(pset.get<std::string>("ShowerCentreOutputLabel")),
    fShowerPCAOutputLabel(pset.get<std::string>("ShowerPCAOutputLabel"))
  {
  }

  ShowerPCADirection::~ShowerPCADirection()
  {
  }

  void ShowerPCADirection::InitialiseProducers()
  {
    InitialiseProduct<std::vector<recob::PCAxis> >(fShowerPCAOutputLabel);
    InitialiseProduct<art::Assns<recob::Shower, recob::PCAxis> >("ShowerPCAxisAssn");
    InitialiseProduct<art::Assns<recob::PFParticle, recob::PCAxis> >("PFParticlePCAxisAssn");
  }

  int ShowerPCADirection::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
      art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder){

    // Get the assocated pfParicle vertex PFParticles
    art::Handle<std::vector<recob::PFParticle> > pfpHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, pfpHandle)){
      throw cet::exception("ShowerPCADirection") << "Could not get the pandora pf particles. \
        Something is not cofingured coreectly Please give the correct pandoa module label. Stopping";
      return 1;
    }

    art::FindManyP<recob::SpacePoint>& fmspp = ShowerEleHolder.GetFindManyP<recob::SpacePoint>(
        pfpHandle, Event, fPFParticleModuleLabel);
    // art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleModuleLabel);

    if (!fmspp.isValid()){
      throw cet::exception("ShowerPCADirection") << "Trying to get the spacepoint and failed. Something is not configured correctly. Stopping ";
      return 1;
    }

    //Get the spacepoints handle and the hit assoication
    art::Handle<std::vector<recob::SpacePoint> > spHandle;
    if (!Event.getByLabel(fPFParticleModuleLabel, spHandle)){
      throw cet::exception("ShowerPCADirection") << "Could not configure the spacepoint handle. Something is configured incorrectly. Stopping";
      return 1;
    }
    art::FindManyP<recob::Hit>& fmh = ShowerEleHolder.GetFindManyP<recob::Hit>(
        spHandle, Event, fPFParticleModuleLabel);
    // art::FindManyP<recob::Hit> fmh(spHandle, Event, fPFParticleModuleLabel);
    if(!fmh.isValid()){
      throw cet::exception("ShowerPCADirection") << "Spacepoint and hit association not valid. Stopping.";
      return 1;
    }

    //Spacepoints
    std::vector<art::Ptr<recob::SpacePoint> > spacePoints_pfp = fmspp.at(pfparticle.key());

    //We cannot progress with no spacepoints.
    if(spacePoints_pfp.size() == 0){return 0;}

    //Find the PCA vector
    TVector3 ShowerCentre;
    recob::PCAxis PCA = CalculateShowerPCA(spacePoints_pfp,fmh,ShowerCentre);
    TVector3 PCADirection = GetPCAxisVector(PCA);

    //Save the shower the center for downstream tools
    TVector3 ShowerCentreErr = {-999,-999,-999};
    ShowerEleHolder.SetElement(ShowerCentre,ShowerCentreErr,fShowerCentreOutputLabel);
    ShowerEleHolder.SetElement(PCA, fShowerPCAOutputLabel);

    //Check if we are pointing the correct direction or not, First try the start position
    if(fUseStartPosition){
      if(!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)){
        throw cet::exception("ShowerPCADirection") << "fUseStartPosition is true but start position is not set. Stopping.";
        return 1;
      }
      //Get the General direction as the vector between the start position and the centre
      TVector3 StartPositionVec = {-999, -999, -999};
      ShowerEleHolder.GetElement(fShowerStartPositionInputLabel,StartPositionVec);

      // Calculate the general direction of the shower
      TVector3 GeneralDir = (ShowerCentre - StartPositionVec).Unit();

      //Calculate the dot product between eigenvector and general direction
      double DotProduct = PCADirection.Dot(GeneralDir);

      //If the dotproduct is negative the Direction needs Flipping
      if(DotProduct < 0){
        PCADirection[0] = - PCADirection[0];
        PCADirection[1] = - PCADirection[1];
        PCADirection[2] = - PCADirection[2];
      }

      //To do
      TVector3 PCADirectionErr = {-999,-999,-999};
      ShowerEleHolder.SetElement(PCADirection,PCADirectionErr,fShowerDirectionOutputLabel);
      return 0;
    }

    //Otherwise Check against the RMS of the shower. Method adapated from EMShower Thanks Mike.
    double RMSGradient = RMSShowerGradient(spacePoints_pfp,ShowerCentre,PCADirection);

    if(RMSGradient < 0){
      PCADirection[0] = - PCADirection[0];
      PCADirection[1] = - PCADirection[1];
      PCADirection[2] = - PCADirection[2];
    }

    //To do
    TVector3 PCADirectionErr = {-999,-999,-999};

    ShowerEleHolder.SetElement(PCADirection,PCADirectionErr,fShowerDirectionOutputLabel);
    return 0;
  }

  //Function to calculate the RMS at segements of the shower and calculate the gradient of this. If negative then the direction is pointing the opposite way to the correct one
  double ShowerPCADirection::RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint> >& sps, TVector3& ShowerCentre, TVector3& Direction){

    //Order the spacepoints
    IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(sps,ShowerCentre,Direction);

    //Get the length of the shower.
    double minProj =IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(sps[0],ShowerCentre,Direction);
    double maxProj =IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(sps[sps.size()-1],ShowerCentre,Direction);

    double length = (maxProj-minProj);
    double segmentsize = length/fNSegments;

    // Create a map from segment (projection) to perpendicular distance
    std::map<int, std::vector<double> > len_segment_map;

    //Split the the spacepoints into segments.
    for(auto const& sp: sps){

      //Get the the projected length
      double proj = IShowerTool::GetLArPandoraShowerAlg().SpacePointProjection(sp,ShowerCentre,Direction);

      //Get the length to the projection
      double perp = IShowerTool::GetLArPandoraShowerAlg().SpacePointPerpendicular(sp,ShowerCentre,Direction,proj);

      int segment = round(proj/segmentsize);
      len_segment_map[segment].push_back(perp);
    }

    int counter = 0;
    double sumx  = 0;
    double sumy  = 0;
    double sumx2 = 0;
    double sumxy = 0;

    //Get the rms of the segments and caclulate the gradient.
    for(auto const& segment: len_segment_map){

      // Require at least 2 space points in a segment
      if (segment.second.size()<2) continue;

      double RMS = CalculateRMS(segment.second);
      //Calculate the gradient using regression
      sumx  += segment.first;
      sumy  += RMS;
      sumx2 += segment.first * segment.first;
      sumxy += RMS * segment.first;
      ++counter;
    }

    return (counter*sumxy - sumx*sumy)/(counter*sumx2 - sumx*sumx);
  }

  double ShowerPCADirection::CalculateRMS(const std::vector<double>& perpVec){
    double sum  = 0;
    for (const auto& perp : perpVec){
      sum = perp*perp;
    }
    return  TMath::Sqrt(sum/(perpVec.size()-1));
  }


  //Function to calculate the shower direction using a charge weight 3D PCA calculation.
  recob::PCAxis ShowerPCADirection::CalculateShowerPCA(std::vector<art::Ptr<recob::SpacePoint> >& sps, art::FindManyP<recob::Hit>& fmh, TVector3& ShowerCentre){

    float TotalCharge = 0;
    float sumWeights = 0;
    float xx = 0;
    float yy = 0;
    float zz = 0;
    float xy = 0;
    float xz = 0;
    float yz = 0;

    //Get the Shower Centre
    ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(sps, fmh, TotalCharge);

    //Normalise the spacepoints, charge weight and add to the PCA.
    for(auto& sp: sps){

      TVector3 sp_position = IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(sp);

      float wht = 1;

      //Normalise the spacepoint position.
      sp_position = sp_position - ShowerCentre;

      if(fChargeWeighted){

        //Get the charge.
        float Charge = IShowerTool::GetLArPandoraShowerAlg().SpacePointCharge(sp,fmh);

        //Get the time of the spacepoint
        float Time = IShowerTool::GetLArPandoraShowerAlg().SpacePointTime(sp,fmh);

        //Correct for the lifetime at the moment.
        Charge *= TMath::Exp((fDetProp->SamplingRate() * Time ) / (fDetProp->ElectronLifetime()*1e3));

        //Charge Weight
        wht *= TMath::Sqrt(Charge/TotalCharge);
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
    matrix <<
      xx, xy, xz,
      xy, yy, yz,
      xz, yz ,zz;

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
    const double eigenValues[3] = {eigenValuesVector(2), eigenValuesVector(1), eigenValuesVector(0)};
    std::vector<std::vector<double> > eigenVectors = {
      { eigenVectorsMatrix(0,2), eigenVectorsMatrix(1,2), eigenVectorsMatrix(2,2) },
      { eigenVectorsMatrix(0,1), eigenVectorsMatrix(1,1), eigenVectorsMatrix(2,1) },
      { eigenVectorsMatrix(0,0), eigenVectorsMatrix(1,0), eigenVectorsMatrix(2,0) }};
    const double avePos[3] = {ShowerCentre[0], ShowerCentre[1], ShowerCentre[2]};

    return  recob::PCAxis(svdOk, nHits, eigenValues, eigenVectors, avePos);
  }

  TVector3 ShowerPCADirection::GetPCAxisVector(recob::PCAxis& PCAxis){

    //Get the Eigenvectors.
    std::vector<double> EigenVector = PCAxis.getEigenVectors()[0];

    return TVector3(EigenVector[0], EigenVector[1],  EigenVector[2]);
  }


  int ShowerPCADirection::AddAssociations(const art::Ptr<recob::PFParticle>& pfpPtr, art::Event& Event,
      reco::shower::ShowerElementHolder& ShowerEleHolder){

    //First check the element has been set
    if(!ShowerEleHolder.CheckElement(fShowerPCAOutputLabel)){
      mf::LogError("ShowerPCADirection: Add Assns") << "PCA not set."<< std::endl;
      return 1;
    }

    int ptrSize = GetVectorPtrSize(fShowerPCAOutputLabel);

    const art::Ptr<recob::PCAxis> pcaPtr = GetProducedElementPtr<recob::PCAxis>(fShowerPCAOutputLabel, ShowerEleHolder, ptrSize-1);
    const art::Ptr<recob::Shower> showerPtr = GetProducedElementPtr<recob::Shower>("shower", ShowerEleHolder);

    AddSingle<art::Assns<recob::Shower, recob::PCAxis> >(showerPtr,pcaPtr,"ShowerPCAxisAssn");
    AddSingle<art::Assns<recob::PFParticle, recob::PCAxis> >(pfpPtr,pcaPtr,"PFParticlePCAxisAssn");

    return 0;
  }
}
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPCADirection)
