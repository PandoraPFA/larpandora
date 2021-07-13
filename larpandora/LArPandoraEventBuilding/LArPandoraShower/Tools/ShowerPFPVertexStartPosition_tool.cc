//############################################################################
//### Name:        ShowerPFPVertexStartPosition                            ###
//### Author:      Dominic Barker (dominic.barker@sheffield.ac.uk          ###
//### Date:        13.05.19                                                ###
//### Description: Tool for finding the start poistion                     ###
//###              methods.                                                ###
//############################################################################

//Framework Includes
#include "art/Utilities/ToolMacros.h"

//LArSoft Includes
#include "lardataobj/RecoBase/Vertex.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Tools/IShowerTool.h"

namespace ShowerRecoTools {

  class ShowerPFPVertexStartPosition : public IShowerTool {

  public:
    ShowerPFPVertexStartPosition(const fhicl::ParameterSet& pset);

    //Calculate the start position
    int CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                         art::Event& Event,
                         reco::shower::ShowerElementHolder& ShowerEleHolder) override;

  private:
    //fcl parameters
    art::InputTag fPFParticleLabel;
    int fVerbose;
    std::string fShowerStartPositionOutputLabel;
    std::string fShowerDirectionInputLabel;
  };

  ShowerPFPVertexStartPosition::ShowerPFPVertexStartPosition(const fhicl::ParameterSet& pset)
    : IShowerTool(pset.get<fhicl::ParameterSet>("BaseTools"))
    , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
    , fVerbose(pset.get<int>("Verbose"))
    , fShowerStartPositionOutputLabel(pset.get<std::string>("ShowerStartPositionOutputLabel"))
    , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  {}

  int
  ShowerPFPVertexStartPosition::CalculateElement(const art::Ptr<recob::PFParticle>& pfparticle,
                                                 art::Event& Event,
                                                 reco::shower::ShowerElementHolder& ShowerEleHolder)
  {

    // Get the assocated pfParicle vertex PFParticles
    auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

    const art::FindManyP<recob::Vertex>& fmv =
      ShowerEleHolder.GetFindManyP<recob::Vertex>(pfpHandle, Event, fPFParticleLabel);
    // art::FindManyP<recob::Vertex> fmv(pfpHandle, Event, fPFParticleLabel);
    if (!fmv.isValid()) {
      throw cet::exception("ShowerPFPVertexStartPosition")
        << "Vertex and PF particle association is somehow not valid. Stopping";
    }

    std::vector<art::Ptr<recob::Vertex>> vtx_cand;
    try {
      vtx_cand = fmv.at(pfparticle.key());
    }
    catch (...) {
      if (fVerbose)
        mf::LogError("ShowerPFPVertexStartPosition") << "PFP-Vertex assan not set, returning";
      return 1;
    }
    //If there is more than one then fail becuase I don't think that this can be the case
    if (vtx_cand.size() != 1) {
      if (fVerbose)
        mf::LogError("ShowerPFPVertexStartPosition")
          << "Wrong number of vertices: " << vtx_cand.size() << ", returning";
      return 1;
    }

    //If there is only one vertex good news we just say that is the start of the shower.
    if (vtx_cand.size() == 1) {
      art::Ptr<recob::Vertex> StartPositionVertex = vtx_cand[0];
      double xyz[3] = {-999, -999, -999};
      StartPositionVertex->XYZ(xyz);
      TVector3 ShowerStartPosition = {xyz[0], xyz[1], xyz[2]};
      TVector3 ShowerStartPositionErr = {-999, -999, -999};
      ShowerEleHolder.SetElement(
        ShowerStartPosition, ShowerStartPositionErr, fShowerStartPositionOutputLabel);
      return 0;
    }

    //If we there have none then use the direction to find the neutrino vertex
    if (ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {

      TVector3 ShowerDirection = {-999, -999, -999};
      ShowerEleHolder.GetElement(fShowerDirectionInputLabel, ShowerDirection);

      const art::FindManyP<recob::SpacePoint>& fmspp =
        ShowerEleHolder.GetFindManyP<recob::SpacePoint>(pfpHandle, Event, fPFParticleLabel);

      //Get the spacepoints handle and the hit assoication
      auto const spHandle = Event.getValidHandle<std::vector<recob::SpacePoint>>(fPFParticleLabel);
      const art::FindManyP<recob::Hit>& fmh =
        ShowerEleHolder.GetFindManyP<recob::Hit>(spHandle, Event, fPFParticleLabel);

      //Get the spacepoints
      std::vector<art::Ptr<recob::SpacePoint>> spacePoints_pfp = fmspp.at(pfparticle.key());

      //Cannot continue if we have no spacepoints
      if (spacePoints_pfp.empty()) { return 0; }

      //Get the Shower Center
      auto const clockData =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(Event);
      auto const detProp =
        art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(Event, clockData);

      TVector3 ShowerCentre = IShowerTool::GetLArPandoraShowerAlg().ShowerCentre(
        clockData, detProp, spacePoints_pfp, fmh);

      //Order the Hits from the shower centre. The most negative will be the start position.
      IShowerTool::GetLArPandoraShowerAlg().OrderShowerSpacePoints(
        spacePoints_pfp, ShowerCentre, ShowerDirection);

      //Set the start position.
      TVector3 ShowerStartPosition =
        IShowerTool::GetLArPandoraShowerAlg().SpacePointPosition(spacePoints_pfp[0]);

      TVector3 ShowerStartPositionErr = {-999, -999, -999};
      ShowerEleHolder.SetElement(
        ShowerStartPosition, ShowerStartPositionErr, fShowerStartPositionOutputLabel);

      return 0;
    }

    if (fVerbose)
      mf::LogWarning("ShowerPFPVertexStartPosition")
        << "Start Position has not been set yet. If you are not calculating the start position "
           "again then maybe you should stop";
    return 0;
  }

}
DEFINE_ART_CLASS_TOOL(ShowerRecoTools::ShowerPFPVertexStartPosition)
