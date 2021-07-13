#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerCheatingAlg.h"

shower::LArPandoraShowerCheatingAlg::LArPandoraShowerCheatingAlg(const fhicl::ParameterSet& pset)
  : fLArPandoraShowerAlg(pset.get<fhicl::ParameterSet>("LArPandoraShowerAlg"))
  , fHitModuleLabel(pset.get<art::InputTag>("HitModuleLabel"))
  , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
  , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
  , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  , fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel"))
{}

std::map<int, const simb::MCParticle*>
shower::LArPandoraShowerCheatingAlg::GetTrueParticleMap() const
{

  const sim::ParticleList& particles = particleInventory->ParticleList();

  std::map<int, const simb::MCParticle*> trueParticles;
  // Make a map of track id to particle
  for (sim::ParticleList::const_iterator particleIt = particles.begin();
       particleIt != particles.end();
       ++particleIt) {
    const simb::MCParticle* particle = particleIt->second;
    trueParticles[particle->TrackId()] = particle;
    //std::cout<<"Particle ID: "<<particle->TrackId()<<" and PDG: "<<particle->PdgCode()<<std::endl;
  }
  return trueParticles;
}

std::map<int, std::vector<int>>
shower::LArPandoraShowerCheatingAlg::GetTrueChain(
  std::map<int, const simb::MCParticle*>& trueParticles) const
{

  // Roll up showers if not already done:
  std::map<int, std::vector<int>> showerMothers;

  // Loop over daughters and find th`e mothers
  for (const auto& particleIt : trueParticles) {
    const simb::MCParticle* particle = particleIt.second;
    const simb::MCParticle* mother = particle;

    if (std::abs(particle->PdgCode()) != 11 && std::abs(particle->PdgCode()) != 22) { continue; }

    // While the grand mother exists and is an electron or photon
    // Note the true mother will skip this loop and fill itself into the map
    while (mother->Mother() != 0 && trueParticles.find(mother->Mother()) != trueParticles.end()) {

      int motherId = mother->Mother();
      if (std::abs(trueParticles[motherId]->PdgCode()) != 11 &&
          std::abs(trueParticles[motherId]->PdgCode()) != 22) {
        break;
      }
      mother = trueParticles[motherId];
    }
    showerMothers[mother->TrackId()].push_back(particle->TrackId());
  }
  return showerMothers;
}

void
shower::LArPandoraShowerCheatingAlg::CheatDebugEVD(
  detinfo::DetectorClocksData const& clockData,
  const simb::MCParticle* trueParticle,
  art::Event const& Event,
  reco::shower::ShowerElementHolder& ShowerEleHolder,
  const art::Ptr<recob::PFParticle>& pfparticle) const
{

  std::cout << "Making Debug Event Display" << std::endl;

  //Function for drawing reco showers to check direction and initial track selection

  // Get run info to make unique canvas names
  int run = Event.run();
  int subRun = Event.subRun();
  int event = Event.event();
  int PFPID = pfparticle->Self();

  // Create the canvas
  TString canvasName = Form("canvas_%i_%i_%i_%i", run, subRun, event, PFPID);
  TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);

  // Initialise variables
  double x = 0;
  double y = 0;
  double z = 0;

  std::vector<art::Ptr<recob::SpacePoint>> showerSpacePoints;
  std::vector<art::Ptr<recob::SpacePoint>> otherSpacePoints;

  auto const hitHandle = Event.getValidHandle<std::vector<recob::Hit>>(fHitModuleLabel);
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits, hitHandle);

  // Get the hits associated with the space points
  const art::FindManyP<recob::SpacePoint> fmsph(hitHandle, Event, fPFParticleLabel);
  if (!fmsph.isValid()) {
    throw cet::exception("LArPandoraShowerCheatingAlg")
      << "Spacepoint and hit association not valid. Stopping.";
  }

  std::map<art::Ptr<recob::SpacePoint>, art::Ptr<recob::Hit>> spacePointHitMap;
  //Get the hits from the true particle
  for (auto hit : hits) {
    int trueParticleID = std::abs(TrueParticleID(clockData, hit));
    std::vector<art::Ptr<recob::SpacePoint>> sps = fmsph.at(hit.key());
    if (sps.size() == 1) {
      art::Ptr<recob::SpacePoint> sp = sps.front();
      if (trueParticleID == trueParticle->TrackId()) { showerSpacePoints.push_back(sp); }
      else {
        otherSpacePoints.push_back(sp);
      }
    }
  }

  if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
    mf::LogError("LArPandoraShowerCheatingAlg")
      << "Start position not set, returning " << std::endl;
    return;
  }
  if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel)) {
    mf::LogError("LArPandoraShowerCheatingAlg") << "Direction not set, returning " << std::endl;
    return;
  }
  if (!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)) {
    mf::LogError("LArPandoraShowerCheatingAlg")
      << "TrackSpacePoints not set, returning " << std::endl;
    return;
  }

  // Get info from shower property holder
  TVector3 showerStartPosition = {-999, -999, -999};
  TVector3 showerDirection = {-999, -999, -999};
  std::vector<art::Ptr<recob::SpacePoint>> trackSpacePoints;

  ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, showerStartPosition);
  ShowerEleHolder.GetElement(fShowerDirectionInputLabel, showerDirection);
  ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel, trackSpacePoints);

  // Create 3D point at vertex, chosed to be origin for ease of use of display
  double startXYZ[3] = {0, 0, 0};
  auto startPoly = std::make_unique<TPolyMarker3D>(1, startXYZ);

  // Get the min and max projections along the direction to know how long to draw
  // the direction line
  double minProj = std::numeric_limits<double>::max();
  double maxProj = -std::numeric_limits<double>::max();

  double x_min = std::numeric_limits<double>::max(), x_max = -std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max(), y_max = -std::numeric_limits<double>::max();
  double z_min = std::numeric_limits<double>::max(), z_max = -std::numeric_limits<double>::max();

  //initialise counter point
  int point = 0;

  // Make 3D points for each spacepoint in the shower
  auto showerPoly = std::make_unique<TPolyMarker3D>(showerSpacePoints.size());
  for (auto spacePoint : showerSpacePoints) {
    TVector3 pos = fLArPandoraShowerAlg.SpacePointPosition(spacePoint) - showerStartPosition;

    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    x_min = std::min(x, x_min);
    x_max = std::max(x, x_max);
    y_min = std::min(y, y_min);
    y_max = std::max(y, y_max);
    z_min = std::min(z, z_min);
    z_max = std::max(z, z_max);

    showerPoly->SetPoint(point, x, y, z);
    ++point;

    // Calculate the projection of (point-startpoint) along the direction
    double proj =
      fLArPandoraShowerAlg.SpacePointProjection(spacePoint, showerStartPosition, showerDirection);

    maxProj = std::max(proj, maxProj);
    minProj = std::min(proj, minProj);
  } // loop over spacepoints

  // Create TPolyLine3D arrays
  double xDirPoints[2] = {minProj * showerDirection.X(), maxProj * showerDirection.X()};
  double yDirPoints[2] = {minProj * showerDirection.Y(), maxProj * showerDirection.Y()};
  double zDirPoints[2] = {minProj * showerDirection.Z(), maxProj * showerDirection.Z()};

  auto dirPoly = std::make_unique<TPolyLine3D>(2, xDirPoints, yDirPoints, zDirPoints);

  point = 0; // re-initialise counter
  auto trackPoly = std::make_unique<TPolyMarker3D>(trackSpacePoints.size());
  for (auto spacePoint : trackSpacePoints) {
    TVector3 pos = fLArPandoraShowerAlg.SpacePointPosition(spacePoint) - showerStartPosition;
    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    trackPoly->SetPoint(point, x, y, z);
    ++point;
  } // loop over track spacepoints

  //  we want to draw all of the PFParticles in the event
  //Get the PFParticles

  auto otherPoly = std::make_unique<TPolyMarker3D>(otherSpacePoints.size());

  // initialise counters
  point = 0; // re-initialise counter

  for (auto const& sp : otherSpacePoints) {
    TVector3 pos = fLArPandoraShowerAlg.SpacePointPosition(sp) - showerStartPosition;
    x = pos.X();
    y = pos.Y();
    z = pos.Z();
    x_min = std::min(x, x_min);
    x_max = std::max(x, x_max);
    y_min = std::min(y, y_min);
    y_max = std::max(y, y_max);
    z_min = std::min(z, z_min);
    z_max = std::max(z, z_max);
    otherPoly->SetPoint(point, x, y, z);
    ++point;
  }

  gStyle->SetOptStat(0);
  TH3F axes("axes", "", 1, x_min, x_max, 1, y_min, y_max, 1, z_min, z_max);
  axes.SetDirectory(0);
  axes.GetXaxis()->SetTitle("X");
  axes.GetYaxis()->SetTitle("Y");
  axes.GetZaxis()->SetTitle("Z");
  axes.Draw();

  otherPoly->SetMarkerStyle(20);
  otherPoly->SetMarkerColor(4);
  otherPoly->Draw();

  // Draw all of the things
  showerPoly->SetMarkerStyle(20);
  showerPoly->Draw();
  trackPoly->SetMarkerStyle(20);
  trackPoly->SetMarkerColor(2);
  trackPoly->Draw();
  startPoly->SetMarkerStyle(21);
  startPoly->SetMarkerSize(2);
  startPoly->SetMarkerColor(3);
  startPoly->Draw();
  dirPoly->SetLineWidth(1);
  dirPoly->SetLineColor(6);
  dirPoly->Draw();

  // Save the canvas. Don't usually need this when using TFileService but this in the alg
  // not a module and didn't work without this so im going with it.
  canvas->Write();
}

int
shower::LArPandoraShowerCheatingAlg::TrueParticleID(detinfo::DetectorClocksData const& clockData,
                                                    const art::Ptr<recob::Hit>& hit) const
{

  double particleEnergy = 0;
  int likelyTrackID = 0;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
  for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
    if (trackIDs.at(idIt).energy > particleEnergy) {
      particleEnergy = trackIDs.at(idIt).energy;
      likelyTrackID = trackIDs.at(idIt).trackID;
    }
  }
  return likelyTrackID;
}

std::pair<int, double>
shower::LArPandoraShowerCheatingAlg::TrueParticleIDFromTrueChain(
  detinfo::DetectorClocksData const& clockData,
  std::map<int, std::vector<int>> const& ShowersMothers,
  std::vector<art::Ptr<recob::Hit>> const& hits,
  int planeid) const
{

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;

  //Find the energy for each track ID.
  std::map<int, double> trackIDToEDepMap;
  std::map<int, double> trackIDTo3EDepMap;
  for (std::vector<art::Ptr<recob::Hit>>::const_iterator hitIt = hits.begin(); hitIt != hits.end();
       ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;

    //Get the plane ID
    geo::WireID wireid = (*hitIt)->WireID();
    int PlaneID = wireid.Plane;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      trackIDTo3EDepMap[std::abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
      if (PlaneID == planeid) {
        trackIDToEDepMap[std::abs(trackIDs[idIt].trackID)] += trackIDs[idIt].energy;
      }
    }
  }

  //Find the energy for each showermother.
  std::map<int, double> MotherIDtoEMap;
  std::map<int, double> MotherIDto3EMap;
  for (std::map<int, std::vector<int>>::const_iterator showermother = ShowersMothers.begin();
       showermother != ShowersMothers.end();
       ++showermother) {
    for (std::vector<int>::const_iterator daughter = (showermother->second).begin();
         daughter != (showermother->second).end();
         ++daughter) {
      MotherIDtoEMap[showermother->first] += trackIDToEDepMap[*daughter];
      MotherIDto3EMap[showermother->first] += trackIDTo3EDepMap[*daughter];
    }
  }

  //Loop over the mothers to find the most like candidate by identifying which true shower deposited the most energy in the hits.
  double maxenergy = -1;
  int objectTrack = -99999;
  for (std::map<int, double>::iterator mapIt = MotherIDto3EMap.begin();
       mapIt != MotherIDto3EMap.end();
       mapIt++) {
    double energy_three = mapIt->second;
    double trackid = mapIt->first;
    if (energy_three > maxenergy) {
      maxenergy = energy_three;
      objectTrack = trackid;
    }
  }

  //If the none of the shower mother deposited no energy then we cannot match this.
  if (maxenergy == 0) { return std::make_pair(-99999, -99999); }

  return std::make_pair(objectTrack, MotherIDtoEMap[objectTrack]);
}
