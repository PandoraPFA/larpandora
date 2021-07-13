#include "larpandora/LArPandoraEventBuilding/LArPandoraShower/Algs/LArPandoraShowerAlg.h"

shower::LArPandoraShowerAlg::LArPandoraShowerAlg(const fhicl::ParameterSet& pset)
  : fUseCollectionOnly(pset.get<bool>("UseCollectionOnly"))
  , fPFParticleLabel(pset.get<art::InputTag>("PFParticleLabel"))
  , fSCEXFlip(pset.get<bool>("SCEXFlip"))
  , fSCE(lar::providerFrom<spacecharge::SpaceChargeService>())
  , fInitialTrackInputLabel(pset.get<std::string>("InitialTrackInputLabel"))
  , fShowerStartPositionInputLabel(pset.get<std::string>("ShowerStartPositionInputLabel"))
  , fShowerDirectionInputLabel(pset.get<std::string>("ShowerDirectionInputLabel"))
  , fInitialTrackSpacePointsInputLabel(pset.get<std::string>("InitialTrackSpacePointsInputLabel"))
{}

//Order the shower hits with regards to their projected length onto
//the shower direction from the shower start position. This is done
//in the 2D coordinate system (wire direction, x)
void
shower::LArPandoraShowerAlg::OrderShowerHits(detinfo::DetectorPropertiesData const& detProp,
                                             std::vector<art::Ptr<recob::Hit>>& hits,
                                             TVector3 const& ShowerStartPosition,
                                             TVector3 const& ShowerDirection) const
{

  std::map<double, art::Ptr<recob::Hit>> OrderedHits;
  art::Ptr<recob::Hit> startHit = hits.front();

  //Get the wireID
  const geo::WireID startWireID = startHit->WireID();

  //Get the plane
  const geo::PlaneID planeid = startWireID.asPlaneID();

  //Get the pitch
  double pitch = fGeom->WirePitch(planeid);

  TVector2 Shower2DStartPosition = {
    fGeom->WireCoordinate(ShowerStartPosition, startHit->WireID().planeID()) * pitch,
    ShowerStartPosition.X()};

  //Vector of the plane
  TVector3 PlaneDirection = fGeom->Plane(planeid).GetIncreasingWireDirection();

  //get the shower 2D direction
  TVector2 Shower2DDirection = {ShowerDirection.Dot(PlaneDirection), ShowerDirection.X()};

  Shower2DDirection = Shower2DDirection.Unit();

  for (auto const& hit : hits) {

    //Get the wireID
    const geo::WireID WireID = hit->WireID();

    if (WireID.asPlaneID() != startWireID.asPlaneID()) { break; }

    //Get the hit Vector.
    TVector2 hitcoord = HitCoordinates(detProp, hit);

    //Order the hits based on the projection
    TVector2 pos = hitcoord - Shower2DStartPosition;
    OrderedHits[pos * Shower2DDirection] = hit;
  }

  //Transform the shower.
  std::vector<art::Ptr<recob::Hit>> showerHits;
  std::transform(OrderedHits.begin(),
                 OrderedHits.end(),
                 std::back_inserter(showerHits),
                 [](std::pair<double, art::Ptr<recob::Hit>> const& hit) { return hit.second; });

  //Sometimes get the order wrong. Depends on direction compared to the plane Correct for it here:
  art::Ptr<recob::Hit> frontHit = showerHits.front();
  art::Ptr<recob::Hit> backHit = showerHits.back();

  //Get the hit Vector.
  TVector2 fronthitcoord = HitCoordinates(detProp, frontHit);
  TVector2 frontpos = fronthitcoord - Shower2DStartPosition;

  //Get the hit Vector.
  TVector2 backhitcoord = HitCoordinates(detProp, backHit);
  TVector2 backpos = backhitcoord - Shower2DStartPosition;

  double frontproj = frontpos * Shower2DDirection;
  double backproj = backpos * Shower2DDirection;
  if (std::abs(backproj) < std::abs(frontproj)) {
    std::reverse(showerHits.begin(), showerHits.end());
  }

  hits = showerHits;
  return;
}

//Orders the shower spacepoints with regards to there perpendicular distance from
//the shower axis.
void
shower::LArPandoraShowerAlg::OrderShowerSpacePointsPerpendicular(
  std::vector<art::Ptr<recob::SpacePoint>>& showersps,
  TVector3 const& vertex,
  TVector3 const& direction) const
{

  std::map<double, art::Ptr<recob::SpacePoint>> OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.
  for (auto const& sp : showersps) {

    // Get the perpendicular distance
    double perp = SpacePointPerpendicular(sp, vertex, direction);

    //Add to the list
    OrderedSpacePoints[perp] = sp;
  }

  //Return an ordered list.
  showersps.clear();
  for (auto const& sp : OrderedSpacePoints) {
    showersps.push_back(sp.second);
  }
}

//Orders the shower spacepoints with regards to there prejected length from
//the shower start position in the shower direction.
void
shower::LArPandoraShowerAlg::OrderShowerSpacePoints(
  std::vector<art::Ptr<recob::SpacePoint>>& showersps,
  TVector3 const& vertex,
  TVector3 const& direction) const
{

  std::map<double, art::Ptr<recob::SpacePoint>> OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.
  for (auto const& sp : showersps) {

    // Get the projection of the space point along the direction
    double len = SpacePointProjection(sp, vertex, direction);

    //Add to the list
    OrderedSpacePoints[len] = sp;
  }

  //Return an ordered list.
  showersps.clear();
  for (auto const& sp : OrderedSpacePoints) {
    showersps.push_back(sp.second);
  }
}

void
shower::LArPandoraShowerAlg::OrderShowerSpacePoints(
  std::vector<art::Ptr<recob::SpacePoint>>& showersps,
  TVector3 const& vertex) const
{

  std::map<double, art::Ptr<recob::SpacePoint>> OrderedSpacePoints;

  //Loop over the spacepoints and get the pojected distance from the vertex.
  for (auto const& sp : showersps) {

    //Get the distance away from the start
    double mag = (SpacePointPosition(sp) - vertex).Mag();

    //Add to the list
    OrderedSpacePoints[mag] = sp;
  }

  //Return an ordered list.
  showersps.clear();
  for (auto const& sp : OrderedSpacePoints) {
    showersps.push_back(sp.second);
  }
}

TVector3
shower::LArPandoraShowerAlg::ShowerCentre(
  std::vector<art::Ptr<recob::SpacePoint>> const& showersps) const
{

  if (showersps.empty()) return TVector3{};

  TVector3 centre_position;
  for (auto const& sp : showersps) {
    TVector3 pos = SpacePointPosition(sp);
    centre_position += pos;
  }
  centre_position *= (1. / showersps.size());

  return centre_position;
}

TVector3
shower::LArPandoraShowerAlg::ShowerCentre(
  detinfo::DetectorClocksData const& clockData,
  detinfo::DetectorPropertiesData const& detProp,
  std::vector<art::Ptr<recob::SpacePoint>> const& showerspcs,
  art::FindManyP<recob::Hit> const& fmh) const
{

  float totalCharge = 0;
  return shower::LArPandoraShowerAlg::ShowerCentre(
    clockData, detProp, showerspcs, fmh, totalCharge);
}

//Returns the vector to the shower centre and the total charge of the shower.
TVector3
shower::LArPandoraShowerAlg::ShowerCentre(detinfo::DetectorClocksData const& clockData,
                                          detinfo::DetectorPropertiesData const& detProp,
                                          std::vector<art::Ptr<recob::SpacePoint>> const& showersps,
                                          art::FindManyP<recob::Hit> const& fmh,
                                          float& totalCharge) const
{

  TVector3 pos, chargePoint = TVector3(0, 0, 0);

  //Loop over the spacepoints and get the charge weighted center.
  for (auto const& sp : showersps) {

    //Get the position of the spacepoint
    pos = SpacePointPosition(sp);

    //Get the associated hits
    std::vector<art::Ptr<recob::Hit>> const& hits = fmh.at(sp.key());

    //Average the charge unless sepcified.
    float charge = 0;
    float charge2 = 0;
    for (auto const& hit : hits) {

      if (fUseCollectionOnly) {
        if (hit->SignalType() == geo::kCollection) {
          charge = hit->Integral();
          //Correct for the lifetime: Need to do other detproperites
          charge *= std::exp((sampling_rate(clockData) * hit->PeakTime()) /
                             (detProp.ElectronLifetime() * 1e3));
          break;
        }
      }
      else {

        //Correct for the lifetime FIX: Need  to do other detproperties somehow
        double Q = hit->Integral() * std::exp((sampling_rate(clockData) * hit->PeakTime()) /
                                              (detProp.ElectronLifetime() * 1e3));

        charge += Q;
        charge2 += Q * Q;
      }
    }

    if (!fUseCollectionOnly) {
      //Calculate the unbiased standard deviation and mean.
      float mean = charge / ((float)hits.size());

      float rms = 1;

      if (hits.size() > 1) {
        rms = std::sqrt((charge2 - charge * charge) / ((float)(hits.size() - 1)));
      }

      charge = 0;
      int n = 0;
      for (auto const& hit : hits) {
        double lifetimecorrection = std::exp((sampling_rate(clockData) * hit->PeakTime()) /
                                             (detProp.ElectronLifetime() * 1e3));
        if (hit->Integral() * lifetimecorrection > (mean - 2 * rms) &&
            hit->Integral() * lifetimecorrection < (mean + 2 * rms)) {
          charge += hit->Integral() * lifetimecorrection;
          ++n;
        }
      }

      if (n == 0) {
        mf::LogWarning("LArPandoraShowerAlg") << "no points used to make the charge value. \n";
      }

      charge /= n;
    }

    chargePoint += charge * pos;
    totalCharge += charge;

    if (charge == 0) {
      mf::LogWarning("LArPandoraShowerAlg") << "Averaged charge, within 2 sigma, for a spacepoint "
                                               "is zero, Maybe this not a good method. \n";
    }
  }

  double intotalcharge = 1 / totalCharge;
  TVector3 centre = chargePoint * intotalcharge;
  return centre;
}

//Return the spacepoint position in 3D cartesian coordinates.
TVector3
shower::LArPandoraShowerAlg::SpacePointPosition(art::Ptr<recob::SpacePoint> const& sp) const
{

  const Double32_t* sp_xyz = sp->XYZ();
  return TVector3{sp_xyz[0], sp_xyz[1], sp_xyz[2]};
}

double
shower::LArPandoraShowerAlg::DistanceBetweenSpacePoints(
  art::Ptr<recob::SpacePoint> const& sp_a,
  art::Ptr<recob::SpacePoint> const& sp_b) const
{
  TVector3 position_a = SpacePointPosition(sp_a);
  TVector3 position_b = SpacePointPosition(sp_b);
  return (position_a - position_b).Mag();
}

//Return the charge of the spacepoint in ADC.
double
shower::LArPandoraShowerAlg::SpacePointCharge(art::Ptr<recob::SpacePoint> const& sp,
                                              art::FindManyP<recob::Hit> const& fmh) const
{

  double Charge = 0;

  //Average over the charge even though there is only one
  std::vector<art::Ptr<recob::Hit>> const& hits = fmh.at(sp.key());
  for (auto const& hit : hits) {
    Charge += hit->Integral();
  }

  Charge /= (float)hits.size();

  return Charge;
}

//Return the spacepoint time.
double
shower::LArPandoraShowerAlg::SpacePointTime(art::Ptr<recob::SpacePoint> const& sp,
                                            art::FindManyP<recob::Hit> const& fmh) const
{

  double Time = 0;

  //Average over the hits
  std::vector<art::Ptr<recob::Hit>> const& hits = fmh.at(sp.key());
  for (auto const& hit : hits) {
    Time += hit->PeakTime();
  }

  Time /= (float)hits.size();
  return Time;
}

//Return the cooordinates of the hit in cm in wire direction and x.
TVector2
shower::LArPandoraShowerAlg::HitCoordinates(detinfo::DetectorPropertiesData const& detProp,
                                            art::Ptr<recob::Hit> const& hit) const
{

  //Get the pitch
  const geo::WireID WireID = hit->WireID();
  const geo::PlaneID planeid = WireID.asPlaneID();
  double pitch = fGeom->WirePitch(planeid);

  return TVector2(WireID.Wire * pitch, detProp.ConvertTicksToX(hit->PeakTime(), planeid));
}

double
shower::LArPandoraShowerAlg::SpacePointProjection(const art::Ptr<recob::SpacePoint>& sp,
                                                  TVector3 const& vertex,
                                                  TVector3 const& direction) const
{

  // Get the position of the spacepoint
  TVector3 pos = shower::LArPandoraShowerAlg::SpacePointPosition(sp) - vertex;

  // Get the the projected length
  return pos.Dot(direction);
}

double
shower::LArPandoraShowerAlg::SpacePointPerpendicular(art::Ptr<recob::SpacePoint> const& sp,
                                                     TVector3 const& vertex,
                                                     TVector3 const& direction) const
{

  // Get the projection of the spacepoint
  double proj = shower::LArPandoraShowerAlg::SpacePointProjection(sp, vertex, direction);

  return shower::LArPandoraShowerAlg::SpacePointPerpendicular(sp, vertex, direction, proj);
}

double
shower::LArPandoraShowerAlg::SpacePointPerpendicular(art::Ptr<recob::SpacePoint> const& sp,
                                                     TVector3 const& vertex,
                                                     TVector3 const& direction,
                                                     double proj) const
{

  // Get the position of the spacepoint
  TVector3 pos = shower::LArPandoraShowerAlg::SpacePointPosition(sp) - vertex;

  // Take away the projection * distance to find the perpendicular vector
  pos = pos - proj * direction;

  // Get the the projected length
  return pos.Mag();
}

//Function to calculate the RMS at segments of the shower and calculate the gradient of this. If negative then the direction is pointing the opposite way to the correct one
double
shower::LArPandoraShowerAlg::RMSShowerGradient(std::vector<art::Ptr<recob::SpacePoint>>& sps,
                                               const TVector3& ShowerCentre,
                                               const TVector3& Direction,
                                               const unsigned int nSegments) const
{

  if (nSegments == 0)
    throw cet::exception("LArPandoraShowerAlg")
      << "Unable to calculate RMS Shower Gradient with 0 segments" << std::endl;

  if (sps.size() < 3) return 0;

  //Order the spacepoints
  this->OrderShowerSpacePoints(sps, ShowerCentre, Direction);

  //Get the length of the shower.
  const double minProj = this->SpacePointProjection(sps[0], ShowerCentre, Direction);
  const double maxProj = this->SpacePointProjection(sps[sps.size() - 1], ShowerCentre, Direction);

  const double length = (maxProj - minProj);
  const double segmentsize = length / nSegments;

  if (segmentsize < std::numeric_limits<double>::epsilon()) return 0;

  std::map<int, std::vector<float>> len_segment_map;

  //Split the the spacepoints into segments.
  for (auto const& sp : sps) {

    //Get the the projected length
    const double len = this->SpacePointProjection(sp, ShowerCentre, Direction);

    //Get the length to the projection
    const double len_perp = this->SpacePointPerpendicular(sp, ShowerCentre, Direction, len);

    int sg_len = std::round(len / segmentsize);
    len_segment_map[sg_len].push_back(len_perp);
  }

  int counter = 0;
  float sumx = 0.f;
  float sumy = 0.f;
  float sumx2 = 0.f;
  float sumxy = 0.f;

  //Get the rms of the segments and caclulate the gradient.
  for (auto const& segment : len_segment_map) {
    if (segment.second.size() < 2) continue;
    float RMS = this->CalculateRMS(segment.second);

    // Check if the calculation failed
    if (RMS < 0) continue;

    //Calculate the gradient using regression
    sumx += segment.first;
    sumy += RMS;
    sumx2 += segment.first * segment.first;
    sumxy += RMS * segment.first;
    ++counter;
  }

  const float denom = counter * sumx2 - sumx * sumx;

  return std::abs(denom) < std::numeric_limits<float>::epsilon() ?
           0 :
           (counter * sumxy - sumx * sumy) / denom;
}

double
shower::LArPandoraShowerAlg::CalculateRMS(const std::vector<float>& perps) const
{

  if (perps.size() < 2) return std::numeric_limits<double>::lowest();

  float sum = 0;
  for (const auto& perp : perps) {
    sum += perp * perp;
  }

  return std::sqrt(sum / (perps.size() - 1));
}

double
shower::LArPandoraShowerAlg::SCECorrectPitch(double const& pitch,
                                             TVector3 const& pos,
                                             TVector3 const& dir,
                                             unsigned int const& TPC) const
{
  const geo::Point_t geoPos{pos.X(), pos.Y(), pos.z()};
  const geo::Vector_t geoDir{dir.X(), dir.Y(), dir.Z()};
  return shower::LArPandoraShowerAlg::SCECorrectPitch(pitch, geoPos, geoDir, TPC);
}
double
shower::LArPandoraShowerAlg::SCECorrectPitch(double const& pitch,
                                             geo::Point_t const& pos,
                                             geo::Vector_t const& dir,
                                             unsigned int const& TPC) const
{

  if (!fSCE || !fSCE->EnableCalSpatialSCE()) {
    throw cet::exception("LArPandoraShowerALG")
      << "Trying to correct SCE pitch when service is not configured" << std::endl;
  }
  // As the input pos is sce corrected already, find uncorrected pos
  const geo::Point_t uncorrectedPos = pos + fSCE->GetPosOffsets(pos);
  //Get the size of the correction at pos
  const geo::Vector_t posOffset = fSCE->GetCalPosOffsets(uncorrectedPos, TPC);

  //Get the position of next hit
  const geo::Point_t nextPos = uncorrectedPos + pitch * dir;
  //Get the offsets at the next pos
  const geo::Vector_t nextPosOffset = fSCE->GetCalPosOffsets(nextPos, TPC);

  //Calculate the corrected pitch
  const int xFlip(fSCEXFlip ? -1 : 1);
  geo::Vector_t pitchVec{pitch * dir.X() + xFlip * (nextPosOffset.X() - posOffset.X()),
                         pitch * dir.Y() + (nextPosOffset.Y() - posOffset.Y()),
                         pitch * dir.Z() + (nextPosOffset.Z() - posOffset.Z())};

  return pitchVec.r();
}

double
shower::LArPandoraShowerAlg::SCECorrectEField(double const& EField, TVector3 const& pos) const
{
  const geo::Point_t geoPos{pos.X(), pos.Y(), pos.z()};
  return shower::LArPandoraShowerAlg::SCECorrectEField(EField, geoPos);
}
double
shower::LArPandoraShowerAlg::SCECorrectEField(double const& EField, geo::Point_t const& pos) const
{

  // Check the space charge service is properly configured
  if (!fSCE || !fSCE->EnableSimEfieldSCE()) {
    throw cet::exception("LArPandoraShowerALG")
      << "Trying to correct SCE EField when service is not configured" << std::endl;
  }
  // Gets relative E field Distortions
  geo::Vector_t EFieldOffsets = fSCE->GetEfieldOffsets(pos);
  // Add 1 in X direction as this is the direction of the drift field
  EFieldOffsets += geo::Vector_t{1, 0, 0};
  // Convert to Absolute E Field from relative
  EFieldOffsets *= EField;
  // We only care about the magnitude for recombination
  return EFieldOffsets.r();
}

void
shower::LArPandoraShowerAlg::DebugEVD(art::Ptr<recob::PFParticle> const& pfparticle,
                                      art::Event const& Event,
                                      reco::shower::ShowerElementHolder const& ShowerEleHolder,
                                      std::string const& evd_disp_name_append) const
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
  if (evd_disp_name_append.length() > 0) canvasName += "_" + evd_disp_name_append;
  TCanvas* canvas = tfs->make<TCanvas>(canvasName, canvasName);

  // Initialise variables
  double x = 0;
  double y = 0;
  double z = 0;

  double x_min = std::numeric_limits<double>::max(), x_max = -std::numeric_limits<double>::max();
  double y_min = std::numeric_limits<double>::max(), y_max = -std::numeric_limits<double>::max();
  double z_min = std::numeric_limits<double>::max(), z_max = -std::numeric_limits<double>::max();

  // Get a bunch of associations (again)
  // N.B. this is a horribly inefficient way of doing things but as this is only
  // going to be used to debug I don't care, I would rather have generality in this case

  auto const pfpHandle = Event.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);

  // Get the spacepoint - PFParticle assn
  art::FindManyP<recob::SpacePoint> fmspp(pfpHandle, Event, fPFParticleLabel);
  if (!fmspp.isValid()) {
    throw cet::exception("LArPandoraShowerAlg") << "Trying to get the spacepoint and failed. Somet\
      hing is not configured correctly. Stopping ";
  }

  // Get the SpacePoints
  std::vector<art::Ptr<recob::SpacePoint>> const& spacePoints = fmspp.at(pfparticle.key());

  //We cannot progress with no spacepoints.
  if (spacePoints.empty()) { return; }

  // Get info from shower property holder
  TVector3 showerStartPosition = {-999, -999, -999};
  TVector3 showerDirection = {-999, -999, -999};
  std::vector<art::Ptr<recob::SpacePoint>> trackSpacePoints;

  //######################
  //### Start Position ###
  //######################
  double startXYZ[3] = {-999, -999, -999};
  if (!ShowerEleHolder.CheckElement(fShowerStartPositionInputLabel)) {
    mf::LogError("LArPandoraShowerAlg") << "Start position not set, returning " << std::endl;
    // return;
  }
  else {
    ShowerEleHolder.GetElement(fShowerStartPositionInputLabel, showerStartPosition);
    // Create 3D point at vertex, chosed to be origin for ease of use of display
    startXYZ[0] = showerStartPosition.X();
    startXYZ[1] = showerStartPosition.Y();
    startXYZ[2] = showerStartPosition.Z();
  }
  auto startPoly = std::make_unique<TPolyMarker3D>(1, startXYZ);

  //########################
  //### Shower Direction ###
  //########################

  double xDirPoints[2] = {-999, -999};
  double yDirPoints[2] = {-999, -999};
  double zDirPoints[2] = {-999, -999};

  //initialise counter point
  int point = 0;

  // Make 3D points for each spacepoint in the shower
  auto allPoly = std::make_unique<TPolyMarker3D>(spacePoints.size());

  if (!ShowerEleHolder.CheckElement(fShowerDirectionInputLabel) &&
      !ShowerEleHolder.CheckElement("ShowerStartPosition")) {
    mf::LogError("LArPandoraShowerAlg") << "Direction not set, returning " << std::endl;
  }
  else {

    // Get the min and max projections along the direction to know how long to draw
    ShowerEleHolder.GetElement(fShowerDirectionInputLabel, showerDirection);

    // the direction line
    double minProj = std::numeric_limits<double>::max();
    double maxProj = -std::numeric_limits<double>::max();

    //initialise counter point
    int point = 0;

    for (auto spacePoint : spacePoints) {
      TVector3 pos = shower::LArPandoraShowerAlg::SpacePointPosition(spacePoint);

      x = pos.X();
      y = pos.Y();
      z = pos.Z();
      allPoly->SetPoint(point, x, y, z);
      ++point;

      x_min = std::min(x, x_min);
      x_max = std::max(x, x_max);
      y_min = std::min(y, y_min);
      y_max = std::max(y, y_max);
      z_min = std::min(z, z_min);
      z_max = std::max(z, z_max);

      // Calculate the projection of (point-startpoint) along the direction
      double proj = shower::LArPandoraShowerAlg::SpacePointProjection(
        spacePoint, showerStartPosition, showerDirection);
      maxProj = std::max(proj, maxProj);
      minProj = std::min(proj, minProj);
    } // loop over spacepoints

    xDirPoints[0] = (showerStartPosition.X() + minProj * showerDirection.X());
    xDirPoints[1] = (showerStartPosition.X() + maxProj * showerDirection.X());

    yDirPoints[0] = (showerStartPosition.Y() + minProj * showerDirection.Y());
    yDirPoints[1] = (showerStartPosition.Y() + maxProj * showerDirection.Y());

    zDirPoints[0] = (showerStartPosition.Z() + minProj * showerDirection.Z());
    zDirPoints[1] = (showerStartPosition.Z() + maxProj * showerDirection.Z());
  }

  auto dirPoly = std::make_unique<TPolyLine3D>(2, xDirPoints, yDirPoints, zDirPoints);

  //#########################
  //### Initial Track SPs ###
  //#########################

  auto trackPoly = std::make_unique<TPolyMarker3D>(trackSpacePoints.size());
  if (!ShowerEleHolder.CheckElement(fInitialTrackSpacePointsInputLabel)) {
    mf::LogError("LArPandoraShowerAlg") << "TrackSpacePoints not set, returning " << std::endl;
    //    return;
  }
  else {
    ShowerEleHolder.GetElement(fInitialTrackSpacePointsInputLabel, trackSpacePoints);
    point = 0; // re-initialise counter
    for (auto spacePoint : trackSpacePoints) {
      TVector3 pos = shower::LArPandoraShowerAlg::SpacePointPosition(spacePoint);

      x = pos.X();
      y = pos.Y();
      z = pos.Z();
      trackPoly->SetPoint(point, x, y, z);
      ++point;
      x_min = std::min(x, x_min);
      x_max = std::max(x, x_max);
      y_min = std::min(y, y_min);
      y_max = std::max(y, y_max);
      z_min = std::min(z, z_min);
      z_max = std::max(z, z_max);
    } // loop over track spacepoints
  }

  //#########################
  //### Other PFParticles ###
  //#########################

  //  we want to draw all of the PFParticles in the event
  //Get the PFParticles
  std::vector<art::Ptr<recob::PFParticle>> pfps;
  art::fill_ptr_vector(pfps, pfpHandle);

  // initialse counters
  // Split into tracks and showers to make it clearer what pandora is doing
  int pfpTrackCounter = 0;
  int pfpShowerCounter = 0;

  // initial loop over pfps to find nuber of spacepoints for tracks and showers
  for (auto const& pfp : pfps) {
    std::vector<art::Ptr<recob::SpacePoint>> const& sps = fmspp.at(pfp.key());
    // If running pandora cheating it will call photons pdg 22
    int pdg = abs(pfp->PdgCode()); // Track or shower
    if (pdg == 11 || pdg == 22) { pfpShowerCounter += sps.size(); }
    else {
      pfpTrackCounter += sps.size();
    }
  }

  auto pfpPolyTrack = std::make_unique<TPolyMarker3D>(pfpTrackCounter);
  auto pfpPolyShower = std::make_unique<TPolyMarker3D>(pfpShowerCounter);

  // initialise counters
  int trackPoints = 0;
  int showerPoints = 0;

  for (auto const& pfp : pfps) {
    std::vector<art::Ptr<recob::SpacePoint>> const& sps = fmspp.at(pfp.key());
    int pdg = abs(pfp->PdgCode()); // Track or shower
    for (auto sp : sps) {
      //TVector3 pos = shower::LArPandoraShowerAlg::SpacePointPosition(sp) - showerStartPosition;
      TVector3 pos = shower::LArPandoraShowerAlg::SpacePointPosition(sp);

      x = pos.X();
      y = pos.Y();
      z = pos.Z();
      x_min = std::min(x, x_min);
      x_max = std::max(x, x_max);
      y_min = std::min(y, y_min);
      y_max = std::max(y, y_max);
      z_min = std::min(z, z_min);
      z_max = std::max(z, z_max);

      // If running pandora cheating it will call photons pdg 22
      if (pdg == 11 || pdg == 22) {
        pfpPolyShower->SetPoint(showerPoints, x, y, z);
        ++showerPoints;
      }
      else {
        pfpPolyTrack->SetPoint(trackPoints, x, y, z);
        ++trackPoints;
      }
    } // loop over sps

    //pfpPolyTrack->Draw();

  } // if (fDrawAllPFPs)

  //#################################
  //### Initial Track Traj Points ###
  //#################################

  auto TrackTrajPoly = std::make_unique<TPolyMarker3D>(TPolyMarker3D(1));
  auto TrackInitTrajPoly = std::make_unique<TPolyMarker3D>(TPolyMarker3D(1));

  if (ShowerEleHolder.CheckElement(fInitialTrackInputLabel)) {

    //Get the track
    recob::Track InitialTrack;
    ShowerEleHolder.GetElement(fInitialTrackInputLabel, InitialTrack);

    if (InitialTrack.NumberTrajectoryPoints() != 0) {

      point = 0;
      // Make 3D points for each trajectory point in the track stub
      for (unsigned int traj = 0; traj < InitialTrack.NumberTrajectoryPoints(); ++traj) {

        //ignore bogus info.
        auto flags = InitialTrack.FlagsAtPoint(traj);
        if (flags.isSet(recob::TrajectoryPointFlagTraits::NoPoint)) { continue; }

        geo::Point_t TrajPositionPoint = InitialTrack.LocationAtPoint(traj);
        TVector3 TrajPosition = {
          TrajPositionPoint.X(), TrajPositionPoint.Y(), TrajPositionPoint.Z()};

        TVector3 pos = TrajPosition;

        x = pos.X();
        y = pos.Y();
        z = pos.Z();
        TrackTrajPoly->SetPoint(point, x, y, z);
        ++point;
      } // loop over trajectory points

      geo::Point_t TrajInitPositionPoint = InitialTrack.LocationAtPoint(0);
      TVector3 TrajPosition = {
        TrajInitPositionPoint.X(), TrajInitPositionPoint.Y(), TrajInitPositionPoint.Z()};
      TVector3 pos = TrajPosition;
      x = pos.X();
      y = pos.Y();
      z = pos.Z();
      TrackInitTrajPoly->SetPoint(TrackInitTrajPoly->GetN(), x, y, z);
    }
  }

  gStyle->SetOptStat(0);
  TH3F axes("axes", "", 1, x_min, x_max, 1, y_min, y_max, 1, z_min, z_max);
  axes.SetDirectory(0);
  axes.GetXaxis()->SetTitle("X");
  axes.GetYaxis()->SetTitle("Y");
  axes.GetZaxis()->SetTitle("Z");
  axes.Draw();

  // Draw all of the things
  pfpPolyShower->SetMarkerStyle(20);
  pfpPolyShower->SetMarkerColor(4);
  pfpPolyShower->Draw();
  pfpPolyTrack->SetMarkerStyle(20);
  pfpPolyTrack->SetMarkerColor(6);
  pfpPolyTrack->Draw();
  allPoly->SetMarkerStyle(20);
  allPoly->Draw();
  trackPoly->SetMarkerStyle(20);
  trackPoly->SetMarkerColor(2);
  trackPoly->Draw();
  startPoly->SetMarkerStyle(21);
  startPoly->SetMarkerSize(0.5);
  startPoly->SetMarkerColor(3);
  startPoly->Draw();
  dirPoly->SetLineWidth(1);
  dirPoly->SetLineColor(6);
  dirPoly->Draw();
  TrackTrajPoly->SetMarkerStyle(22);
  TrackTrajPoly->SetMarkerColor(7);
  TrackTrajPoly->Draw();
  TrackInitTrajPoly->SetMarkerStyle(22);
  TrackInitTrajPoly->SetMarkerColor(4);
  TrackInitTrajPoly->Draw();

  // Save the canvas. Don't usually need this when using TFileService but this in the alg
  // not a module and didn't work without this so im going with it.
  canvas->Write();
}
