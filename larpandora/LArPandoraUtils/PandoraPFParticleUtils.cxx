/**
*
* @file dunereco/AnaUtils/PandoraPFParticleUtils.cxx
*
* @brief Utility containing helpful functions for end users to access information about PFParticles
*/

#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

namespace lar_pandora
{

std::vector<art::Ptr<anab::T0>> PandoraPFParticleUtils::GetT0(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return PandoraPFParticleUtils::GetAssocProductVector<anab::T0>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<anab::CosmicTag>> PandoraPFParticleUtils::GetCosmicTag(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return PandoraPFParticleUtils::GetAssocProductVector<anab::CosmicTag>(pParticle,evt,label,label); 
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::PFParticle>> PandoraPFParticleUtils::GetChildParticles(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    auto theseParticles = evt.getHandle<std::vector<recob::PFParticle>>(label);
    bool success = theseParticles.isValid();
    
    if (!success)
    {   
        mf::LogError("LArPandora") << " Failed to find product with label " << label << " ... returning empty vector" << std::endl;
        return std::vector<art::Ptr<recob::PFParticle>>();
    }
    
    std::vector<art::Ptr<recob::PFParticle>> children;

    for (unsigned int iPart = 0; iPart < theseParticles->size(); ++iPart)
    {     
        if (theseParticles->at(iPart).Parent() == pParticle.key())
        {
            art::Ptr<recob::PFParticle> pChild(theseParticles,iPart);
            children.push_back(pChild);
        }
    }

    return children;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> PandoraPFParticleUtils::GetHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{    
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<art::Ptr<recob::Cluster>> theseClusters = PandoraPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<art::Ptr<recob::Hit>> theseHits;
    for (const art::Ptr<recob::Cluster> pCluster : theseClusters)
    {
      const std::vector<art::Ptr<recob::Hit>> tempHits = PandoraPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
      theseHits.insert(theseHits.end(),tempHits.begin(),tempHits.end());
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> PandoraPFParticleUtils::GetViewHits(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label, const unsigned short &view)
{
    // There isn't a direct association between PFParticles and hits, so we go via clusters
    const std::vector<art::Ptr<recob::Cluster>> theseClusters = PandoraPFParticleUtils::GetAssocProductVector<recob::Cluster>(pParticle,evt,label,label);

    std::vector<art::Ptr<recob::Hit>> theseHits;
    for (const art::Ptr<recob::Cluster> pCluster : theseClusters)
    { 
        const std::vector<art::Ptr<recob::Hit>> tempHits = PandoraPFParticleUtils::GetAssocProductVector<recob::Hit>(pCluster,evt,label,label);
        for(const art::Ptr<recob::Hit> pHit : tempHits)
        {
            if (pHit->View() == view)
            {
                theseHits.push_back(pHit);
            }
        }
    }
    return theseHits;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> PandoraPFParticleUtils::GetSpacePoints(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return PandoraPFParticleUtils::GetAssocProductVector<recob::SpacePoint>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Track> PandoraPFParticleUtils::GetTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    return PandoraPFParticleUtils::GetAssocProduct<recob::Track>(pParticle,evt,particleLabel,trackLabel);
}    

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Shower> PandoraPFParticleUtils::GetShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    return PandoraPFParticleUtils::GetAssocProduct<recob::Shower>(pParticle,evt,particleLabel,showerLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Vertex> PandoraPFParticleUtils::GetVertex(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    return PandoraPFParticleUtils::GetAssocProduct<recob::Vertex>(pParticle,evt,particleLabel,particleLabel);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::Slice> PandoraPFParticleUtils::GetSlice(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> pMetadata = PandoraPFParticleUtils::GetMetadata(pParticle,evt,label);
    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    unsigned int sliceIndex;
    std::map<std::string,float>::iterator mapItr = metaMap.find("SliceIndex");

    if (mapItr == metaMap.end())   
    {
        throw cet::exception("Pandora") << " PandoraPFParticleUtils::GetSlice --- No associated slice found";
    }
    else
    {
        sliceIndex = mapItr->second;
    }

    return PandoraPFParticleUtils::GetProductVector<recob::Slice>(evt,label).at(sliceIndex);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<larpandoraobj::PFParticleMetadata> PandoraPFParticleUtils::GetMetadata(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &label)
{
    return PandoraPFParticleUtils::GetAssocProduct<larpandoraobj::PFParticleMetadata>(pParticle,evt,label,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFParticleUtils::HasTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    // This function needs to fail if GetTrack would fail
    const std::vector<art::Ptr<recob::Track>> theseTracks = PandoraPFParticleUtils::GetAssocProductVector<recob::Track>(pParticle,evt,particleLabel,trackLabel);

    return !theseTracks.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFParticleUtils::HasShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    // This function needs to fail if GetShower would fail
    const std::vector<art::Ptr<recob::Shower>> theseShowers = PandoraPFParticleUtils::GetAssocProductVector<recob::Shower>(pParticle,evt,particleLabel,showerLabel);

    return !theseShowers.empty();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFParticleUtils::IsTrack(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &trackLabel)
{
    if (!PandoraPFParticleUtils::HasTrack(pParticle, evt, particleLabel, trackLabel))
        return false;

    return lar_pandora::LArPandoraHelper::IsTrack(pParticle);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFParticleUtils::IsShower(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel, const std::string &showerLabel)
{
    if (!PandoraPFParticleUtils::HasShower(pParticle, evt, particleLabel, showerLabel))
        return false;

    return lar_pandora::LArPandoraHelper::IsShower(pParticle);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFParticleUtils::IsClearCosmic(const art::Ptr<recob::PFParticle> &pParticle, const art::Event &evt, const std::string &particleLabel)
{
    const art::Ptr<larpandoraobj::PFParticleMetadata> pMetadata = PandoraPFParticleUtils::GetMetadata(pParticle,evt,particleLabel);

    std::map<std::string,float> metaMap = pMetadata->GetPropertiesMap();

    return metaMap.find("IsClearCosmic") != metaMap.end();
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraPFParticleUtils::IsNeutrino(const art::Ptr<recob::PFParticle> &pParticle)
{
    const int pdg = pParticle->PdgCode();
    return ((std::abs(pdg) == 12) || (std::abs(pdg) == 14) || (std::abs(pdg) ==16));
}

} // namespace lar_pandora


