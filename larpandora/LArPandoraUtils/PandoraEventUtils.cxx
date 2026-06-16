/**
*
* @file larpandora/LArPandoraUtils/PandoraEventUtils.cxx
*
* @brief Utility containing helpful functions for end users to access products from events
*/

#include "larpandora/LArPandoraUtils/PandoraEventUtils.h"
#include "larpandora/LArPandoraUtils/PandoraPFParticleUtils.h"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace lar_pandora
{

std::vector<art::Ptr<recob::PFParticle>> PandoraEventUtils::GetPFParticles(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<recob::PFParticle>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Track>> PandoraEventUtils::GetTracks(const art::Event &evt, const std::string &label)
{
    mf::LogWarning("Pandora") << " Please note: accessing PFParticle tracks through this method is not the recommended workflow.\n"
                                 << " Please use PandoraEventUtils::GetPFParticles and access the tracks with PandoraPFParticleUtils::GetTrack."
                                 << std::endl;

    return PandoraEventUtils::GetProductVector<recob::Track>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Shower>> PandoraEventUtils::GetShowers(const art::Event &evt, const std::string &label)
{
    mf::LogWarning("Pandora") << " Please note: accessing PFParticle showers through this method is not the recommended workflow.\n"
                                 << " Please use PandoraEventUtils::GetPFParticles and access the tracks with PandoraPFParticleUtils::GetShower."
                                 << std::endl;

    return PandoraEventUtils::GetProductVector<recob::Shower>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Vertex>> PandoraEventUtils::GetVertices(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<recob::Vertex>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Hit>> PandoraEventUtils::GetHits(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<recob::Hit>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Wire>> PandoraEventUtils::GetWires(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<recob::Wire>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::SpacePoint>> PandoraEventUtils::GetSpacePoints(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<recob::SpacePoint>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::Slice>> PandoraEventUtils::GetSlices(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<recob::Slice>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<simb::MCTruth>> PandoraEventUtils::GetMCTruths(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<simb::MCTruth>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<simb::MCParticle>> PandoraEventUtils::GetMCParticles(const art::Event &evt, const std::string &label)
{
    return PandoraEventUtils::GetProductVector<simb::MCParticle>(evt,label);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

std::vector<art::Ptr<recob::PFParticle>> PandoraEventUtils::GetClearCosmics(const art::Event &evt, const std::string &label)
{
    std::vector<art::Ptr<recob::PFParticle>> theseParticles = PandoraEventUtils::GetProductVector<recob::PFParticle>(evt,label);

    std::vector<art::Ptr<recob::PFParticle>> theseCosmics;

    // We only want primary cosmic rays
    for (art::Ptr<recob::PFParticle> pParticle : theseParticles)
    {
        if (!pParticle->IsPrimary())
        {
            continue;
        } 

        if (PandoraPFParticleUtils::IsClearCosmic(pParticle, evt, label))
        {
            theseCosmics.push_back(pParticle);
        }
    }

    return theseCosmics;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

art::Ptr<recob::PFParticle> PandoraEventUtils::GetNeutrino(const art::Event &evt, const std::string &label)
{
    if (!HasNeutrino(evt,label))
    {
        throw cet::exception("Pandora") << "PandoraEventUtils::GetNeutrino --- No neutrino found";
    }
    
    art::Ptr<recob::PFParticle> pNeutrino;
    std::vector<art::Ptr<recob::PFParticle>> particles = PandoraEventUtils::GetPFParticles(evt,label);
    for (art::Ptr<recob::PFParticle> pParticle : particles)
    {
        if (PandoraPFParticleUtils::IsNeutrino(pParticle))
        {
            pNeutrino = pParticle;
            break; 
        }
    }
    return pNeutrino;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

bool PandoraEventUtils::HasNeutrino(const art::Event &evt, const std::string &label)
{
    bool hasNeutrino = false;
    std::vector<art::Ptr<recob::PFParticle>> particles = PandoraEventUtils::GetPFParticles(evt,label);
    for (art::Ptr<recob::PFParticle> pParticle : particles)
    {
        if(!pParticle->IsPrimary())
        {
            continue;
        }
        if (PandoraPFParticleUtils::IsNeutrino(pParticle))
        {
            hasNeutrino = true;
            break;
        }
    }
    return hasNeutrino;
}

} // namespace lar_pandora

