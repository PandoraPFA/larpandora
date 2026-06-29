////////////////////////////////////////////////////////////////////////
/// \file    ShowerVarManager.h
/// \brief   A class to manage the Ivysaurus shower variables input 
/// \author  Isobel Mawby - i.mawby1@lancaster.ac.uk
////////////////////////////////////////////////////////////////////////

#ifndef SHOWERMANAGER_H
#define SHOWERMANAGER_H

#include <vector>
#include <string>
#include <map>

#include "TVector3.h"

#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"

namespace ivysaurus
{

class ShowerVarManager
{
  public:
  class ShowerVars
  {
    public:
      ShowerVars();

      bool GetIsNormalised() const;
      std::pair<float, bool> GetDisplacement() const;
      std::pair<float, bool> GetDCA() const;
      std::pair<float, bool> GetTrackStubLength() const;
      std::pair<float, bool> GetNuVertexAvSeparation() const;
      std::pair<float, bool> GetNuVertexChargeAsymmetry() const;

      void SetIsNormalised(const bool isNormalised);
      void SetDisplacement(const float displacement);
      void SetDCA(const float dca);
      void SetTrackStubLength(const float trackStubLength);
      void SetNuVertexAvSeparation(const float nuVertexAvSeparation);
      void SetNuVertexChargeAsymmetry(const float nuVertexChargeAsymmetry);
      
    private:
      bool m_isNormalised;
      std::pair<float, bool> m_displacement;
      std::pair<float, bool> m_DCA;
      std::pair<float, bool> m_trackStubLength;
      std::pair<float, bool> m_nuVertexAvSeparation;
      std::pair<float, bool> m_nuVertexChargeAsymmetry;
  };

    ShowerVarManager(const fhicl::ParameterSet& pset);
    ~ShowerVarManager();

    bool EvaluateShowerVars(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, ShowerVarManager::ShowerVars &showerVars) const;
    void NormaliseShowerVars(ShowerVars &showerVars) const;

  private:
    bool GetParentEndpoint(const art::Event &evt, const art::Ptr<recob::PFParticle> &pfparticle, TVector3 &parentEnd) const;
    void FillDisplacement(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const;
    void FillTrackStub(const art::Event &evt, const art::Ptr<recob::Shower> shower, 
        ShowerVarManager::ShowerVars &showerVars) const;
    void FillAvSeparation(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::PFParticle> &pfparticle, 
        const art::Ptr<recob::Shower> &shower, ShowerVarManager::ShowerVars &showerVars) const;
    void FillChargeAsymmetry(const art::Event &evt, const TVector3 &parentEndpoint, const art::Ptr<recob::Shower> &shower, 
        ShowerVarManager::ShowerVars &showerVars) const;
    float GetViewChargeAsymmetry(const art::Event &evt, const TVector3 &parentEndpoint, const TVector3 &showerStart, 
        const std::vector<art::Ptr<recob::Hit>> &viewHits, const IvysaurusUtils::PandoraView &pandoraView) const;
    float NormaliseShowerVar(const std::pair<float, bool> &inputShowerVar, const float mean, const float std) const;    

    std::string m_recoModuleLabel;
    std::string m_trackModuleLabel;    
    std::string m_showerModuleLabel;
    std::string m_hitModuleLabel;
    float m_displacementMean;
    float m_displacementStd;
    float m_DCAMean;
    float m_DCAStd;
    float m_trackStubLengthMean;
    float m_trackStubLengthStd;
    float m_nuVertexAvSeparationMean;
    float m_nuVertexAvSeparationStd;
    float m_nuVertexChargeAsymmetryMean;
    float m_nuVertexChargeAsymmetryStd;
};

/////////////////////////////////////////////////////////////

inline bool ShowerVarManager::ShowerVars::GetIsNormalised() const
{
    return m_isNormalised;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetDisplacement() const
{
    return m_displacement;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetDCA() const
{
    return m_DCA;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetTrackStubLength() const
{
    return m_trackStubLength;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetNuVertexAvSeparation() const
{
    return m_nuVertexAvSeparation;
}

/////////////////////////////////////////////////////////////

inline std::pair<float, bool> ShowerVarManager::ShowerVars::GetNuVertexChargeAsymmetry() const
{
    return m_nuVertexChargeAsymmetry;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetDisplacement(const float displacement) 
{
    m_displacement.first = displacement;
    m_displacement.second = true;    
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetDCA(const float dca)
{
    m_DCA.first = dca;
    m_DCA.second = true;
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetTrackStubLength(const float trackStubLength)
{
    m_trackStubLength.first = trackStubLength;
    m_trackStubLength.second = true;    
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexAvSeparation(const float nuVertexAvSeparation)
{
    m_nuVertexAvSeparation.first = nuVertexAvSeparation;
    m_nuVertexAvSeparation.second = true;    
}

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetNuVertexChargeAsymmetry(const float nuVertexChargeAsymmetry)
{
    m_nuVertexChargeAsymmetry.first = nuVertexChargeAsymmetry;
    m_nuVertexChargeAsymmetry.second = true;    
}
    

/////////////////////////////////////////////////////////////

inline void ShowerVarManager::ShowerVars::SetIsNormalised(const bool isNormalised)
{
    m_isNormalised = isNormalised;
}

}

#endif  // SHOWERVARMANAGER_H
