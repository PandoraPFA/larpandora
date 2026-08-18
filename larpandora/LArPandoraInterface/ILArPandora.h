/**
 *  @file   larpandora/LArPandoraInterface/ILArPandora.h
 *
 *  @brief  Interface class for LArPandora producer modules, which reconstruct recob::PFParticles from recob::Hits
 */

#ifndef I_LAR_PANDORA_H
#define I_LAR_PANDORA_H 1

#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace recob {
  class Hit;
  class OpHit;
}
namespace pandora {
  class Pandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_pandora {

  typedef std::map<int, art::Ptr<recob::Hit>> IdToHitMap;
  typedef std::map<int, art::Ptr<recob::OpHit>> IdToOpHitMap;

  /**
 *  @brief  ILArPandora class
 */
  class ILArPandora : public art::EDProducer {
  public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset the parameter set
     */
    ILArPandora(fhicl::ParameterSet const& pset);

    /**
     *  @brief  Destructor
     */
    virtual ~ILArPandora();

  protected:
    /**
     *  @brief  Create pandora instances
     */
    virtual void CreatePandoraInstances() = 0;

    /**
     *  @brief  Configure pandora instances
     */
    virtual void ConfigurePandoraInstances() = 0;

    /**
     *  @brief  Delete pandora instances
     */
    virtual void DeletePandoraInstances() = 0;

    /**
     *  @brief  Create pandora input hits, mc particles etc.
     *
     *  @param  evt the art event
     *  @param  idToHitMap to receive the populated pandora hit id to art hit map
     */
    virtual void CreatePandoraInput(art::Event& evt, IdToHitMap& idToHitMap) = 0;

    /**
     *  @brief  Process pandora output particle flow objects
     *
     *  @param  evt the art event
     *  @param  idToHitMap the pandora hit id to art hit map
     */
    virtual void ProcessPandoraOutput(art::Event& evt, const IdToHitMap& idToHitMap) = 0;

    /**
     *  @brief  Run all associated pandora instances
     */
    virtual void RunPandoraInstances() = 0;

    /**
     *  @brief  Reset all associated pandora instances
     */
    virtual void ResetPandoraInstances() = 0;

    const pandora::Pandora* m_pPrimaryPandora; ///< The address of the primary pandora instance
  };

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline ILArPandora::ILArPandora(fhicl::ParameterSet const& pset)
    : EDProducer(pset), m_pPrimaryPandora(nullptr)
  {}

  //------------------------------------------------------------------------------------------------------------------------------------------

  inline ILArPandora::~ILArPandora() {}

} // namespace lar_pandora

#endif // #ifndef I_LAR_PANDORA_H
