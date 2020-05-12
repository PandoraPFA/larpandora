/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraShowerCreation_module.cc
 *
 *  @brief  module for lar pandora shower creation
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

#include "larreco/Calorimetry/LinearEnergyAlg.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"

#include <memory>

namespace lar_pandora
{

class LArPandoraShowerCreation : public art::EDProducer
{
public:
    explicit LArPandoraShowerCreation(fhicl::ParameterSet const &pset);

    LArPandoraShowerCreation(LArPandoraShowerCreation const &) = delete;
    LArPandoraShowerCreation(LArPandoraShowerCreation &&) = delete;
    LArPandoraShowerCreation & operator = (LArPandoraShowerCreation const &) = delete;
    LArPandoraShowerCreation & operator = (LArPandoraShowerCreation &&) = delete;

    void produce(art::Event &evt) override;

private:
    /**
     *  @brief  Build a recob::Shower object
     *
     *  @param  id the id code for the shower
     *  @param  larShowerPCA the lar shower pca parameters extracted from pandora
     *  @param  vertexPosition the shower vertex position
     */
    recob::Shower BuildShower(const int id, const lar_content::LArShowerPCA &larShowerPCA, const pandora::CartesianVector &vertexPosition) const;

    /**
     *  @brief  Build a recob::PCAxis object
     *
     *  @param  larShowerPCA the lar shower pca parameters extracted from pandora
     */
    recob::PCAxis BuildPCAxis(const lar_content::LArShowerPCA &larShowerPCA) const;

    std::string     m_pfParticleLabel;              ///< The pf particle label
    bool            m_useAllParticles;              ///< Build a recob::Track for every recob::PFParticle

    // TODO When implementation lived in LArPandoraOutput, it contained key building blocks for calculation of shower energies per plane.
    // Now functionality has moved to separate module, will require reimplementation (was deeply embedded in LArPandoraOutput structure).
    // const calo::LinearEnergyAlg    *m_pShowerEnergyAlg;       ///< The address of the shower energy algorithm
};

DEFINE_ART_MODULE(LArPandoraShowerCreation)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

#include "art/Persistency/Common/PtrMaker.h"

#include "canvas/Utilities/InputTag.h"

#include "larcore/Geometry/Geometry.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PCAxis.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Vertex.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <iostream>

namespace lar_pandora
{

LArPandoraShowerCreation::LArPandoraShowerCreation(fhicl::ParameterSet const &pset) :
    EDProducer{pset},
    m_pfParticleLabel(pset.get<std::string>("PFParticleLabel")),
    m_useAllParticles(pset.get<bool>("UseAllParticles", false))
{
    produces< std::vector<recob::Shower> >();
    produces< std::vector<recob::PCAxis> >();
    produces< art::Assns<recob::PFParticle, recob::Shower> >();
    produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
    produces< art::Assns<recob::Shower, recob::Hit> >();
    produces< art::Assns<recob::Shower, recob::PCAxis> >();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPandoraShowerCreation::produce(art::Event &evt)
{
    std::unique_ptr< std::vector<recob::Shower> > outputShowers( new std::vector<recob::Shower> );
    std::unique_ptr< std::vector<recob::PCAxis> > outputPCAxes( new std::vector<recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower> > outputParticlesToShowers( new art::Assns<recob::PFParticle, recob::Shower> );
    std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis> > outputParticlesToPCAxes( new art::Assns<recob::PFParticle, recob::PCAxis> );
    std::unique_ptr< art::Assns<recob::Shower, recob::Hit> > outputShowersToHits( new art::Assns<recob::Shower, recob::Hit> );
    std::unique_ptr< art::Assns<recob::Shower, recob::PCAxis> > outputShowersToPCAxes( new art::Assns<recob::Shower, recob::PCAxis> );

    const art::PtrMaker<recob::Shower> makeShowerPtr(evt);
    const art::PtrMaker<recob::PCAxis> makePCAxisPtr(evt);

    int showerCounter(0);

    // Organise inputs
    PFParticleVector pfParticleVector, extraPfParticleVector;
    PFParticlesToSpacePoints pfParticlesToSpacePoints;
    PFParticlesToClusters pfParticlesToClusters;
    LArPandoraHelper::CollectPFParticles(evt, m_pfParticleLabel, pfParticleVector, pfParticlesToSpacePoints);
    LArPandoraHelper::CollectPFParticles(evt, m_pfParticleLabel, extraPfParticleVector, pfParticlesToClusters);

    VertexVector vertexVector;
    PFParticlesToVertices pfParticlesToVertices;
    LArPandoraHelper::CollectVertices(evt, m_pfParticleLabel, vertexVector, pfParticlesToVertices);

    for (const art::Ptr<recob::PFParticle> pPFParticle : pfParticleVector)
    {
        // Select shower-like pfparticles
        if (!m_useAllParticles && !LArPandoraHelper::IsShower(pPFParticle))
            continue;

        // Obtain associated spacepoints
        PFParticlesToSpacePoints::const_iterator particleToSpacePointIter(pfParticlesToSpacePoints.find(pPFParticle));

        if (pfParticlesToSpacePoints.end() == particleToSpacePointIter)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "No spacepoints associated to particle ";
            continue;
        }

        // Obtain associated clusters
        PFParticlesToClusters::const_iterator particleToClustersIter(pfParticlesToClusters.find(pPFParticle));

        if (pfParticlesToClusters.end() == particleToClustersIter)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "No clusters associated to particle ";
            continue;
        }

        // Obtain associated vertex
        PFParticlesToVertices::const_iterator particleToVertexIter(pfParticlesToVertices.find(pPFParticle));

        if ((pfParticlesToVertices.end() == particleToVertexIter) || (1 != particleToVertexIter->second.size()))
        {
            mf::LogDebug("LArPandoraShowerCreation") << "Unexpected number of vertices for particle ";
            continue;
        }

        // Copy information into expected pandora form
        pandora::CartesianPointVector cartesianPointVector;
        for (const art::Ptr<recob::SpacePoint> spacePoint : particleToSpacePointIter->second)
            cartesianPointVector.emplace_back(pandora::CartesianVector(spacePoint->XYZ()[0], spacePoint->XYZ()[1], spacePoint->XYZ()[2]));

        double vertexXYZ[3] = {0., 0., 0.};
        particleToVertexIter->second.front()->XYZ(vertexXYZ);
        const pandora::CartesianVector vertexPosition(vertexXYZ[0], vertexXYZ[1], vertexXYZ[2]);

        // Call pandora "fast" shower fitter
        try
        {
            // Access centroid of shower via this method
            const lar_content::LArShowerPCA initialLArShowerPCA(lar_content::LArPfoHelper::GetPrincipalComponents(cartesianPointVector, vertexPosition));

            // Ensure successful creation of all structures before placing results in output containers, remaking LArShowerPCA with updated vertex
            const pandora::CartesianVector &centroid(initialLArShowerPCA.GetCentroid());
            const pandora::CartesianVector &primaryAxis(initialLArShowerPCA.GetPrimaryAxis());
            const pandora::CartesianVector &secondaryAxis(initialLArShowerPCA.GetSecondaryAxis());
            const pandora::CartesianVector &tertiaryAxis(initialLArShowerPCA.GetTertiaryAxis());
            const pandora::CartesianVector &eigenvalues(initialLArShowerPCA.GetEigenValues());

            // Project the PFParticle vertex onto the PCA axis
            const pandora::CartesianVector projectedVertexPosition(centroid - primaryAxis.GetUnitVector() * (centroid - vertexPosition).GetDotProduct(primaryAxis));

            // By convention, principal axis should always point away from vertex
            const float testProjection(primaryAxis.GetDotProduct(projectedVertexPosition - centroid));
            const float directionScaleFactor((testProjection > std::numeric_limits<float>::epsilon()) ? -1.f : 1.f);

            const lar_content::LArShowerPCA larShowerPCA(centroid, primaryAxis * directionScaleFactor, secondaryAxis * directionScaleFactor, tertiaryAxis * directionScaleFactor, eigenvalues);
            const recob::Shower shower(LArPandoraShowerCreation::BuildShower(showerCounter++, larShowerPCA, projectedVertexPosition));
            const recob::PCAxis pcAxis(LArPandoraShowerCreation::BuildPCAxis(larShowerPCA));
            outputShowers->emplace_back(shower);
            outputPCAxes->emplace_back(pcAxis);
        }
        catch (const pandora::StatusCodeException &)
        {
            mf::LogDebug("LArPandoraShowerCreation") << "Unable to extract shower pca";
            continue;
        }

        // Output objects
        art::Ptr<recob::Shower> pShower(makeShowerPtr(outputShowers->size() - 1));
        art::Ptr<recob::PCAxis> pPCAxis(makePCAxisPtr(outputPCAxes->size() - 1));

        HitVector hitsInParticle;
        LArPandoraHelper::GetAssociatedHits(evt, m_pfParticleLabel, particleToClustersIter->second, hitsInParticle);

        // Output associations, after output objects are in place
        util::CreateAssn(*this, evt, pShower, pPFParticle, *(outputParticlesToShowers.get()));
        util::CreateAssn(*this, evt, pPCAxis, pPFParticle, *(outputParticlesToPCAxes.get()));
        util::CreateAssn(*this, evt, *(outputShowers.get()), hitsInParticle, *(outputShowersToHits.get()));
        util::CreateAssn(*this, evt, pPCAxis, pShower, *(outputShowersToPCAxes.get()));
    }

    mf::LogDebug("LArPandora") << "   Number of new showers: " << outputShowers->size() << std::endl;
    mf::LogDebug("LArPandora") << "   Number of new pcaxes:  " << outputPCAxes->size() << std::endl;

    evt.put(std::move(outputShowers));
    evt.put(std::move(outputPCAxes));
    evt.put(std::move(outputParticlesToShowers));
    evt.put(std::move(outputParticlesToPCAxes));
    evt.put(std::move(outputShowersToHits));
    evt.put(std::move(outputShowersToPCAxes));
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::Shower LArPandoraShowerCreation::BuildShower(const int id, const lar_content::LArShowerPCA &larShowerPCA, const pandora::CartesianVector &vertexPosition) const
{
    const pandora::CartesianVector &showerLength(larShowerPCA.GetAxisLengths());
    const pandora::CartesianVector &showerDirection(larShowerPCA.GetPrimaryAxis());

    const float length(showerLength.GetX());
    const float openingAngle(larShowerPCA.GetPrimaryLength() > 0.f ? std::atan(larShowerPCA.GetSecondaryLength() / larShowerPCA.GetPrimaryLength()) : 0.f);
    const TVector3 direction(showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ());
    const TVector3 vertex(vertexPosition.GetX(), vertexPosition.GetY(), vertexPosition.GetZ());

    // TODO
    const TVector3 directionErr;
    const TVector3 vertexErr;
    const std::vector<double> totalEnergyErr;
    const std::vector<double> dEdx;
    const std::vector<double> dEdxErr;
    const std::vector<double> totalEnergy;
    const int bestplane(0);

    return recob::Shower(direction, directionErr, vertex, vertexErr, totalEnergy, totalEnergyErr, dEdx, dEdxErr, bestplane, id, length, openingAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

recob::PCAxis LArPandoraShowerCreation::BuildPCAxis(const lar_content::LArShowerPCA &larShowerPCA) const
{
    const pandora::CartesianVector &showerCentroid(larShowerPCA.GetCentroid());
    const pandora::CartesianVector &showerDirection(larShowerPCA.GetPrimaryAxis());
    const pandora::CartesianVector &showerSecondaryVector(larShowerPCA.GetSecondaryAxis());
    const pandora::CartesianVector &showerTertiaryVector(larShowerPCA.GetTertiaryAxis());
    const pandora::CartesianVector &showerEigenValues(larShowerPCA.GetEigenValues());

    const bool svdOK(true); ///< SVD Decomposition was successful
    const double eigenValues[3] = {showerEigenValues.GetX(), showerEigenValues.GetY(), showerEigenValues.GetZ()}; ///< Eigen values from SVD decomposition
    const double avePosition[3] = {showerCentroid.GetX(), showerCentroid.GetY(), showerCentroid.GetZ()}; ///< Average position of hits fed to PCA

    std::vector< std::vector<double> > eigenVecs = { /// The three principle axes
        { showerDirection.GetX(), showerDirection.GetY(), showerDirection.GetZ() },
        { showerSecondaryVector.GetX(), showerSecondaryVector.GetY(), showerSecondaryVector.GetZ() },
        { showerTertiaryVector.GetX(), showerTertiaryVector.GetY(), showerTertiaryVector.GetZ() }
    };

    // TODO
    const int numHitsUsed(100); ///< Number of hits in the decomposition, not yet ready
    const double aveHitDoca(0.); ///< Average doca of hits used in PCA, not ready yet
    const size_t iD(util::kBogusI); ///< Axis ID, not ready yet

    return recob::PCAxis(svdOK, numHitsUsed, eigenValues, eigenVecs, avePosition, aveHitDoca, iD);
}

} // namespace lar_pandora
