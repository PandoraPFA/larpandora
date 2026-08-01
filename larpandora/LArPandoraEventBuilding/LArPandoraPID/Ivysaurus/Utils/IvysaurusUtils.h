/**
 *  @file   larpandora/LArPandoraEventBuilding/LArPandoraPID/Ivysaurus/Utils/IvysaurusUtils.h
 *
 *  @brief  Helpers used by the Ivysaurus PID
 */
#ifndef IVYSAURUS_UTILS_H
#define IVYSAURUS_UTILS_H
// ART
#include "art/Framework/Principal/Event.h"
// LArSoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
// ROOT
#include "TVector3.h"
//C++
#include <vector>

namespace IvysaurusUtils
{
    enum PandoraView {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W};

    /**
     *  @brief  Project a 3D coordinate, or LArSoft hit coordinate into an input Pandora 2D view
     *
     *  @param inputPosition the input coordinate
     *  @param pandoraView the input Pandora 2D view
     *
     *  @return the projected 2D position
     */        
    const TVector3 ProjectIntoPandoraView(const TVector3 &inputPosition3D, const PandoraView pandoraView);

    /**
     *  @brief  Obtain the Pandora view of a LArSoft hit
     *
     *  @param hit the input LArSoft hit
     *
     *  @return the Pandora view
     */ 
    const PandoraView GetPandora2DView(const art::Ptr<recob::Hit> &hit);

    /**
     *  @brief  Obtain the position and width of a LArSoft hit in Pandora coordinates
     *
     *  @param evt the art event
     *  @param hit the input LArSoft hit
     *  @param hitType the input Pandora coordinates
     *  @param[out] pandoraHitPosition the output position in pandora coordinates
     *  @param[out] width the output width in pandora coordinates
     */     
    void ObtainPandoraHitPositionAndWidth(const art::Event &evt, const art::Ptr<recob::Hit> hit, 
        const PandoraView hitType, TVector3 &pandoraHitPosition, float &width);

    /**
     *  @brief  Obtain Pandora U-coordinate from input y+z coordinates
     *
     *  @param yCoord the input y-coordinate
     *  @param zCoord the input z-coordinate
     *
     *  @return the output U-coordinate
     */
    float YZToU(const float yCoord, const float zCoord);

    /**
     *  @brief  Obtain Pandora V-coordinate from input y+z coordinates
     *
     *  @param yCoord the input y-coordinate
     *  @param zCoord the input z-coordinate
     *
     *  @return the output V-coordinate
     */    
    float YZToV(const float yCoord, const float zCoord);

    /**
     *  @brief  Obtain Pandora W-coordinate from input y+z coordinates
     *
     *  @param yCoord the input y-coordinate
     *  @param zCoord the input z-coordinate
     *
     *  @return the output W-coordinate
     */
    float YZToW(const float yCoord, const float zCoord);

    /**
     *  @brief  Determine the completeness of a MCParticle specified by an input TrackID
     *          from a group of input hits 
     *
     *  @param clockData the clock data service
     *  @param selectedHits the input group of hits
     *  @param eventHits the hits of the event
     *  @param trackID the trackID of the target MCParticle
     *
     *  @return the completeness
     */    
    double CompletenessFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits, 
        const std::vector<art::Ptr<recob::Hit>> &eventHits, const int trackID);

    /**
     *  @brief  Determine the purity of a MCParticle specified by an input TrackID
     *          from a group of input hits 
     *
     *  @param clockData the clock data service
     *  @param selectedHits the input group of hits
     *  @param trackID the trackID of the target MCParticle
     *
     *  @return the purity
     */        
    double HitPurityFromTrueParticleID(detinfo::DetectorClocksData const& clockData, const std::vector<art::Ptr<recob::Hit>> &selectedHits,
        const int trackID);

    /**
     *  @brief  Apply the trapezium rule to estimate the integral of Gaussian between two limits with a given step size
     *
     *  @param limitA the lower limit
     *  @param limitB the upper limit
     *  @param mean the mean of the gaussian
     *  @param std the standard deviation of the gaussian 
     *  @param stepSize the step size used by the integral estimator
     *
     *  @return the estimated integral
     */    
    float IntegrateGaussian(const float limitA, const float limitB, const float mean, const float std, const float stepSize);

    /**
     *  @brief  Implementation of the trapezium rule to calculate a segment of the integral of a defined Gaussian
     *
     *  @param limitA the lower limit
     *  @param limitB the upper limit
     *  @param mean the mean of the gaussian
     *  @param std the standard deviation of the gaussian 
     *
     *  @return the estimated integral
     */    
    float TrapeziumRule(const float lowerLimit, const float upperLimit, const float mean, const float std);
}

#endif // IVYSAURUS_UTILS_H

