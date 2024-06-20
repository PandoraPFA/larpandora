#ifndef LAR_PANDORA_GET_DETECTOR_TYPE_H
#define LAR_PANDORA_GET_DETECTOR_TYPE_H 1

#include <memory>

namespace lar_pandora {

  class LArPandoraDetectorType;

  namespace detector_functions {

    /**
         *  @brief  Factory class that returns the correct detector type interface
         *
         *  @result The detector type interface
         */
      std::unique_ptr<LArPandoraDetectorType> GetDetectorType();

  } // namespace detector_functions

} // namespace lar_pandora

#endif
