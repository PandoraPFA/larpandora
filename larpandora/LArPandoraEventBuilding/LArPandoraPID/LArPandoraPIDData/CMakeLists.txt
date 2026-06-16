cet_make_library(
  SOURCE
    PandoraPIDResult.cxx

  LIBRARIES
    PUBLIC
      ROOT::GenVector
      lardataobj::RecoBase
      canvas::canvas
)

install_headers()
install_source()