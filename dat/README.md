## What should be in this directory?
Under `dat` directory, there should be some files:
* `detector_specs.cfg` ... used when outside LArSoft to define detector-specific parameters
* `flashmatch.cfg` ... serves as an example configuration file for toy MC sample and for instruction purpose.
* Photon library data file ... used when outside LArSoft but this file is not in git repository as it is huge. If you want to know where, try instantiating `phot::PhotonVisbilityService` through `flashmatch::DetectorSpecs`. It spits out an error message with the full download path.


