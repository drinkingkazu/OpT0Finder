## What is Flash Matching?
In our detector we assume particle interactions produce signals in both TPC and PMT. *Typically* an interaction signal reconstructed from PMT signal is represented by `recob::OpFlash` and a corresponding representation from the TPC side is a cluster of particles that are produced from the same interaction (or the same root particle, such as a neutrino = neutrino interaction, a cosmic ray, etc.). 

`recob::OpFlash` carries two information: a timing in the detector clock + photo-electron (PE) spectrum over PMTs, which is an estimate of photo-electron counts detected by each PMT. Although this is not a strict definition, we consider `recob::OpFlash` to represent one interaction without a granularity of particles in the same interaction. This is a good assumption as a photon detection timing resolution is several to 10 nano-seconds (prompt light time constant is ~6ns, and time of flight is several ns while particles finish traveling in similar time scale). In order to match with a `recob::OpFlash`, TPC interaction unit should be a set of particles that consists one interaction.

That's a long introduction... but finally! A flash matching is a reconstruction process in which reconstructed interactions (which unit is discussed above) from TPC and PMT are matched :)

### What is OpT0Finder?
OpT0Finder is a mini-framework for developing and running flash matching algorithms. The top-level program execution is handled by `flashmatch::FlashMatchManager`. If you want to just try running the code and learn by examples, you can skip this section and maybe look at [this] toy-simulation example. In below, I ignore `flashmatch` namespace when specifying C++ class under that scope.

In OpT0Finder, TPC and PMT interactions are defined as `QCluster_t` and `Flash_t` respectively. So `FlashMatchmanager` is a code to find the right match between sets of `QCluster_t` and `Flash_t`. Though this is not a strict requirement, in general the matching is done by producing a _hypothesis_ `Flash_t`, which is an expected PE per PMT based on `QCluster_t`, and comparing the similarity between the hypothesized and reconstructed (input, provided) `Flash_t`. In doing this, `FlashMatchManager` execute 5 types of algorithms.
* TPC interaction filter ... optional, can exclude some `QCluster_t` at run time before running matching
* PMT interaction filter ... optional, can exclude some `Flash_t` at run time before running matching
* Match prohibit ... optional, can exclude a certain `QCluster_t` and `Flash_t` from matching
* Hypothesis algorithm ... required, used to generate a hypothesis `Flash_t` from `QCluster_t`.
* Match algorithm ... required, used to match hypothesis `Flash_t` and reconstructed (input) `Flash_t`.

### Big picture: how does `FlashMatchManager` work?
`FlashMatchManager` works in 4 big steps:
0. Configure parameters
1. Register a list of `QCluster_t` = TPC interactions
2. Register a list of `Flash_t` = PMT interactions
3. Run flash matching
That's it! So, in principle, only 4 function calls.

You can look at `dat/flashmatch.cfg` as an example. By now you should be able to understand the `FlashMatchmanager` configuration. If you are interested in learning about algorithms, you can find README and the source code under `flashmatch/Algorithms` directory. For `FlashMatchManager` itself, the code is under `flashmatch/Base` directory.


