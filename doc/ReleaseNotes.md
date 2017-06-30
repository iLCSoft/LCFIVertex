# v00-07-04

* 2017-06-20 Andre Sailer ([PR#4](https://github.com/iLCSoft/LCFIVertex/pull/4))
  - Split processors into their own library as some tools are used in LCFIPlus, ilcsoft/LCFIVertex#2
  - Drop old boost tarball, rely in cvmfs/system installation of boost that is now needed elsewhere as well, ilcsoft/LCFIVertex#3
  - TState: replace usage of gear with MarlinUtil to get bfield value

* 2017-06-15 Andre Sailer ([PR#1](https://github.com/iLCSoft/LCFIVertex/pull/1))
  - Fixed warnings

# v00-07-03


# v00-07-01
* same code as  v00-07 but correct version number ( was  v00-06-02 in v00-07)
  
# v00-07
F.Gaede
* export external boost include directory in LCFIVertexConfig.cmake
* fixed installation of vertex_lcfi directory when build with external boost
* allow for using external boost 

# v00-06-02
* fixed duplicate namespace in forward declaration (F.Gaede)

# v00-06-01
* fixed cmake scripts for out-of-source installations

# v00-06
* made Ptcalc() a public static function so it can be used in analysis code
* fixed logical comparison bug in twotrackpid.cpp
* relaxed maximum number of iterations in vertexfuncmaxfinderclassicstepper.cpp
* added cmake configuration files in order to compile/link against LCFIVertex
* moved initialization of _dataSet variable in FlavourTag.cc to init() in order to fix problem of ProcessRunHeader() being skipped when SkipNEvents parameter is set in a steering file

# v00-05
* simplified cmake scripts (now uses the new ILCUTIL package )
* bug fix for building documentation in SL5.5

# v00-04-pre
tomohiko
* for gcc-4.4 compatibility: the boost library was updated to v1.44; other minor patches were added.

# v00-03-01
* same as v00-03- made compatible for MacOS

# v00-03
grimes
* src/RPCutProcessor.cc: Modified the Monte Carlo particle ID cut so that it can use a ReconstructedParticle to MCParticle relation collection as well as the current Track to MCParticle.  Also added a check on the type of objects the relation collection refers to and from.

engels
* CMakeLists.txt, LCFIVertexConfig.cmake.in: added 32 bit compatibility option

harderk
* src/TrueAngularJetFlavourProcessor.cc: prevent true jet flavour confusion due to anomalous PDG codes with more than five digits

walsh
* src/RPCutProcessor.cc: Remove ReconstuctedParticle V0s from selected Jets.
* src/NeuralNetTrainer.cc: bug fix: using D0Significance2 instead of Z0Significance in parts of the code.
* vertex_lcfi/src/lciointerface.cpp: line 165: cerr -> cout
* vertex_lcfi/src/lciointerface.cpp: Bug fix: only RPs with a track are used to make a Jet in jetFromLCIORP.

jeffery
* doc/zvtop.pdf, vertex_lcfi/zvtop/include/doc.h: Put zvtop docs in pdf

devetak
* src/TrueAngularJetFlavourProcessor.cc: Bug... repetition is bad but errors are worse...

harderk
* src/ConversionTagger.cc: remove LCCollection::setSubset() call that seems to cause problems with re-reading the V0Veto collection in separate Marlin jobs
* diagnostics/: include/V0Performance.h, src/V0Performance.cc: improve treatment of SimTrackerHit count
* src/ConversionTagger.cc: adjust default values and histogram ranges according to cut tuning results so far
* include/ConversionTagger.h, src/ConversionTagger.cc: replace entire vertex reconstruction by new MarlinUtil HelixClass code
* src/ConversionTagger.cc: improved (3d) calculation of helix distances

# v00-02-06-dev 
harderk
* src/ConversionTagger.cc: clean up hardcoded numbers such as particle masses
* diagnostics/src/V0Performance.cc: fix bug that would cause segmentation fault when dealing with corrupted input collections
lastovic
* vertex_lcfi/src/TState.cpp: Magnetic field taken from Gear.

# v00-02-05-dev
walsh
* src/DSTCollectionProcessor.cc:- Put try-catch blocks in the correct order;- Define pNames.size() and mcNames.size() as integers in the initialisation loops of fv and mcv; - FTSelectedJets -> Jets as the default for _JetCollectionName

martin
* src/RPCutProcessor.cc: Comment out warning line: " MCParticles related to this track! No cuts performed on MC PDG code. which appears too much to act as any warning!
* src/LCFIAIDAPlotProcessor.cc: Bug fix in purity printout for BC tag.	Thanks to Sonja for spotting this!

harderk
* src/ConversionTagger.cc: exploit momentum asymmetry between Lambda decay proton and pion
* diagnostics/src/V0Performance.cc: another Lambda plot
* diagnostics/src/V0Performance.cc: add histogram with momentum distributions of Lambda decay products

# v00-02-04-dev
lynch
* include/DSTAIDAPlotProcessor.h, include/DSTCollectionProcessor.h, 
* src/DSTAIDAPlotProcessor.cc, src/DSTCollectionProcessor.cc,
* src/DSTPlotProcessor.cc, steering/ftfordst.xml: changed DST
* Processors so that PID info is added to the existing collection added DSTAIDAPlotProcessor

harderk

* include/ConversionTagger.h, src/ConversionTagger.cc: add capability to enable or disable tagging of conv,K0,Lambda0 separately
* include/ConversionTagger.h, src/ConversionTagger.cc: reject candidates coming from vicinity of the IP
* src/ConversionTagger.cc: bug fix: existing conv/V0 candidates (e.g. found by PandoraPFA prior to running ConversionTagger) were not added to the LCCollection containing all conv/V0 candidates.
* diagnostics/: include/V0Performance.h, src/V0Performance.cc: count and print total number of conv/V0 candidates separately for conv, K0, Lambda. this improves diagnostics to tune ConversionTagger cuts separately for the different candidate classes
* diagnostics/: include/V0Performance.h, src/V0Performance.cc: add capability to look for conversions only in selected LCCollections. remove some debug output
* src/ConversionTagger.cc: reduce verbosity a bit

# v00-02-03-dev
2nd release for full-chain tests for ILD mass reconstruction.
Changes to the way the flavour tag collections are written out for DST

* TrueAngularJetProcessor + NeuralNetTrainer + PlotProcessor + AIDAPlotProcessor: changed so that the MC Truth info is stored in one LCFloatVec, in the order TrueJetFlavour, TruePDGCode, TrueHadronCharge, TruePartonCharge (Erik Devetak, Victoria Martin)
    
* 2 new processors created in order to add FT and MC Truth info to the jets using ParticleID objects This will reduce the number of collections for the mass reconstruction DSTs:

* DSTCollectionProcessor : adds 2 ParticleID algorithms to each jet: "LCFIFlavourTag": contains the flavour tagging info,  "MCTruth": contains the Monte Carlo truth info

* DSTPlotProcessor : a rewrite of PlotProcessor, which uses the DST collections  from the ParticleID objects instead of the full collections. also contains a function that will compare the parameters from both collection types and print out a message if they are different. (Clare Lynch)

# v00-02-03-pre
release for full-chain tests for ILD mass reconstruction, fixing a number of small issues found since last release and adding further functionality:
 
* steering files updated to include steering file running FullLDCTracking, PandoraPFA (instead of cheattracks+jetfind); based on current version of steering file for mass reconstruction from ILD mailing list (Clare Lynch)
 
* LCFIAIDAPlotProcessor: now compatible with RAIDA;    (Victoria Martin) corrected: "ptmasscorrection" --> "ptcorrectedmass"  (Roberval Walsh)
 
* RPCutProcessor: small bug fix: this processor should not require GEAR unless hadronic interaction removal is selected  (Kristian Harder)
 
* ZVRES: small bug fix in LCFIVertex/vertex_lcfi/zvtop/src/vertexfinderclassic.cpp to ensure consistency with ZVTOP paper  (Ben Jeffery)
 
* ZVKIN: return IP if no vertices found; sort found vertices wrt decay length from IP, in ascending order (Ben Jeffery)
 
* flavour tag inputs: code added to obtain new parameters for calculation of joint probability  from fit; currently only available if the package is run with JAIDA (Erik Devetak)

* neural network code:   minor fixes: ensured code compiles with gcc 4.2.1; (Kristian Harder, Dave Bailey) missing "const" added; at request of CALICE added class to determine the importance of the input variables  for a given data set and trained network based on the approach used in the TMVA package, see http://tmva.sourceforge.net/docu/TMVAUsersGuide.pdf p66 (Dave Bailey)

* NOTE: with this version, a new dependence on MarlinUtil is introduced; this became necessary for the ConversionTagger processor, currently under development- don't use this processor yet!  (Kristian Harder)

# v00-02-02
* LCFIVertex/vertex_lcfi/algo/src/pereventipfitter.cpp changed to use Kalman filter for vertex fitting (reduces run time for PerEventIPFitter.cc by about a factor 100) (Kalman filter developed by S. Gorbunov, I. Kisel, Vertex Package interface: Tomas Lastovicka)

* RPCutProcessor: for suppression of hadronic interactions using MC information this  processor is now updated to support arbitrary detector geometries (Kristian Harder)

* LCFIAIDAPlotProcessor: new processor providing diagnostics for flavour tag (Victoria Martin) Please note -BUILD_WITH="AIDAJNI" must be specified to cmake to build LCFIAIDAPlotProcessor
* LCFIVertex/macro/FlavourTagInputsOverlay.py: python script to create plots from output of LCFIAIDAPlotProcessor (Victoria Martin)
* LCFIVertex/macro/MakePurityVsEfficiencyRootPlot.C: root macro to plot purity vs efficiency graphs that can be written out by PlotProcessor, if compiling with root option; permits results from different runs (e.g. with different settings) to be compared (Kristian Harder) 

* vertex charge calculation: correction of algorithm to agree with the procedure developed with fast MC (Erik Devetak)
* separate VertexChargeProcessor, storing results in separate LCFloatVec collections; this used to be done in FlavourTagInputsProcessor (Erik Devetak)

# v00-01-01
* same code as v00-01 the first official release except that cmake build is supported