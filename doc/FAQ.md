# Frequently Asked Questions

## Problems during the installation

Please note that we strongly encourage people to use the cmake installation
procedures described in the [README
section](../README.md#installation-via-cmake-recommended). Installation via
Makefile has only been tested on few systems, and might not be cross-platform
compatible.

#### Why is ROOT not found by cmake?

Include the ROOT binary path to the `PATH` variable and make sure it is
exported. Only then it can be picked up by cmake correctly. As implemented in
[FindROOT.cmake](../cmake/FindROOT.cmake), the ROOT location is determined by
cmake by searching for the `root-config` executable.

#### Why is my local BAT version not found by cmake?

The KLFitter cmake configuration allows an integrated build of the BAT library,
as described in the [README](../README.md#installation-via-cmake-recommended).
If you use a local version of BAT (no download via cmake) and it cannot be
found, make sure to set the `$BATINSTALLDIR` variable and export it. This
variable is used by cmake to locate the library.

#### What to do when linkage against ROOT/BAT fails ...

We have encountered problems with the compiler versions when linking against
ROOT or BAT, in particular if any of the two libraries was compiled with a
different version than the one cmake uses. If you set up a custom compiler for
the build process, maybe cmake doesn't pick this compiler up correctly. This can
be tested by checking whether `whereis c++` and `which c++` point to the same
binary. To tell cmake to use the latter one, you can export the compiler
location:

```sh
> export CXX=`which c++`
> export CC=`which gcc`
```

## Questions/problems with the implementation

#### What is the difference between the top mass parameter and the hadronic/leptonic top mass?

This question refers to `LikelihoodTopLeptonjets`: The top mass parameter is the
central value of the top mass Breit-Wigner distribution. The hadronic and
leptonic top masses are the invariant masses of the three particles from the
hadronic and leptonic top decay, respectively. The hadronic and leptonic top
masses are assumed to be Breit-Wigner distributed around a central value, i.e.
the top mass parameter, in the definition of the likelihood.

#### In which units should I add objects to KLFitter?

KLFitter handles all four-vectors in GeV. This means also that all values in
`PhysicsConstants` need to be set in GeV.

#### Why do I get a message `WARNING : A particle with negative mass...`?

Particles that are added to the KLFitter must not have negative masses. Negative
masses lead to bad fits with `NaN` outputs. It is a feature that once such a bad
fit was performed, all fits in the following events also become bad.

#### Why are the fit results of all events starting from one event `NaN`?

Please check that you do not pass particles with negative masses to KLFitter.

#### What does KLFitter do if more than 2 jets are b-tagged and the b-tagging veto is used?

This question refers to `LikelihoodTopLeptonjets`, but it also applies to other
likelihoods where you could pass more b-tagged jets to KLFitter than b-quarks
are present in the model: With the "b-tagging veto" methods (such as `kVeto` or
`kVetoNoFit`), events with > 2 b-tags are not compatible with the underlying
model you are testing with KLFitter. Hence, all event probabilities will be zero
(`kVeto` and related) or all permutations will be removed, so that none remains
(`kVetoNoFit` and related).

#### Does KLFitter need at least one tagged b-jet if you use one of the b-tagging methods?

No.

#### Why do I get a `vector out of range` error when switching from an "e+jets" event to a "mu+jets" event (or vice versa)?

This question refers to `LikelihoodTopLeptonjets` (and related likelihoods):
`LikelihoodTopLeptonJets` is testing a hypothesis for the underlying model,
which can either be "e+jets" or "mu+jets". You can set this via
`LikelihoodTopLeptonJets::SetLeptonType(LeptonType leptontype)`. The default is
`kElectron`. When switching from an "e+jets" to a "mu+jets" event (or vice
versa), you must also change the underlying model. If you do not do this,
`LikelihoodTopLeptonJets` will ask for the first lepton of the other type and
will throw a `vector out of range` exception.

#### Why is the number of permutations not reduced when I use the b-tagging veto `kVeto`?

In this option, we do not want to remove permutations due to the b-tagging
approach. Only the event probability is affected by weighting it with either 0
or 1, but not the number of permutations, although b-tags veto certain
permutations. To reduce the number of permutations use `kVetoNoFit`.

#### Why do I get a `index out of range` error when I use the `kVetoNoFit` b-tagging veto?

This is most likely due to events with more than 2 b-tagged jets.

#### What is the number of permutations in `LikelihoodTopLeptonsJet` if I pass n jets to KLFitter?

The number of permutations, N, is:

N = (n!) / (2 * (n-4)!)

where n is the number of considered jets. The factor 2 accounts for the
indistinguishable light quarks and the factor (n-4)! removes all invariant
permutations of the jets that are not used in the fit.

#### I have 2520 permutations in one event. What did I do wrong?

This answer is specific to `LikelihoodTopLeptonjets`, but it applies to other
likelihoods in a similar way: I means that you consider by far too many (in this
example 10) jets in the kinematic fit, i.e. you probably pass all jets in the
event to KLFitter. Besides of the calculation of the permutation table (for 2520
permutation) taking a lot of time, it most probably does not make sense to fit
more than 5 or 6 jets. With more jets passed to KLFitter, it is very likely that
by chance you will find any kind of event with a lepton and missing transverse
momentum matching somehow the lepton+jet decay topology and an identification of
the correct permutation is unlikely.

#### How do I improve speed when using large numbers of jets as input to KLFitter?

Say, you want to use 6 jets in your analysis, but still want to consider more
jets within the possible permutations, say up to 8. Then, you can give an
additional argument to `KLFitter::SetParticles(KLFitter::Particles * particles,
int nPartonsInPermutations = -1)` and set `nPartonsInPermutations` to 6, while
you pass up to 8 jets to KLFitter. The calculation of the permutations will not
take into account all 8! jet permutations, but efficiently only calculate those
needed.
