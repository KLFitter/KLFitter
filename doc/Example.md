**DESCRIPTION OF STANDALONE EXAMPLE. SHOULD PROVIDE LINK TO [KLFitterExtras](https://gitlab.cern.ch/KLFitter/KLFitterExtras). THE KLFitterExtras README.md SHOULD POINT BACK TO THE KLFitter README.md AND THIS FILE EXPLICITLY. I HAVE ALSO MOVED HERE THE DESCRIPTION OF THE EXAMPLE-SPECIFIC OUTPUT FROM WhatIsKLF.md.**



##### Output

The output in the provided examples comes in form of a ROOT file which
contains the following trees:

* `TreeMeasured`: the measured particles.
* `TreeSelected`: the particles selected for the fitting.
* `TreeModel`: the model particles (e.g., including also the top quark
  and other intermediate particles).
* `TreeTruth`: the Monte Carlo truth information (if avialable).
* `TreeMap`: maps containing the mapping between the measured and the
  selected particles.

The TreeModel contains the results from the fit. The most important
stored variables are:

* `EventNumber`: the event number.
* `N_permutations`: the number of permutations.
* `best_permutation[N_permutations]`: a list of indices ordered by the
  EventProbability. The best permutation has the index
  best_permutation[0], etc.
* `EventProbability[N_permutations]`: the event probability for each
  permutation.
* `LogLikelihood[N_permutations]`: the log Likelihood for each
  permutation.
* `ConvergenceStatusBit[N_permutations]`: the convergence status bit
  word for each permutation.
* `x_E[N_permutations]`, `x_pt`: x is a particle known in the model
  (i.e., initial state or intermediate particle, etc.). The `_E`, ...
  describe the kinematics for each permutation.
* `par_x[N_permutations]`: x is the best fit parameter for each
  permutation.
* `parerr_x[N_permutations]`: x is the estimated uncertainty on the
  best fit parameter for each permutation.
