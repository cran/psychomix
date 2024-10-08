# psychomix 1.1-9

* Achim Zeileis takes over maintenance from Hannah Frick.

* Updated example reference output for CRAN checks.


# psychomix 1.1-8

* Removed functionality relating to the `mRm` package which got removed
  from CRAN.


# psychomix 1.1-7

* Set `RNGversion("3.5.0")` for reproducibility of previous results due
  to the changes in `sample()` for R version 3.6.0.


# psychomix 1.1-6

* Load the `lattice` package in the vignettes in preparation to changes
  in the `flexmix` package.


# psychomix 1.1-5

* Conditional registration of `effect()` and `allEffects()` methods when `effects`
  package is loaded.


# psychomix 1.1-4

* Added a new function `mptmix()` for MPT mixture models, also known as
  latent-class MPT models. The design follows that of `raschmix()` and
  `btmix()`. Internally, the `flexmix` driver `FLXMCmpt()` is called to fit the
  finite mixture model. The function has not yet been fully tested and
  may change in future versions.


# psychomix 1.1-3

* Updated the example for `btmix()`.


# psychomix 1.1-2

* The manuscript "Rasch Mixture Models for DIF Detection: A Comparison
  of Old and New Score Specifications." has now been published in
  Educational and Psychological Measurement, 75(2), 208-234.
  [doi:10.1177/0013164414536183](https://doi.org/10.1177/0013164414536183).
  A preprint version is included as `vignette("scores", package = "psychomix")`.

* Improved functionality for `raschmix()` to allow for item response
  data in the `itemresp` class.

* Improved axis labeling in `plot()` method for `raschmix` objects.

* In the `plot()` method `nchar(..., type = "width")` is now used to determine
  the default abbreviation.

* If suggested packages are needed internally, these are only called
  with `::` semantics and not `require()`d anymore.


# psychomix 1.1-1

* Adapted `raschmix()` to work with both the old psychotools version
  0.2-0 and the new 0.3-0.
  
* Updated the `"scores"` vignette which is now also accepted for publication
  in Educational and Psychological Measurement.


# psychomix 1.1-0

* Improved functionality for `raschmix()` to allow for differences 
  between components in terms of identified parameters.

* Improved function `raschmix()` to leverage new functionality from 
  the `flexmix` package: Parameter estimates from the previous M-step 
  can now be used for initialization.

* Function `raschmix()` can now model the score distribution to be 
  equal across all components (for both a `"saturated"` and a `"meanvar"` 
  specification of the score model).


# psychomix 1.0-0

* Official first stable release of the `raschmix()` functionality in the
  package, accompanying the manuscript "Flexible Rasch Mixture Models
  with Package psychomix" by Frick, Strobl, Leisch, and Zeileis,
  published in the Journal of Statistical Software 48(7). See
  `citation("psychomix")` for details.


# psychomix 0.2-1

* For increased numerical stability the default minprior control
  parameter in `raschmix()` is now 0.05 (as in `flexmix`) and not 0
  (as in the previous `psychomix` version).

* Revised `vignette("raschmix", package = "psychomix")`. Specifically,
  there is a discussion of how the `FLXMCrasch()` can be used directly
  with `flexmix()` or `stepFlexmix()` from the `flexmix` package.

* Improved function `simRaschmix()` to allow for a flexible specification
  of the data generating process. 
  
* Added an `effectsplot()` function that leverages the `effects` package
  for visualizing the effects of the concomitant variables (if any) in
  the mixture model. This has not yet been fully tested and may change
  in future versions.
  
* Added a new function `btmix()` for Bradley-Terry mixture models. The
  design follows that of `raschmix()` rather closely. Based on `btReg.fit()`
  from package `psychotools`, there is a `flexmix` driver called `FLXMCbtreg()`.
  The `btmix()` function is a convenience interface calling `stepFlexmix()`
  with the `FLXMCbtreg()` driver. This has not yet been fully tested and
  may change in future versions.


# psychomix 0.1-1

* First CRAN release of new `psychomix` package for fitting
  psychometric mixture models based on `flexmix` infrastructure. At the
  moment only Rasch mixture models are implemented in various flavors:
  with/without concomitant variables, different parametrizations
  of the score distribution (saturated vs. mean/variance specification).
  See `vignette("raschmix", package = "psychomix")` for details.
