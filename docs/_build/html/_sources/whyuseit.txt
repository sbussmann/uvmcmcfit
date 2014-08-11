Why Use ``uvmcmcfit``?
**********************

Here are three reasons why you might be interested in using ``uvmcmcfit``.

Markov Chain Monte Carlo
------------------------

``uvmcmcfit`` uses D. Foreman-Mackey's `emcee
<http://dan.iel.fm/emcee/current/>`_ to do Markov Chain Monte Carlo (MCMC).
`emcee <http://dan.iel.fm/emcee/current/>`_ is an affine-invariant MCMC
ensemble sampler that is efficient at identifying the best-fit model and the
uncertainties on the parameters of the model. 

Interferometry
--------------

Interferometers measure visibilities, not images.  In many situations (e.g.,
strong gravitational lensing), it is advantageous to measure the goodness of
fit from the visibilities rather than the deconvolved images.  

``uvmcmcfit`` includes a pure-Python adaptation of Miriad's ``uvmodel`` task
that allows it to generate simulated visibilities given observed visibilities
and a model image.

Gravitational Lensing
---------------------

Many of the brightest galaxies that have been observed by interferometers are
gravitationally lensed by an intervening galaxy or group of galaxies along the
line of sight.  Understanding the intrinsic properties of the lensed galaxies
requires a lens model.

``uvmcmcfit`` includes a simple ray-tracing routine that allows it to account
for both strongly lensed systems (where multiple images of the lensed galaxy
are detected) and weakly lensed systems (where only a single image of the
lensed galaxy is detected).

