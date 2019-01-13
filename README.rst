Radial Velocity Follow-up Calculator
====================================

This repository contains python code to calculate the number of radial velocity measurements required to detect the semi-amplitude of a transiting exoplanet at a desired detection significance. 

These calculations can be performed in the general case of correlated RV noise and in the white noise limit which is less computationally expensive but is often not an accurate approximation to real empirical RV time-series. Calculations of the number of RV measurements and total observing time are useful for planning RV follow-up observations and how to distribute observing time among a set of targeted transiting planetary systems. 

Additionally, a web-based tool to perform equivalent calculations can be found `here <http://maestria.astro.umontreal.ca/rvfc/>`_. 

The model and calculator are described in detail in this `paper <https://arxiv.org/abs/1807.01263>`_.

Basic Usage
-----------

Begin by cloning this repository and navigating to its base directory via::

	git clone https://github.com/r-cloutier/RVFollowupCalculator.git
	cd RVFollowupCalculator

Because of the numerous parameters required to run calculations with the RVFC (i.e. planetary, stellar, spectrograph, RV noise sources), its basic usage is best demonstrated with an example given in this `ipython notebook <https://github.com/r-cloutier/RVFollowupCalculator/blob/master/example_nRV_calculation.ipynb>`_ which can be opened in its interactive mode via::

	jupyter notebook example_nRV_calculation.ipynb

Alternatively, if the user does not have `jupyter <https://jupyter.org/>`_ installed, use::

	ipython notebook example_nRV_calculation.ipynb

Citation
--------

If you use this code please cite the following paper::

  @ARTICLE{2018AJ....156...82C,
  author = {{Cloutier}, Ryan and {Doyon}, Ren{\'e} and {Bouchy}, Francois and {H{\'e}brard}, Guillaume},
  title = "{Quantifying the Observational Effort Required for the Radial Velocity Characterization of TESS Planets}",
  journal = {\aj},
  keywords = {methods: analytical, planets and satellites: detection, planets and satellites: fundamental parameters, techniques: radial velocities, Astrophysics - Earth and Planetary Astrophysics},
  year = 2018,
  month = Aug,
  volume = {156},
  eid = {82},
  pages = {82},
  doi = {10.3847/1538-3881/aacea9},
  archivePrefix = {arXiv},
  eprint = {1807.01263},
  primaryClass = {astro-ph.EP},
  adsurl = {https://ui.adsabs.harvard.edu/\#abs/2018AJ....156...82C},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }
