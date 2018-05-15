# Denoising_GCVwavaletF
=======================

This repository contains MATLAB scripts and sample data for applying denoising method presented in: 
"Mousavi, S. M., and C. A. Langston (2017). Automatic Noise-Removal/Signal-Removal Based on the
General-Cross-Validation Thresholding in Synchrosqueezed domains, and its application on earthquake data,
Geophysics.82(4), V211-V227 doi: 10.1190/geo2016-0433.1"

`demo.m` includes all info you need to know for running the code 

## Paper
(https://www.researchgate.net/publication/315917772_Automatic_noise-removalsignal-removal_based_on_general_cross-validation_thresholding_in_synchrosqueezed_domain_and_its_application_on_earthquake_data)

## Talk 
(https://earthquake.usgs.gov/contactus/menlo/seminars/1093)

## Abstract 
Recorded seismic signals are often corrupted by noise. An automatic noise attenuation method for single-channel seismic data is presented, based upon high-resolution time-frequency analysis. Synchrosqueezing is a time-frequency reassignment method aimed at sharpening a time–frequency picture. Noise can be distinguished from the signal and attenuated more easily in this reassigned domain. The threshold level is estimated using a general cross validation approach that does not rely on any prior knowledge about the noise level. Efficiency of thresholding has been improved by adding a pre-processing step based on Kurtosis measurement and a post-processing step based on adaptive hard-thresholding. The proposed algorithm can either attenuate the noise (either white or colored) keeping the signal or remove the signal and keep the noise. Hence, it can be used in either normal denoising applications or pre-processing in ambient noise studies. We test the performance of the proposed method on synthetic, microseismic, and earthquake seismograms.

## A Short Description 
Seismic data recorded by surface arrays are often contaminated by unwanted noise. In many conventional seismic methods, 
the reliability of the seismic data and accuracy of parameter extraction, such as onset time, polarity, and amplitude, 
are directly affected by the background noise level. As a result, the accuracy of event location and other attributes 
derived from seismic traces are also influenced by the noise content. Therefore, there is a great need for developing 
suitable procedures that improve signal-to-noise ratios allowing for robust seismic processing. In this presentation, 
I introduce four different methods for automatic denoising of seismic data. These methods are based on the time-frequency 
thresholding approach. The efficiency and performance of the thresholding-based method for seismic data have been improved 
significantly. Proposed methods are automatic and data driven in the sense that all the filter parameters for denoising are 
dynamically adjusted to the characteristics of the signal and noise. These algorithms are applied to single channel data 
analysis and do not require large arrays of seismometers or coherency of arrivals across an array. Hence, they can be applied
to every type of seismic data and can be combined with other array based methods. Results show these methods can improve 
detection of small magnitude events and accuracy of arrival time picking.

![This figure illustrates the process of adding noise to the synthetic seismogram](Fig1.png) 
This figure illustrates the process of adding noise to the synthetic seismogram

![Results on synthetic seismogram](Fig2.png) 
Results on synthetic seismogram

![Results on real seismogram](Fig3.png)
Results on real seismogram


