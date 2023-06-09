# locust-dev
Locally Optimized Clustering for Unsteady Structure Tracking

Hi!
This repository contains sample and test code for the unsteady feature tracking algorithm associated with the following manuscript:

	Rowell, Colin R., Jellinek, A. Mark, Gilchrist, Johan T., (in review) Tracking volcanic plume thermal evolution and eruption source unsteadiness in ground-based thermal imagery using spectral-clustering.

This repository contains only a subset of all scripts and function used in the above manuscript. The contents focus on essential input and test data sets for the core elements of the workflow, specifically feature tracking and atmospheric profile removal. Other components of the workflow such as pre-processing will be added as needed and as time allows. The first section below, a summary of the workflow and major functions also includes information on which components of are included here. The next section after that also lists some important data and code dependencies that will be required to run certain elements of the workflow, such as the Optical Flow Analysis Toolbox (Sun et al., 2014) and sample data sets that are too large to include in this repository.



=============== CONTENTS ===============
	(1) Data availability, data and code dependencies
	()	Summary of workflow and major functions (what's included, directory structure)
	()	Data and code dependencies
	()	Short descriptions of core functions
	()	Notes for using the feature tracking algorithm


=============== DATA AVAILABILITY AND DEPENDENCIES ===============

Curated test data:

Main dataset for the manuscript:



=============== SUMMARY OF FULL WORKFLOW AND MAJOR FUNCTIONS ===============
--> WHAT'S INCLUDED SO FAR

1 - fully included and ready for use, test data included in repository
T - included and ready for use, test data available seperately
! - available, subject to dependencies
x - not currently included

x	(1) 	Irbis to Matlab data conversion
x	(2) 	Image registration (stabilization)
x	(3) 	Image projection mapping
x	(4) 	Plume masking (modified from plumeTracker, Bombrun et al., 2018)
x	(5.1) 	Re-gridding image frames (x,z)
x 	(5.2) 	Re-sample gridded frames 3D data cube, regularly spaced in time
x 	(6) 	Get 2D velocity fields with Optical Flow Analysis
x 	(7) 	Atmospheric profile fitting and removal
x 	(8) 	Column source time-series retrieval
x 	(9) 	Create time-averaged images
x 	(10) 	Feature tracking
x 	(11) 	Virtual source estimation and power-law fitting

--> DIRECTORY STRUCTURE

sabancayaScripts/
	Contains project input/control scripts and various calculations specific to the manuscript. 
		- "project" files are useful references for all workflow input/output. 
		- "structureTracking" files contain input parameters used for all tracked column structures as part of the manuscript. Subsets of these are highlighted that apply to provided test data sets
		- ""

=============== NOTES FOR USING PULSETRACKER ALGORITHM ===============

Parameters 'memoryN' and 'uMax' are automatically calculated (and reported in command line output) when running pulseTracker, but the process is somewhat slow. If running the tracker repeatedly to test tracks, it is best to specify these two values as input after their initial calculation in order skip calculation time in succeeding runs.

COMMON ISSUES WITH pulseTrack

The tracked cluster grows too large to encompass more than the feature of interest.
	-> Try increasing any of 'Tpercentile', 'minClust','maxClust', or slightly reducing 'uTol','winSzRatio'

The tracked cluster shrinks or is too small to cover the feature of interest
	-> Try decreasing any of 'Tpercentile', 'minClust','maxClust', or slightly increasing 'uTol','winSzRatio'

The tracked cluster wanders off the feature of interest to track something else.
	-> Try increasing the prior regularization 'lambda'
	-> Try specifying/changing the initial prior mask to better target the feature of interest

The tracked cluster lags behind the feature of interest, eventually falling off entirely.
	-> Check/adjust 'uTol' - this is intimately linked with the other internal parameters 'uMax' and 'pxTol' by the relation (pxTol = uMax*uTol*dt/dz), where dt and dz are the frame time step (in seconds) and pixel grid spacing (in meters), respectively. 'pxTol' is always allowed to be at least 1, but 2 is often better.
    -> Check the above parameters as well!
    -> Try adjusting the clusterWeights or priorWeights if need be.
    clusterWeights = [X Z Temp. U W]
    priorWeights   = [Temp. W Distance Area]
