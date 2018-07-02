# DivideTools
Set of matlab based tools for analyzing drainage divide stability built on top of TopoToolbox, refer to paper published in Earth and Planetary Science Letters for additional details: https://www.sciencedirect.com/science/article/pii/S0012821X18302292. If you use these functions for a publication, please cite the linked paper.

# Using DivideTools
The starting point for analysis is running DivideStability. The other functions in this repository generally require inputs from DivideStability to run.

# Errors Relating to GRIDobj Instatances Not Aligning
If you get an error when running AcrossDivide that two GRIDobjs do not align, the most likely cause of this is that your reprojected your DEM before supplying it to DivideStability but did not specify a cellsize and instead let your GIS program of choice choose a cellsize. This usually results in a cellsize with a long string of decimals that can be rounded differently during different processing steps. To avoid this error, it is strongly encouraged that when reprojecting your DEM that you set a cellsize to a whole number. You can do this directly in Matlab using the resample function included with TopoToolbox, for example:

DEM=resample(DEM,round(DEM.cellsize),'bicubic');
