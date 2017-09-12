# ADEMS code v2

Changes from previous version:…

The tracking software consists of a series of functions for opening data and tracking particles, as well as a number of analysis functions.

Useful to know: There is a ‘show_output’ option that can be used to view graphs and manually advance at each stage of the algorithm.

Cursor_mode: set =1 and user can manually specify where spots are. The code will return intensity values at this point over the whole time series.

**FinalCode12**: This is the main tracking program and is a function of image_label. It returns spot arrays for each channel and frame_average. It uses:

* **extractImageSequence**: extracts user set frames from tif specified by image_label
* **LaserOn2**: Calculates where the first illuminated frame is based on maximum intensity
* **FrameAverage**: Calculates a frame average/summation over set no. frames
* **findSpots3**: thresholds the image with Otsu’s method to find candidate spots
* **findSpotCentre2**: performs iterative Gaussian masking to find spot centre and total intensity
* **fit2DgaussianFixedCenter**: fits a constrained 2D Gaussian to find sigma_x/y and central intensity
* **MergeCoincidentSpots3**:  iteratively finds spots which are too close together and averages their centres, assigns junk values to now redundant spot numbers
* **LinkSpots3**: links spots into trajectories based on proximity

**FinalCode12Preloaded**: Is a version of final code for when you have already loaded the image data (for example if you have opened and manipulated it already)

**TrackAll5**: loops over date stored in folder hierarchy: strain\date\sample\cell(field of view) and tracks largest tif (over certain threshold). Saves data to two analysis directories. Runs in parallel with error handling

**IntensityPlot**: generates plots from the data file and frame average

**MoviePlot2**: reads in spot data and tif data and plots trajectories on each frame to make into a movie

**SpotBaseline**: reads in spot data and image data,  uses last x/y co-ords of each trajectory as baseline cords and calculates 10 frames intensity values, stores it back in a new spot array. Use this to calculate zero values for seeing step bleaching.

**Segmentation**: scripts for segmenting cells

**Useful Scripts**: various scripts including kernel density function plots, plotting reconstructions, hot plots, cluster analysis, chung kennedy filtering etc,

SpotsCh1/2 is an array, each row contains the information for a spot found in an image frame in the series. The columns contain the following information:

1. X coordinate (pixels)
2. Y coordinate (pixels)
3. Clipping_flag (a switch, please  ignore)
4. Mean local background pixel intensity
5. Total spot intensity, background corrected
6. X spot sigma width
7. Y spot sigma width
8. Spot central intensity (peak intensity in a fitted Gaussian)
9. Frame number the spot was found in
10. Trajectory number, spots in the same trajectory have the same trajectory number
11. Signal to noise ratio
12. Frame in which laser exposure began

