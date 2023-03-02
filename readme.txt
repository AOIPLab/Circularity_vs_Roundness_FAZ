Circularity vs. Roundness for real and simulated FAZs
Created by: Jenna Grieshop
Date: 6/17/2020 - 9/3/2020
MCW AOIP Dr. Carroll

This program creates simulated circles, ellipses, and polygons that can
be measured with methods that an FAZ would be. It also reads in
actual FAZ segmentations to measure them. Various metrics are produced
including area, perimeter, roundness, circularity, etc.

This script was developed in MATLAB R2019b and was connected to ImageJ using MIJ/FIJI
To set up MIJ steps were followed from this site: http://bigwww.epfl.ch/sage/soft/mij/ under the Download and install heading
Lines 556 and 557 may need to be adjusted depending on your matlab path

This script uses functions ExactMinBoundCircle and max_inscribed_circle from Matlab File Exchange and are contained in the Circumscribed Circle Folder
(Tolga Birdal (2023). Maximum Inscribed Circle using Distance Transform (https://www.mathworks.com/matlabcentral/fileexchange/30805-maximum-inscribed-circle-using-distance-transform), MATLAB Central File Exchange. Retrieved February 28, 2023.)
(Anton Semechko (2023). Exact minimum bounding spheres and circles (https://github.com/AntonSemechko/Bounding-Spheres-And-Circles), GitHub. Retrieved February 28, 2023.)

Best practice is to run each of the sections as you need them individually (Actual FAZs, Circle Family, Ellipse Families, Hexagon, Octagon).
Always run the first section first to set up paths. Follow the instructions in the title of the file explorer.

Needed to run the Actual FAZ section:
	- .tif images of the FAZ
	- .csv file with the segmentation coordinates
	- a batch/LUT .csv formated with the subject ID, Eye, Axial Length, and Device Number


Function Descriptions:

Ellipse Family Function
 - elliFam(maj, min, AR1, AR2, folder)
 - maj: major axis
 - min: minor axis
 - AR1: first number of aspect ratio
 - AR2: second number of aspect ratio
 - folder: folder to store results in

Theoretical and Mathematical Calculations Function
 - calc(sides,shape,radius, major, minor, coords) *NOTE: if a parameter does not apply for shape, plug in 1)*
 	- sides: number of sides the shape has
 	- shape: number referencing the shape that is being analyzed (1 - circle, 2 - ellipse, 3 - hexagon, 4 - octagon, 5 - actual FAZ)
 	- radius: radius of circle; use as sidelength for polygons
 	- major: major axis of ellipse
 	- minor: minor axis of ellipse
 	- coords: coordinates for the shape

RegionProps and inscribed/circumscribed circles Function
 - RP(coords, filename)
 	- coords: coordinates for the shape
 	- filename: filename for the mask to be saved as

Connect to ImageJ Function
 - imageJ(shape,rNum,major,minor,sL,name) *NOTE: if a parameter does not apply for shape, plug in 1)*
 	- shape: number referencing the shape that is being analyzed (1 - circle, 2 - ellipse, 3 - hexagon, 4 - octagon, 5 - actual FAZ)
 	- rNum: radius
 	- major: major axis
 	- minor: minor axis
 	- sL: side length for polygons
 	- name: name of actual FAZ
