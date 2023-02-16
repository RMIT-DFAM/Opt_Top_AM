# Opt_Top_AM
Optimal Topology for Additive Manufacture

Introduction
This package of scripts and functions builds a planar optimal geometry and modifies to be manufacturable with additive methods. A full description is available in the publication:
Optimal topology for additive manufacture: A method for enabling additive manufacture of support-free optimal structures [Int J. Materials and Design 63 (2014) 678–690]

Scripts & Functions Include:
Topology Optimisation & AM Orientation
-	topmain.m : Main code (THIS IS RAN TO OBTAIN RESULTS)
-	top88.m: 88 Line topology optimisation code
-	top.m: 99 Line topology optimisation code by Ole Sigmund with update
-	topCantiDist.m: 99 Line topology by Ole Sigmund
-	smoothedge.m: Defines constrained and loaded edges from a cell array of checkered
boundaries and return a cell array with smoothed boundaries, except for
constrained and loaded ones.
-	rotopt.m: This function takes user defined inputs for initial, step & final orientation and rotates part to be manufactured, which then allows user to see optimal orientation for minimal support structure during additive manufacture.
-	nodes.m: This function is used to create boundary node coordinate from a binary matrix input starting from the south west corner of the ‘2D Planar Grid’.
-	igesout.m:  .igs Converter for points, lines and Nurbs curves and surfaces
-	gcolor.m: Function creates gradient plots for structure based on user defined feasible angle of manufacture
-	ConvertToBinary.m: Converts a full matrix of elements into a binary matrix, by comparing every value with user given threshold and print a grayscale colour map of the binary matrix.
-	Comma2point.m: replaces all occurrences of comma (",") with point (".") in a text-file.
-	Buildextsupport.m: This function builds support structures for external boundaries.
Additional functions:
All other code comes from ‘addaxis’ by Harry Lee available on Mathworks File Exchange, link : (https://au.mathworks.com/matlabcentral/fileexchange/9016-addaxis )
It is used to add multiple axis to a plot.
Instructions
All that is required is running the topmain.m script; user will be prompted to input variables such as number of elements in the x & y direction, volume fraction etc.





