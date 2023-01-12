# TomographyAnimation

Matlab codes for creating tomographic simulation videos such as seen in the video "X-ray tomography 101: introduction" (https://youtu.be/7pdRq4XLT90).

Descriptions of functions:

Xray_parallel_plot.m
Creates a simulation showing parallel-beam tomographic measurement rotating around an object. The object is visible; the user can freely choose the object as a grayscale image of size 1024x1024.

Xray_parallel_mystery_plot.m
Same as Xray_parallel_plot.m, but hiding the object behing a round patch with a question mark.

Xray_parallel_BP_plot.m
Builds up (unfiltered) back-projection image from the tomographic data.

Xray_parallel_FBPrev_plot.m
Demonstrates filtered back-projection (FBP) reconstruction. The resulting video needs to be inverted in time. 
