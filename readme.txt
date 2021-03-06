###################################################################                                                 
#    fast Robust Matrix Completion (fRMC) for background subtraction                             
#
#    Behnaz Rezaei (brezaei@ece.neu.edu)
#    Augment Cognition Lab (ACLab) http://www.northeastern.edu/ostadabbas/
#    Copyright, 2017                        
#                                                     
###################################################################

1. Introduction

If using the fast Robust Matrix Completion code, please cite the following work in any resulting publication:

@inproceedings{rezaei2017background,
  title={Background Subtraction via Fast Robust Matrix Completion},
  author={Rezaei, Behnaz and Ostadabbas, Sarah},
  booktitle={Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition},
  pages={1871--1879},
  year={2017}
}

###################################################################

2. Installation

This code is written for Matlab (tested with versions 2016b, 2017a) and requires the Matlab Image Processing and Computer Vision Toolbox. 

###################################################################

3. Getting Started

The code is quite compact, as is the method. Start with Demo_fRMC.m which applies the fRMC on a sample video in Data folder and finally saves the example results of background subtraction as a video in the Data/result_fRMC folder as well as the elapsed time. It saves the background and foreground as .mat data so you can use then for visualization or other purposes. You may also want to update the morphological post-processing for better results for your own video. If the video is long you would better to run this algorithm block-wise, the recommended video size is 500 frames.
This algorithm is verified on the two datasets: SABS2012(http://www.vis.uni-stuttgart.de/en/research/information-visualisation-and-visual-analytics/visual-analytics-of-video-data/sabs.html),
                                                BMC2012(http://bmc.iut-auvergne.com/)

###################################################################

4. Contents

Code:
   infaceExtFrankWolf	- Apply in-face extended frank-wolf algorithm for low-rank matrix completion. if you are not sure about gamma1 and gamma2 just give [] as argument.
   Demo_fRMC      	- Demo demonstrating fRMC on a sample video from BMC2012 dataset 

Other:
    readme.txt  	- This file.

###################################################################

5. History / ToDo

Version 1.0 (19/11/2017)
 - initial version
