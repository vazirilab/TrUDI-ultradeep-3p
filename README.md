# TrUDI-ultradeep-3p
Analysis code for TrUDI ultra-deep three-photon imaging


## System requirements 

* OS: Windows 10 or later; macOS 15 or later
 
* MATLAB: 2024b or later 

* MATLAB toolbox: Image Processing Toolbox, Parallel Computing Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox, Cellpose for Microscopy Segmentation

* Demo tested on a workstation running Windows 11 with a 24-core CPU, 128 GB RAM. Installed MATLAB version 2024b with all packages above. 

### External libraries

* CaImAn-MATLAB: https://github.com/flatironinstitute/CaImAn-MATLAB

* NoRMCorre: https://github.com/flatironinstitute/NoRMCorre

* Cellpose MATLAB: https://www.mathworks.com/help/medical-imaging/ref/cellpose.html

* ScanImageTiffReader (OS dependent): https://docs.scanimage.org/Appendix/ScanImage%2BTiffReader.html


## Installation

Extract the code in a separate folder. 

Ensure the extracted folder dependencies above.

Add the paths to the libraries and functions located in the External and Utilities folders.

Typical installation time: < 1 min.


## Brief description of main processing files

* preprocess_scanimage_rawdata.m – Reads and prepares raw data from ScanImage for subsequent pipeline processing.

* motion_correction_batch.m – Performs motion correction.

* motion_correction_plots.m – Visualizes motion analysis and corrected images.

* caiman_roi_detection_batch.m – Runs CaImAn ROI detection.

* caiman_roi_refine_batch.m – Refines detected ROIs from CaImAn analysis based on specified ROI criteria.

* caiman_results_plots.m – Visualizes CaImAn analysis traces and ROIs.

* cellpose_segmentation_stack_data.m – Performs anatomical segmentation of stack data using Cellpose.

* cellpose_segmentation_refine_results.m – Refines results from Cellpose analysis.

* noise_estimate_from_data.m – Estimates noise in the experimental data and fits it with a noise model.


## Demo usage

The demo data contains a small experimental dataset. The demo_run.m files demonstrate CaImAn analysis steps. It also adds all the necessary paths to the included packages. Each step should take less than 1 min. 

1. Run preprocess_scanimage_rawdata.m to prepare raw data with default parameters. Expected output: New subfolder with pre-processed data as a TIF file.

2. Run motion_correction_plots.m to visualize motion analysis and corrected images. Expected output: Figures with motion-corrected images and a plot of the motion correction metric.

3. Run motion_correction_batch.m to apply motion correction with default parameters. Expected output: New subfolder with motion-corrected data as TIF and MAT file.

4. Run caiman_roi_detection_batch.m to identify active ROIs on motion-corrected TIF data folder. Expected output: New subfolder with detected ROIs data for the next step.

5. Run caiman_roi_refine_batch.m to refine detected ROIs from CaImAn analysis based on specified ROI criteria on the detected ROI data folder. Expected output: New subfolder with refined ROIs data.

6. Run caiman_results_plots.m to visualize CaImAn analysis traces and ROIs. Expected output: Figures with ROI overlay on the image and time-trace plots.

