# XASanalysis-APSpumpprobe

These scripts process the raw data from pump-probe X-ray absorption experiments at the Argonne National Laboratory's Advanced Photon Source executed in 2017. The experiments were performed by the LSU group at EPFL, and relevant and resultant papers will be be linked on this page after publication.

A short description of the scripts:
- ProcessAPSdata : This script is the primary script for processing data. It loads raw data files and applies necessary corrections and normalizations to the data.
- produce_I0_averages : Called within ProcessAPSdata, this function produces averaged I_{0} values for X-ray absorption spectra in order to reduce introducing additional noise from I_{0} correction.
- normalize_data : Using a region of interest where data exists for all runs, this scales runs so they can be combined without increasing apparent noise.
- move_values : Rearranges columns of imported data to match the 2D size and arrangement of the preceding loaded runs.
- extract_detectors : Parses the raw data files to make a list of the detectors operating during a run, which may change as hardware is changed during the experiment.
- find_detector_indices : Uses detector labels from raw files to determine the identity of data columns which may change from one scan to another.
- correct_baseline : Normalizes to zero all scan points that should result in zero X-ray fluorescence intensity, such as before time zero or far in the pre-edge spectrum.
monochromator_lowpass : Not currently used, but is designed to improve noise that is significantly higher than the resolution of the monochromator.
