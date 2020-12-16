# Hemodynamic-response-estimation-from-fNIRS

Accurate Hemodynamic Response Estimation by Removal of Stimulus-Evoked Superficial Response in fNIRS Signals.

The code addresses the problem of hemodynamic response estimation when task-evoked extra-cerebral components are present in functional near-infrared spectroscopy (fNIRS) signals. These components might bias the hemodynamic response estimation. Therefore careful and accurate denoising of data is needed.
We propose a dictionary-based algorithm to process every single event-related segment of the acquired signal for both long separation and short separation channels. Stimulus-evoked components and physiological noise are modeled through two distinct waveform dictionaries (contained in the .mat file "Dictionaries concentration").
For each segment, after the removal of the physiological noise component in each channel, a template is employed to estimate stimulus-evoked responses in both channels. Then, the estimate from the short-separation channel is employed to correct for the evoked superficial response and refine the hemodynamic response estimate from the long-separation channel.

To execute the code run "main.m". 
