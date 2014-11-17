MEA-toolbox
===========

2011-2014 Aleksander Domanski PhD (adomansk@exseed.ed.ac.uk) (Centre for Integrative Physiology, Biomedical Sciences University of Edinburgh)

Software for multichannel neuroscience data analysis.
- C++/MATLAB tools for analysis, visualization and causal modeling of in vitro multi-electrode array electrophysiology data 

- ENGINEERING GOAL: These software tools were developed to streamline the study cortical network activity primarily on planar     8x8 electrode arrays. 

- SCIENTICIC GOAL: Comparing group data between different models of Autism/Intellectual Disability, developmental time-series     analysis.


- Accepts sinals from Panasonic/MED64 MEA system (.MED), and Multichannel systems (.MCD - via Neuroshare API) recordings.
- Can incorporate simultaneously acquired single-cell signals (<=dual channel patch-clamp recordings, population calcium imaging   (2-Photon or epifluorescence) or single unit data) via Ephus .xsg or Molecular Devices (Axon) .abf files. Robust against        sampling jitter and sample-rate mismatches

- Can be configured to analyse network responses to episodic stimulation (e.g. patterned thalamocortial stimulation), or          autonomously detect spontaneous network activity.

Some highlights:

- Importing, Data archiving, inter-genotype statistical comparisons
- Data filtering: Local Field Potential (LFP) and Multi-Unit Analysis
- Inter-trial averaging or comparison
- Current-Source Density (CSD) analysis
- Visulaisation GUI (LFP, MUA, CSD)
- Movie rendering 
- Unit sorting
- Spike-triggered averaging
- MUA rate analysis
- Network activity propagation analysis (PCA, Granger Causaliy, Möbius transform, Information theory tools)
- Multichannel time (peri-event histograms, cross-correlation) and frequency domain (Wavelet spectrograms, Multi-taper Fourier    spectrograms, coherence) analysis

- Incorporates several published toolboxes - dependency details and acknowledgements to follow.
- Documentation to follow, data structures dependencies are detailed in ./instructions/data structures.docx
