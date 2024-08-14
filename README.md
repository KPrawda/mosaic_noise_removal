# Mosaic - a method to remove non-stationary noise from remeated sweeo measurements
<br>
Repository for the manuscript _Non-stationary Noise Removal from Repeated Sweep Measurements_ by Karolina Prawda, Sebastian Schlecht and Vesa Välimäki. 

Submitted to JASA Express Letters on 30.04.2024

**Abstract**: <br>
Acoustic measurements using sine sweeps are prone to background noise and non-stationary disturbances. Repeated measurements can be averaged to improve the resulting signal-to-noise ratio. However, averaging leads to poor rejection of non-stationary high-energy disturbances and, in the case of a time-variant environment, causes attenuation at high frequencies. This paper proposes a robust method to combine repeated sweep measurements using across-measurement median filtering in the time~frequency domain. The method, called Mosaic, successfully rejects non-stationary noise, suppresses background noise, and is more robust towards time variation than averaging. The proposed method allows high-quality measurement of impulse responses in a noisy environment.  

<br><br>

![Mosaic non-stationary noise removal method](https://github.com/KPrawda/mosaic_noise_removal/blob/main/Mosaic.PNG)
The Mosaic-T method uses median filtering in time domain to remove non-stationary noise event that do not overlap in time. None of the sweeps needs to be completely devoid of non-stationary noise, but every time sample must be clean in at least half of the total repetitions. 

<br><br>

![Mosaic non-stationary noise removal method](https://github.com/KPrawda/mosaic_noise_removal/blob/main/Mosaic-TF.png)
The Mosaic-TF method uses median filtering in time-frequency domain to remove non-stationary noise event that may overlap either in time or frequency. None of the sweeps needs to be completely devoid of non-stationary noise, but every time/frequency tile must be clean in at least half of the total repetitions. 

<br><br>

Please see the folders for sound examples and figures.

**Cite as:**<br>
@article{prawda2024_mosaic,<br>
author = {Prawda, Karolina  and Schlecht, Sebastian J.  and Välimäki, Vesa},<br>
title = {Non-stationary Noise Removal from Repeated Sweep Measurements},<br>
journal = {J. Acoust. Soc. Am. Express Lett.},<br>
volume = {},<br>
number = {},<br>
pages = {},<br>
year = {2024},<br>
doi = {}}<br>
 
