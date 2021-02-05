# creaCoroDeteccion


Semi-automatic Centerline Extraction Framework

Efficiently obtaining a reliable coronary artery centerline from computed tomography angiography data is relevant in clinical practice. 

This project aims to implement centerline extraction semi-automatic methods. It is based on two seed points (start and end) placed by the radiologist. Then, It extracts interactively the centerline using a medialness image [1][2] and a minimal-cost path algorithm [3]. All code has been developed entirely in C++, using ITK and VTK libraries and the medical imaging framework CreaTools [4]. Currently, further medialness techniques are being implemented in order to do a comparison between them. 


# References:

[1] Esteban Correa-Agudelo, Leonardo Flórez-Valencia, Maciej Orkisz, Claire Mouton, Eduardo E. Dávila Serrano, Marcela Hernández Hoyos, "A Modular Workflow Architecture for Coronary Centerline Extraction in Computed Tomography Angiography Data", ICCVG 2014, Warsaw, Poland. (doi:10.1007/978-3-319-11331-9_19)

[2] CreaTools http://www.creatis.insa-lyon.fr/site/en/CreaTools_home

Further information:

http://www.creatis.insa-lyon.fr/

Master Thesis: Chapter 3
