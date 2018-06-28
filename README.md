# PIR_LiTianyi
the code, note and report of PIR 'Deep Learning pour la création de modèles réduits'

The 'Code' file contains the MATLAB code i used in the project. 

The 'Report' file contains the new result report that i made to my prof. 

The 'final report' file contains the final report both in tex and pdf. 

The 'Bibliography_research' file contains the PPT and report about my biblio research. 

The 'PPT' file contains the PPT of presentation, sorry i haven't change the diagram. 

In the 'code' file, 'Read_dat.m' are used to read the '.dat' file from ABAQUS. But when i change the load combinations in ABAQUS,i am not succeed in writing the max stress and their position in to '.dat'. So i just write the result directly into my matlab code. But using 'Read_dat.m' can be much easier if there are results in '.dat' file.

The code are written in the basis of the code provided in http://gaussianprocess.org/gpml/
Please make sure you have download this code library before launching my code.

Advices for students who continue this projet:
There are 3 books which are very helpful in this projet: 

1.Engineering Design via Surrogate Modelling:A Practical Guide
This book is the most comprehensive book. Not only the Grassian Process are inculded, but also other different surrogate modeling method. It also discutes the sampling method, multi-output problem and infill criteria. (almost everything about the surrogate modeling) . But the code library is not as convenient as 'Gaussian Processes for Machine Learning' (in my opinion).

2.Gaussian Processes for Machine Learning
This book mainly provide the code library. It has different explination of Grausiaan Process from different views. But the sampling method as well as multi-output problem are not discussed. (but scaling problems are discussed !)

3.Thèse of Ankit CHIPLUNKAR
The code in this book are written in the basis of code library provided in 'Engineering Design via Surrogate Modelling:A Practical Guide'. This book gives a very brief explanation of GP prinple. So i advise you to learn GP prinple from this book.(the language is much easier to understand !)  But the scaling and  multi-output problem are discussed in this book.
