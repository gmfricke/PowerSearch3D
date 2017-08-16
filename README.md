# PowerSearch3D

This program reads in observations of biological swarm movement stored in CSV files. The expected contents of the files are positions and time for each moving entity. The program can create movement patterns that start at the same locations as the observed entities and bounded by their movement extrema. The generated movement patterns can be powerlaw (Levy-like, Brownian, CRW, etc). Using these model movement patterns the program provides data on the number of idealised targets encountered by the idealised movement. The point is to be able to compare the efficiency of searcher-target contact given observed movement data and compare that to idealised movement models.

The program requires QT4 and Ubuntu 14.04. 

You will need to install the following packages:

'''
sudo apt install libboost-dev
sudo apt install freeglut3-dev
'''
