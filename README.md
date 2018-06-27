# learning-ds-tutorial
This package includes demo scripts and a GUI simulation for learning stable non-linear Dynamical Systems (DS) from demonstrations using SEDS [1] and LPV-DS [2,3] approaches developed in LASA-EPFL.

## Installation Instructions
This package contains a set of submodules. In order to clone the package, just run the following run:
```
git clone git@github.com:epfl-lasa/RSS2018Tutorial.git
```
Or if the ssh access has not been set up on your PC 
```
git clone https://github.com/epfl-lasa/RSS2018Tutorial.git
```
After cloning one must intialize/download the submodules with the following commands:
```
cd ~./learning-ds-tutorial
git submodule init
git submodule update
```
In case you want to update the submodules to their latest version, you can do so with the following command:
```
git submodule update --remote
```

[![](https://github.com/nbfigueroa/learning-ds-tutorial/blob/master/img/GUI_2.png)](https://www.youtube.com/watch?v=5fQO9Oluih0)

*A guided video explaining how to use the GUI, can be found in this [link](https://www.youtube.com/watch?v=5fQO9Oluih0)*

**References**     
[1] Khansari Zadeh, S. M. and Billard, A. (2011) Learning Stable Non-Linear Dynamical Systems with Gaussian Mixture Models. IEEE Transaction on Robotics, vol. 27, num 5, p. 943-957.    
[2] Mirrazavi Salehian, S. S. (2018) Compliant control of Uni/ Multi- robotic arms with dynamical systems. PhD Thesis.  
[3] Figueroa, N and Billard, A. ".." [Under Review]

**Contact**: [Nadia Figueroa](http://lasa.epfl.ch/people/member.php?SCIPER=238387) (nadia.figueroafernandez AT epfl dot ch)
