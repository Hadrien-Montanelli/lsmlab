%% Add Gypsilab functions to PATH:

addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/miscellaneous')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openBmm')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openDom')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openEbd')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openFem')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openFfm')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openHmx')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openMmg')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openMsh')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openOpr')
addpath('/Users/montanelli/Dropbox/HM/WORK/MATLAB/gypsilab/openRay')

%% Preferences:

clc, close all
clf
format short 
format compact

%% Plotting:

lw = 2; fs = 30; ms = 18;
set(0, 'defaultlinelinewidth', lw)
set(0, 'defaultaxesfontsize', fs)
set(0, 'defaultlinemarkersize', ms)
LW = 'linewidth'; 
FS = 'fontsize'; 
MS = 'MarkerSize'; 
HV = 'HandleVisibility';
clr = parula(300);
parulaB = [0 0.447000000000000 0.741000000000000];
parulaY = [0.929000000000000 0.694000000000000 0.125000000000000];
parulaR = [0.850000000000000 0.325000000000000 0.098000000000000];
parulaM = [0.494000000000000 0.184000000000000 0.556000000000000];
parulaG = [0.466000000000000 0.674000000000000 0.188000000000000];