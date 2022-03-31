clear;
clc;

Yeast_MP_sucrose_function(1,0,250*1e6*4.6e-11,0,101);
Yeast_MP_sucrose_scratp_function(1,0,250*1e6*4.6e-11,0,102);

load data101
load data102

plot(T101,Y101(:,2),T102,Y102(:,2))