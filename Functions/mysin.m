
function y = mysin (n,f,p)
%DESCRIPTION:
% Generate sin time series.
%
%INPUT:
%- n: Number of datapoints
%- f: Number of periods in n points
%- p: Phase in units of pi
%
%OUTPUT:
%- y: vector of n points containing the sin wave

t = 0:n-1;
y = sin(2*pi*f*t/n+p*pi);
