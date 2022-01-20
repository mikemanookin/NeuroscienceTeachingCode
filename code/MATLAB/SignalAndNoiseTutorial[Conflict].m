% SignalAndNoiseTutorial.m
% by Mike Manookin 11/2019

x = 0:2;
y = [0.81 0.17 0.02 0.64 0.28 0.08];

% Let's define anonymous functions.

% This would be for the single-photon condition.
pfun1 = @(m,k)([(m.^k * exp(-m))./factorial(k),((2*m).^k * exp(-2*m))./factorial(k)]);

% This would be for the coincident model:
pfun2 = @(m,k)([(m.^k * exp(-m))./factorial(k)+(m.^(k-1) * exp(-m))./factorial(k-1),((2*m).^k * exp(-2*m))./factorial(k)+((2*m).^(k-1) * exp(-2*m))./factorial(k-1)]);

m = 0.21;

% Single-photon model.
p1 = 0;


