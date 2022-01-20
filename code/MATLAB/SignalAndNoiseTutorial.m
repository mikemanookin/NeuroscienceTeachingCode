% SignalAndNoiseTutorial.m
% by Mike Manookin 11/2019

x = 0 : 10;

p1 = poisspdf(x, 1);


p2 = poisspdf(x, 2);

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



% Let's generate data from Equation 10 from Baylor et al 1979:
m = 0.53;
a = 1.2;
sigma0 = 0.15;
sigma1 = 0.3;
x = -1:0.1:3;

p = zeros(1, length(x));
for j = 1 : length(x)
    r = x(j);
    for k = 0:3
        p(j) = p(j) +  (exp(-m)*m^k)/factorial(k) .* 1./(2*pi*(sigma0^2+k*sigma1^2)).^0.5 .* exp(-(r - k*a).^2 / (2*(sigma0^2 + k*sigma1^2)));
    end
end

% Divide by the number of k's you've iterated over..
p = p / 4;

figure(1); clf;
plot(x,p)

