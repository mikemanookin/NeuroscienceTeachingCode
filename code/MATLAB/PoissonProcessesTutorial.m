% PoissonProcessesTutorial.m
% by Mike Manookin 11/2019


% Before we get into Poisson distributions and the processes they describe, 
% we need to first introduce Gaussian distributions. The reason is that 
% Poisson distributions are just a special type of Gaussian distribution. 
% This is the equation for a Gaussian distribution:

% f(x) = 1/sqrt(2*pi*sigma^2) * exp(-(x - mu).^2/(2*sigma^2))
% mu is the mean of the distribution

% sigma is the standard deviation of the distribution

% It's alright if equations aren't your thing. The thing to remember is 
% that you have some distribution of your data and that there is a 
% mean/average (mu) and another parameter describing how spread out the 
% distribution is---the standard deviation (sigma). The following code 
% allows us to play with those two parameters and see how it affects the 
% Gaussian distribution. Go ahead and change 'mu' and 'sigma' and then run 
% the code to see how the distribution changes.

% A Poisson distribution is a special type of Gaussian distribution in which the
% mean and variance are equal. Many physical and biological phenomena are well-
% described by Poisson distributions.

% Let's make our own Gaussian function. We can define the function as a
% separate piece of code or define an 'anonymous function' in the current
% code. I'll do the latter to show you how it works. (It can come in really
% handy.)
gaussfun = @(p,x)(1/sqrt(2*pi*p(2)^2) * exp(-(x - p(1)).^2/(2*p(2)^2)));

% We can now pass x-axis values (x) and parameters for the Gaussian (p) to
% the anonymous function and get our output. 
% Note: 
%   p(1) = mean
%   p(2) = standard deviation

% First, let's try to unpack the Poisson distribution a bit by making a generic
% Gaussian distribution with a mean (mu) and standard deviation (sigma). Choose
% whatever numbers you like.
mu = 2.0;
sigma = 2.6;

% Define the x-axis.
x = -10:0.1:10;

% Pass the parameters to our anonymous function:
y = gaussfun([mu,sigma],x);

figure(100); clf;
hold on
plot(x,y);
hold off;
xlabel('Gaussian distribution'); 

%%
% Instead of using a function, we can also generate a bunch of random
% numbers that have a Gaussian distribution using MATLAB's built-in 'randn'
% function. If you want more information about any function type 'help
% functionName' or 'doc functionName' into the command window (e.g.,
% >> doc randn

% Generate 10^5 random numbers.
y2 = sigma * randn(1,1e5) + mu;

% Let's plot the histogram of the numbers we generated and then plot the
% function as well.

figure(100); clf;
hold on
histogram(y2, 50, 'Normalization', 'probability'); axis tight;
plot(x, y/max(y)*max(ylim), 'b', 'LineWidth', 2); % Need to scale to the same maximum value as the histogram
hold off;
xlabel('Gaussian distribution'); ylabel('probability');

%%
% Go ahead and play around with the parameters for the Gaussian distribution. 
% What you'll notice as you change the parameters, is that the values along 
% the x-axis change, but the general, bell-curve, shape stays the same.

% A Poisson distribution is just a special case in which the mean and 
% variance are equal. The standard deviation ( sigma ) is the square root 
% of the variance and, since the mean and variance are equal, it is also 
% the square root of the mean ( mu ).

% Define the mean/variance.
mu = 3.0; 

% The standard deviation (sigma) is just the square root of this value.
sigma = sqrt(mu);

% Let's do the same thing now for our Poisson distribution.
y = gaussfun([mu,sigma],x);
y2 = sigma * randn(1,1e5) + mu;

% Let's plot the histogram of the numbers we generated and then plot the
% function as well.

figure(100); clf;
hold on
histogram(y2, 50, 'Normalization', 'probability'); axis tight;
plot(x, y/max(y)*max(ylim), 'b', 'LineWidth', 2); % Need to scale to the same maximum value as the histogram
hold off;
xlabel('Gaussian distribution'); ylabel('probability');

%
% All that's happened is that we've constrained the Gaussian distribution 
% to have only one free parameter instead of two. Wasn't that easy?!

% So, why did I make you learn this then? Well, it turns out that there are 
% many physical processes that are well described with Poisson statistics. 
% Here are a few examples that are relevant to you and me:
%
%   1) In physics, the frequency of photon arrivals. 
%   2) Neural spike counts (not exactly due to the refractory period)...
%   3) Frequency of vesicle release at synapses

% Let's look at concrete example. Photon arrivals at a detector obey 
% Poisson statistics so let's simulate a distribution of photon counts 
% based on a mean values. Now, if we're talking about photon counts in, 
% say, a 100 ms window, we aren't going to observe any counts that are less 
% than zero, so we can clip the negative values out.

mu = 5.0;
sigma = sqrt(mu);

y2 = sigma * randn(1,1e5) + mu;

% Photon counts are discrete (integers), so let's round.
y2 = round(y2);

% We can't have negative photon counts, so set them to zero.
y2(y2 < 0) = 0;

% Let's get the underlying function too.
x = min(y2) : 0.1 : max(y2);
y = gaussfun([mu,sigma],x);

figure(100); clf;
hold on
histogram(y2, max(y2)+1, 'Normalization', 'probability'); axis tight;
plot(x(x>=0), y(x>=0)/max(y)*max(ylim), 'b', 'LineWidth', 2); % Need to scale to the same maximum value as the histogram
hold off;
xlabel('Gaussian distribution'); ylabel('probability');

% That's all there is to it. We can see the probability (y-axis) that we 
% observe x number of photons (x-axis). This type of visualization can be 
% very helpful as we try to understand neural systems, etc. There are many 
% examples of this, but for the time being I would like you to become 
% comfortable with this concept. So, please play around with the code to 
% familiarize yourself with how things change as we shift the parameters 
% around.

% Gaussian functions will keep coming up, not just in understanding 
% neuroscience concepts, but also in your data analysis. For example, many
% of the tools that simplify statistical calculations (such as t-tests) 
% assume that your data conforms to a Gaussian distribution. If you're not
% certain that this is the case, then you need to use other methods. 

%% Cumulative distribution functions
% So far in this tutorial, we have explored Gaussian/Poisson probability
% density functions (PDF). Now, we want to turn our attention to the 
% cumulative distribution functions (CDF) for these distributions.

% So what is a cumulative distribution function?
%
% To answer that question, let's back up and give an example of a Gaussian
% distribution. The classic example, of course, is the distribution of
% IQ's. The mean IQ is, by definition, 100 and the standard deviation is
% 15. Let's make a little IQ PDF/distribution.
mu = 100;
sigma = 15;

x = 0:200;
y = gaussfun([mu,sigma],x);
y = y/sum(y); % Area under the curve must equal 1.

figure(100); clf;
hold on
plot(x, y, 'b', 'LineWidth', 2); % Need to scale to the same maximum value as the histogram
hold off;
xlabel('distribution of IQ scores'); ylabel('probability');

%%
% Alright, let's now look at the probability of a person having an IQ of
% exactly 120.
p120 = y(x == 120)

% For our example, that turns out to be about 1%... It turns out that the
% 120 threshold is pretty important... In Malcolm Gladwell's book
% 'Outliers' he cites studies indicating that up to 120, IQ is a reliable
% predictor of a person's success, but above that threshold, it no longer
% correlates with success. The argument is that the IQ test captures
% something about intelligence, but there are aspects of intellect that it
% does not assess adequately. I digress, but it's a fun read. 

% OK, so what if we wanted to know what percentage of the population have
% an IQ less than 120? How would we calculate that? Well, knowing the value
% at 120, isn't particularly informative... This is where the cumulative
% distribution function becomes useful.

% To find the probability of someone having an IQ less than 120, we need to
% take the definite integral from [0,119]. This calculates the area under
% the distribution (i.e., cumulative distribution function) from  0 to 119. 
cumPr120 = sum(y(x <= 120))

% That gives us a value of 0.90, meaning that 90% of the population has an
% IQ less than 120. We could do the same thing and get the probability that
% a random person on the street has an IQ less than or equal to any value.
% One way to calculate this would be in a 'for' loop:

cumPr = zeros(size(x)); % Good practice to pre-allocate memory...
for k = 1 : length(x) % This is looping through each value in x...
    cumPr(k) = sum(y(x <= x(k)));
end

% Let's look at the result.
figure(100); clf;
hold on
plot(x, cumPr, 'b', 'LineWidth', 2); 

% Let's highlight the value at 120 for fun.
plot(120, cumPr(x == 120), 'ko');
plot(120*ones(1,2), [0 cumPr(x == 120)],'k--');
plot([0 120], cumPr(x == 120)*ones(1,2),'k--');

hold off;
xlabel('IQ score'); ylabel('cumulative probability');

% I've highlighted the value at 120 to illustrate what we discussed above.
% This is the cumulative distribution function (CDF) for the Gaussian
% distribution that we looked at earlier. That's all there is to it!!

%%
% It turns out that there are several ways to calculate the CDF. The for
% loop example above is one way, and I would like to show you two other
% ways of doing the calculation. Both ways are easier than the for loop
% example:

% 1) 'cumsum': Matlab has a built-in function called 'cumsum' that does
% exactly what we did in the for loop, but using this function just means
% fewer lines of code that we have to write. Let's do it!
cumPrCumSum = cumsum(y);

% Let's plot the output from our for loop above and the new way to make
% sure they align.

figure(100); clf;
hold on
plot(x, cumPr, 'b', 'LineWidth', 2); 
plot(x, cumPrCumSum, 'r', 'LineWidth', 2); 
hold off;
xlabel('IQ score'); ylabel('cumulative probability');

% Yep. We only see one curve, because the two curves are identical. Here's
% the second way.

% 2) 'normcdf': This function generates a Gaussian (also called 'Normal')
% cumulative distribution function based on mean and standard deviation.
% Check it out.
cumPrFun = normcdf(x, mu, sigma);

% Let's check this as well.
figure(100); clf;
hold on
plot(x, cumPr, 'b', 'LineWidth', 2); 
plot(x, cumPrFun, 'r', 'LineWidth', 2); 
hold off;
xlabel('IQ score'); ylabel('cumulative probability');

% No big deal, right?! You're just adding up probabilities. Note: if you're
% dealing with a Poisson distribution, the function is 'poisscdf'.

%% Fitting a cumulative distribution function

% One thing you'll be doing a lot of in your research is trying to model
% your data. Imagine that you have some data that you want to present. The
% first thing to do is show the data. Let's make some pretend data of some
% IQ scores:
dataHist = sigma * randn(1,1e3) + mu;
[data, binEdges] = histcounts(dataHist, 100,'Normalization', 'cdf');
x = (binEdges(1:end-1) + binEdges(2:end))/2; % Get the center of the bins.

% Alright, we've created a cumulative distribution of IQ scores. Now, let's
% fit the data. I'll show you two different functions that I typically use
% for nonlinear fitting...

% Another anonymous function.
fitfun = @(p,x)(normcdf(x,p(1),p(2)));

% 1) 'nlinfit'
params = nlinfit(x,data,fitfun,[100 15]);

% Let's take a look at our fitted parameters:
params

% Hopefully, they're pretty close to the mu and sigma values we used to
% generate our pretend data.

% 2) 'lsqcurvefit'
params2 = lsqcurvefit(fitfun,[100 15],x,data);
params2

% Let's see how good our fits were.
figure(100); clf;
hold on
plot(x, data, 'b.'); 
plot(x, fitfun(params,x), 'r', 'LineWidth', 1); 
hold off;
xlabel('IQ score'); ylabel('cumulative probability');

% There! We fit the data. That's all there is to it. The two fitting
% function each have their own advantages and I'm happy to discuss those,
% but for now, I just wanted you to see the process in action.

