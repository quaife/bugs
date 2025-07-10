% File: bugsN.m
% Author: Bryan Quaife (Florida State University, Department of
% Scientific Computing)
% Contact: bquaife@fsu.edu
% Date: 2025-07-09

% Description:
% Run a Monte-Carlo simulation to compute probabilities that N bugs
% moving on the boundary of the unit circle will coalesce (all meet at a
% single point).

% Reference:
% Josh Briley and Bryan Quaife, "N Bugs on a Circle", Proceedings
% of the Royal Society A, 2025, under review.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% flag for plotting
iplot = ~true;
% time step size and number of time steps. 1e-2 works good for N<=100
dt = 1e-2;
% number of samples to estimate probability
nsamples = 1000;
% range of number of bugs to consider
N = (2:2:100);

% probability that bugs coalesce
probCoalesce = zeros(size(N));
% maximum number of steps required until it was clear the bugs would
% either cycle or coalesce for each value of N
maxsteps = zeros(size(N));
% counter for number of times it cycles for each value of N
ncycle = zeros(size(N));
% counter for number of times it coalesces for each value of N
ncoalesce = zeros(size(N));

% space for saving the total number of steps required to detect
% coalescing or cycling bugs for each value of N and each sample
totalsteps = zeros(numel(N),nsamples);

% loop over number of bugs
for k = 1:numel(N);
  % number of bugs
  nbugs = N(k);

  % build a matrix so that we can quickly look up index of all bugs
  % except itself
  ind_others = zeros(nbugs,nbugs-1);
  for ell = 1:nbugs
    ind_others(ell,:) = [(1:ell-1) (ell+1:nbugs)];
  end

  % set of random initial locations with the first bug always at (1,0)
  thetaAll = 2*pi*rand(nsamples,nbugs);
  thetaAll(:,1) = 0;

  % loop over the number of samples
  for i = 1:nsamples;
    % initialize that bugs have not converged
    converged = false;
    % grab initial condition
    theta = thetaAll(i,:);

    % plot the initial condition if iplot == true
    if iplot
      figure(1); clf; hold on
      plot(exp(1i*linspace(0,2*pi,1000)),'k','linewidth',3)
      plot(exp(1i*theta),'r.','markersize',40)
      axis equal
      axis(1.1*[-1 1 -1 1])
      pause
    end

    % keep track of number of steps to reach one of the steady states
    nsteps = 0;
    % while it is unclear which state the bugs will reach
    while ~converged
      nsteps = nsteps + 1;
      maxsteps(k) = max(maxsteps(k),nsteps);
      % difference between angles, but mod 2*pi
      omega = mod(diff([theta theta(1)]),2*pi);
      % compute velocity which is either +1 or -1, or zero if two bugs
      % are at the exact same location. All depends if the difference
      % between two bugs' angles, mod 2*pi, are greater than or less
      % than pi
      vel = (omega < pi) - (omega > pi);
      % if all velocities are the same, bugs are in a cycle
      if (norm(vel - 1) == 0 || norm(vel + 1) == 0)
        % increase ncycle counter
        ncycle(k) = ncycle(k) + 1;
        % set converged flag to true so that we don't do any more time
        % steps
        converged = true;
        % record total number of steps taken
        totalsteps(k,i) = nsteps;
        % break out of the while loop so that we don't next check if the
        % bugs are cycle
        break
      end
      for ell = 1:nbugs
        % if maximum distance between one of the bugs and all others is
        % less than pi, then it will inevitably coalesce
        if (max(mod(theta(ind_others(ell,:)) - theta(ell),2*pi)) < pi)
          % increase ncoalesce counter
          ncoalesce(k) = ncoalesce(k) + 1;
          % set converged flag to true so that we don't do any more time
          % steps
          converged = true;
          % record total number of steps taken
          totalsteps(k,i) = nsteps;
          % break out so that we don't look at the other bugs and
          % compare with all but itself
          break
        end
      end

      % if it is still unclear if bugs will cycler or coalesce, update
      % their locations
      theta = mod(theta + dt*vel,2*pi);

      % plot the current time step if desired
      if iplot
        figure(1); clf; hold on
        plot(exp(1i*linspace(0,2*pi,1000)),'k','linewidth',3)
        plot(exp(1i*theta),'r.','markersize',40)
        axis equal
        axis(1.1*[-1 1 -1 1])
        pause();
      end
    end
  end

  % probability of coalescing
  probCoalesce(k) = ncoalesce(k)/nsamples;
  % plot the probability of coalescing up to the current iterate for the
  % number of bugs
  clf;
  plot(N(1:k),probCoalesce(1:k),'b-','linewidth',4)
  xlim([N(1) N(end)])
  ylim([0 1])
  set(gca,'fontsize',20)
  xlabel('Number of Bugs','fontsize',24)
  ylabel('Probability of Coalescing','fontsize',24)
  pause(1e-2);
end

% fit the probability of coalescing to a linear curve in loglog space
p = polyfit(log10(N),log10(probCoalesce),1);
hold on
% plot the algebraic representation for the probability of coalescing
plot(N,10.^(polyval(p,log10(N))),'k--','linewidth',2)

h = legend('Monte Carlo','Algebraic Fit');
set(h,'box','off')
set(h,'location','northeast')
set(h,'fontsize',24)

