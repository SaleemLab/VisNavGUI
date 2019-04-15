function [d,best_i] = dprime(V,class_idx,noise_idx)
% DPRIME
% Prepares indices of classified points (class_idx). If no noise_idx is
% specified, takes each class and treats all points outside of class as noise.
% Otherwise, treats class_idx as signal_idx.
% Runs DPRIMECALC on each of these permutations, returns the best d value.
%
% MOVED TO DPRIMECALC.M:
% Calculates d', comparing the selected cluster against all other datapoints
% USAGE: d = dprime(V,signal_idx) OR dprime(V,signal_idx,noise_idx)
%    V = clustered data points, in 3D
%    signal_idx = indices of points in V that are part of the signal
%     noise_idx = indices of points in V that are part of the noise
% If noise_idx is not provided, function will assume all non-signal indices
%    are noise indices
%
% SEE ALSO: DPRIMECALC, HOWMANYCLASSES
%
% change log:
%   created by Andrew Zaharia 2009-01-07
%   AZ 2009-01-22: split into dprime & dprimecalc
%   AZ 2009-07-08: changed so that all points not in class in question are
%      counted as noise, rather than all points in other classes

d = 0;

if nargin == 4
   signal_idx = class_idx{1};
   
   d = dprimecalc(V,signal_idx,noise_idx);
   
elseif nargin < 4
   if size(class_idx,2) > 1
      
      k = howManyClasses(class_idx);
      
      clusters = 1:k;
      for i = clusters
         signal_idx = class_idx{i};
         %all clusters, not including a
         noise_idx = setdiff(1:length(V),signal_idx)';
         
         d_temp = dprimecalc(V,signal_idx,noise_idx);
         if d_temp >= d
            d            = d_temp;
%             signal_idx = signal_idx;
%             noise_idx =  noise_idx;
            best_i = i;
         end
      end
      
   elseif size(class_idx,2) == 1
      signal_idx = class_idx;
      noise_idx = 1:length(V);
      noise_idx(signal_idx) = [];
      noise_idx = noise_idx';
      
      d = dprimecalc(V,signal_idx,noise_idx);
   end
end

