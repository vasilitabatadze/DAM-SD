function [arg] = CC(x,n,alfa)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arg=gamma(n*alfa+1)/gamma(n*alfa+alfa)*x.^(n*alfa+alfa-1);
end

