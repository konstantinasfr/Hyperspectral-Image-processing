function flg=iswithin(x,lo,hi)
% returns T for values within range of input
% SYNTAX:
%  [log] = iswithin(x,lo,hi)
%      returns T for x between lo and hi values, inclusive
flg= (x>=lo) & (x<=hi);