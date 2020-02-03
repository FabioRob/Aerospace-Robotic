function y = zero_to_360(x)
% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
%
% This function reduces an angle to lie in the range
% 0 - 360 degrees.
%
% x - the original angle in degrees
% y - the angle reduced to the range 0 - 360 degrees
%
% User M-functions required: none
% ------------------------------------------------------------
if x >= 360
x = x - fix(x/360)*360;
elseif x < 0
x = x - (fix(x/360) - 1)*360;
end
y = x;
return