function [Y] = crop_image(X,n)

m = size(X);

% find center
x_center = floor(m/2);
hw       = floor(n/2);
x_left   = x_center - hw;
x_right  = x_left + n - 1;
 
Y = X(x_left(1):x_right(1),x_left(2):x_right(2),x_left(3):x_right(3));


end

