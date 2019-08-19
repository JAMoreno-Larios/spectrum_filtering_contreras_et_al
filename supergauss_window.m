function w = supergauss_window(axis, width, center, sigma, power)
    % Generates a super-Gaussian window centered at center.
    w = exp(-(((axis - center)/width).^2 / 2 / sigma ^ 2) .^ power);
    w(abs(axis - center) > width/2) = 0;
end