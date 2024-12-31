function dy_dx = derivative(f, x, dx)
    y2 = f(x + dx);
    y1 = f(x - dx);
    dy_dx = (y2 - y1) ./ (2 * dx);
end