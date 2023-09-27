
function y = G_model(params, x, n_waves)
    y = zeros(size(x));

    for i = 1:n_waves
        A = params((i - 1) * 3 + 1);
        B = params((i - 1) * 3 + 2);
        C = params((i - 1) * 3 + 3);
        y = y + A * sin(B * x + C);
    end
end