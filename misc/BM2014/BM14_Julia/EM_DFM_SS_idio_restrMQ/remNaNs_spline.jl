

function remNaNs_spline(X, options)

    X = copy(X);
    T, N = size(X);
    k = options.k;
    indNaN = isnan.(X);

    if options.method == 1 # replace all the missing values
            for i = 1:N  
                x = X[:, i];
                nas = findall(indNaN[:,i]);
                x[nas].= nanmedian(x);
                x_MA = filt(ones(2*k+1)/(2*k+1), 1, [x[1]*ones(k);x;x[end]*ones(k)]);
                x_MA = x_MA[2*k+1:end];
                x[nas] = x_MA[nas];
                X[:,i] = x;
            end
    elseif options.method == 2 # replace missing values after removing leading and closing zeros
            rem1 = dropdims(sum(indNaN, dims = 2), dims = 2) .> N*0.8;
            nanLead = cumsum(rem1) .== 1:T;
            nanEnd = cumsum(rem1[end:-1:1]) .== 1:T;
            nanEnd = nanEnd[end:-1:1];
            nanLE = nanLead .| nanEnd;
            X = X[.!nanLE,:] # X[nanLE,:] = [];
            indNaN = isnan.(X);
            for i = 1:N  
                x = X[:,i];
                isnanx = isnan.(x);
                wcc = findall(.~isnanx);
                t1 = minimum(wcc);
                t2 = maximum(wcc);
                spline_interpolator = CubicSpline(x[wcc], wcc);
                x[t1:t2] = spline_interpolator.(t1:t2);
                isnanx = isnan.(x);
                x[isnanx] .= nanmedian(x);
                x_MA = filt(ones(2*k+1)/(2*k+1), 1, [x[1]*ones(k);x;x[end]*ones(k)]);
                x_MA = x_MA[2*k+1:end];
                x[isnanx] = x_MA[isnanx];
                X[:,i] = x;
            end
    elseif options.method == 3 #only remove rows with leading and closing zeros
            rem1 = dropdims(sum(indNaN, dims = 2), dims = 2) .== N;
            nanLead = cumsum(rem1) .== 1:T;
            nanEnd = cumsum(rem1[end:-1:1]) .== 1:T;
            nanEnd = nanEnd[end:-1:1];
            nanLE = nanLead .| nanEnd;
            X = X[.!nanLE,:] # X[nanLE,:] = [];
            indNaN = isnan.(X);
    elseif options.method == 4 # remove rows with leading and closing zeros & replace missing values
        rem1 = dropdims(sum(indNaN, dims = 2), dims = 2) .== N;
        nanLead = cumsum(rem1) .== 1:T;
        nanEnd = cumsum(rem1[end:-1:1]) .== 1:T;
        nanEnd = nanEnd[end:-1:1];
        nanLE = nanLead .| nanEnd;
        X = X[.!nanLE,:] # X[nanLE,:] = [];
        indNaN = isnan.(X);
            for i = 1:N  
                x = X[:,i];
                isnanx = isnan.(x);
                wcc = findall(.~isnanx);
                t1 = minimum(wcc);
                t2 = maximum(wcc);
                spline_interpolator = CubicSpline(x[wcc], wcc);
                x[t1:t2] = spline_interpolator.(t1:t2);
                isnanx = isnan.(x);
                x[isnanx] .= nanmedian(x);
                x_MA = filt(ones(2*k+1)/(2*k+1), 1, [x[1]*ones(k);x;x[end]*ones(k)]);
                x_MA = x_MA[2*k+1:end];
                x[isnanx] = x_MA[isnanx];
                X[:,i] = x;
            end
    elseif options.method == 5 #replace missing values  
            for i = 1:N  
                x = X[:,i];
                isnanx = isnan.(x);
                wcc = findall(.~isnanx);
                t1 = minimum(wcc);
                t2 = maximum(wcc);
                spline_interpolator = CubicSpline(x[wcc], wcc);
                x[t1:t2] = spline_interpolator.(t1:t2);
                isnanx = isnan.(x);
                x[isnanx] .= nanmedian(x);
                x_MA = filt(ones(2*k+1)/(2*k+1), 1, [x[1]*ones(k);x;x[end]*ones(k)]);
                x_MA = x_MA[2*k+1:end];
                x[isnanx] = x_MA[isnanx];
                X[:,i] = x;
            end
    end
    return X, indNaN
end
    