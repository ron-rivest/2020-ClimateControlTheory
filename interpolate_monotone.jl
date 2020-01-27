# interpolate_monotone.jl
# Ronald L. Rivest
# January 26, 2019

# Use monotone cubic splines to interpolate a curve fitting given data.

"""
    interpolate_monotone_curve(x, y)

    Input: x -- list of n distinct x-coords, 
                e.g. [x_1, ..., x_n] = [2020, 2050, 2100]
                assumed to be in increasing order
           y -- list of n target values, 
                e.g. [y_1, ..., y_n] = [1.1, 1.7, 2.0]
           
    Output: curve -- function (piecewise cubic spline) mapping from x to y

    Discussion:
           The desired curve is monotonic cubic interpolation
               https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
           The "talk" section re correctness of step 5 is also relevant & necessary.

    Refs:
    https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
    http://bl.ocks.org/niclasmattsson/7bceb05fba6c71c78d507adae3d29417    

    Usage example:
	f = interpolate_monotone_curve([0, 1, 2, 3, 4], [0, 1, 4, 9, 16]);
	for x in 0:0.5:4.0
	    val = f(x)
            println(x, " squared is about ", val)
        end
"""

# Define basic cubic Hermite splines

h00(t) = (1+2t)*(1-t)^2
h10(t) = t*(1-t)^2
h01(t) = t^2 * (3-2t)
h11(t) = t^2 * (t-1)


function interpolate_monotone_curve(x, y)

    # deal with length issues
    @assert length(x)==length(y)
    n = length(x)
    if n == 0
        return function (x); return 0; end
    elseif n == 1
        return function (x); return y[1]; end
    end

    # 0. sort x and y similarly so x is in nondecreasing order
    p = sortperm(x)
    x = x[p]
    y = y[p]

    # 1. get consecutive differences and slopes
    dx = diff(x)
    dy = diff(y)
    delta = dy ./ dx      # note: assumes xs are distinct so no zerodiv!

    # 2. Compute m

    m = fill(0.0, n)      # slopes
    for k in 2:n-1
        if delta[k-1]*delta[k] <= 0.0 
            m[k] = 0
        else
            m[k] = (delta[k-1]+delta[k])/2.0
        end
    end
    m[1] = delta[1]
    m[n] = delta[n-1]


    # 3. Check for flat spline segment
    alpha = fill(0.0, n)
    beta = fill(0.0, n)
    for k in 1:n-1
        if delta[k] == 0
            m[k] = 0
            m[k+1] = 0
        else
            # 4. local max or mins in data
            alpha[k] = m[k]/delta[k]
            beta[k] = m[k+1]/delta[k]
            if alpha[k] < 0
                m[k] = 0
            end
            if beta[k] < 0
                m[k+1] = 0
            end
            # 5. Prevent overshoot and ensure monotonicity
            # This is buggy, but it is what the wikipedia article seems to say was OK
            # if alpha[k] > 3 || beta[k] > 3
            #     m[k] = 3 * delta[k]
            # end
            # Here is fix, from "talk" post by Berland
            if alpha[k]^2 + beta[k]^2 > 9
                tau = 3 * (alpha[k]^2 + beta[k]^2)^(-0.5)
                m[k] = tau * alpha[k] * delta[k]
                m[k+1] = tau * beta[k] * delta[k]
            end
        end
    end

    # return interpolant function -- what y is predicted for input u?
    return function (u)
        # Give correct result if u too small or too large
        if u>=x[end]
            return y[end]
        elseif u<=x[1]
            return y[1]
        end
        # Search for the interval [x_k, x_{k+1}] that u is in
        low = 1
        high = n
        while low <= high
            mid = Int64(floor((low+high)/2))
            if x[mid] < u;   low = mid + 1  end
            if x[mid] > u;   high = mid - 1 end
            if x[mid] == u;  return y[mid] end
        end
        k = max(1, high)

        # Interpolate
        t = (u-x[k]) / dx[k]
        return y[k]*h00(t) + dx[k]*m[k]*h10(t) +
               y[k+1]*h01(t) + dx[k]*m[k+1]*h11(t)
    end # of function f(u)
end

################################################################

function test_interpolate_monotone_curve()
    f = interpolate_monotone_curve([1,2,3,4], [1,4,9,9])
    for u in 0:0.5:5.0
        val = f(u)
        println("Square of ", u, " is about ", val)
    end
end

# test_interpolate_monotone_curve()

################################################################
## The following illustrates non-monotonicity bug if
## step 5 is not implemented correctly.
## The correct implementation is from "talk" section of Wikipedia.

function demo_bug()
    xs = [0, 0.30, 0.5]
    ys = [0, 0.05, 0.5]
    f = interpolate_monotone_curve(xs, ys)
    ### BUGGY -- THIS f SHOULD BE MONOTONE BUT ISN'T!
    for i in 0.0:0.01:0.5
        println(i, " ",f(i))
    end
end

# demo_bug()

