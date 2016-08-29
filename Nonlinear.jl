using BandedMatrices
using ForwardDiff
# y'' = F(x, y, y''),    x ∈ (a, b)
# a₀y(a) + b₀y'(a) = c₀
# a₁y(b) + b₁y'(b) = c₁
function BVP_nonlinear(order::Int, f::Function,
                       a::Number,  b::Number,
                       a₀::Vector, a₁::Vector;
                       N=Int((b-a)*15),   y=0,
                       tol = 1e-3,  MaxIt=Int(1e3))
    i = 2:N
    h = (b-a) / N
    x = [a:h:b;]
    if y==0
        slope = (a₁[end]/a₁[1] - a₀[end]/a₀[1]) / (b-a)
        y = convert(Array{Float64}, slope*x+a₀[end]/a₀[1])
    end

    ty = typeof(a₁[end])
    A = zeros(ty, N+1, 5)
    F = zeros(ty, N+1)

    ĥ = 1 / h^order

    A[1, 1]     = a₀[1]
    A[end, 5]   = a₁[1]
    F[1]   = a₀[end]
    F[end] = a₁[end]

    for d in 1:5
        for c in 1:order-1
            A[1,   d] += FPSF[c, d]   * a₀[c+1]  / h^c
            A[end, d] -= FPSF[c, 6-d] * a₁[c+1]  / h^c
        end
    end

    for b in i
        for j in 1:5
            A[b, j] = FPSC[order, j] * ĥ
        end
    end

    M = bzeros(N+1, N+1, 4, 4)
    J = bzeros(N+1, N+1, 4, 4)

    M[1, 1:4] = A[1, 1:4][:]
    M[2, 2:6] = A[2, 1:5][:]
    M[end-1, end-5:end-1]  = A[end-1, 1:5][:]
    M[end, end-3:end]      = A[end, 2:5][:]

    for i in 3:N-1
        M[i, i-2:i+2] = A[i, 1:5][:]
    end

    y_n1 = similar(y)
    updateF!(F, f, x, y, N, h, order)
    JG = 1.
    k = 0

    while norm(JG) >= tol
        compJ!(J, f, x, y, h, N, M, order)
        G = M * y - F
        JG = J \ -G
        y_n1 = JG + y
        updateF!(F, f, x, y_n1, N, h, order)
        y = y_n1
            k+=1
        if k > MaxIt
            warn("Max iterations reached")
            break
        end
    end
    return M, J, y, F, linspace(a, b, N+1)
end

function updateF!(F::Vector, f::Function, x::Vector,
                  y::Vector, N::Int, h::Number, order::Int)
    if order == 1
        for i = 2:N
            F[i] = f(x[i], y[i])
        end
    elseif order == 2
        F[2] = f(x[2], y[2], (y[3]-y[1])/2h)
        F[N] = f(x[N], y[N], (y[N+1]-y[N-1])/2h)
        for i = 3:N-1
            F[i] = f(x[i], y[i],
            (FPSC[1, 1]*y[i-2] + FPSC[1, 2]*y[i-1] + FPSC[1, 3]*y[i] + FPSC[1, 4]*y[i+1] + FPSC[1, 5]*y[i+2])/h)
        end
    elseif order == 3
        F[2] = f(x[2], y[2], (y[2]-y[1])/2h, (y[1]-2y[2]^2+y[3])/h^2)
        F[N] = f(x[N], y[N], (y[N+1]-y[N-1])/2h, (y[N-1]-2y[N]^2+y[N+1])/h^2)
        for i = 3:N-1
            F[i] = f(x[i], y[i],
            (FPSC[1, 1]*y[i-2] + FPSC[1, 2]*y[i-1] + FPSC[1, 3]*y[i] + FPSC[1, 4]*y[i+1] + FPSC[1, 5]*y[i+2])/h,
            (FPSC[2, 1]*y[i-2] + FPSC[2, 2]*y[i-1] + FPSC[2, 3]*y[i] + FPSC[2, 4]*y[i+1] + FPSC[2, 5]*y[i+2])/h^2)
        end
    elseif order == 4
        F[2] = f(x[2], y[2], (y[2]-y[1])/2h, (y[1]-2y[2]^2+y[3])/h^2,
            (FPSF[3, 1]*y[2] + FPSF[3, 2]*y[3] + FPSF[3, 3]*y[4] + FPSF[3, 4]*y[5] + FPSF[3, 4]*y[5])/h^3)
        F[N] = f(x[N], y[N], (y[N+1]-y[N-1])/2h, (y[N-1]-2y[N]^2+y[N+1])/h^2,
            (-FPSF[3, 5]*y[N] - FPSF[3, 4]*y[N-1] - FPSF[3, 3]*y[N-2] - FPSF[3, 2]*y[N-3] - FPSF[3, 1]*y[N-4])/h^3)
        for i = 3:N-1
            F[i] = f(x[i], y[i],
            (FPSC[1, 1]*y[i-2] + FPSC[1, 2]*y[i-1] + FPSC[1, 3]*y[i] + FPSC[1, 4]*y[i+1] + FPSC[1, 5]*y[i+2])/h,
            (FPSC[2, 1]*y[i-2] + FPSC[2, 2]*y[i-1] + FPSC[2, 3]*y[i] + FPSC[2, 4]*y[i+1] + FPSC[2, 5]*y[i+2])/h^2,
            (FPSC[3, 1]*y[i-2] + FPSC[3, 2]*y[i-1] + FPSC[3, 3]*y[i] + FPSC[3, 4]*y[i+1] + FPSC[4, 5]*y[i+2])/h^3)
        end
    else
        error("This solver dose not support 5th-order or higher order ODEs yet.")
    end
    return F
end

function compJ!(J::BandedMatrices.BandedMatrix, f::Function,
                x::Vector, y::Vector, h::Number, N::Int,
                M::BandedMatrices.BandedMatrix, order::Int)

    if order == 1

            fc = (xys) ->  f(xys[1], xys[2])

    elseif order == 2

            fc = (xys) ->  f(xys[1], xys[4],
                dot(FPSC[1, :], xys[2:6])/h)

            # ff = (xys) -> f(xys[1], xys[2], xys[3]-xys[1]/2h)

            # fb = (xys) -> f(xys[1], xys[2], xys[3]-xys[1]/2h)

            ff = (xys) ->  f(xys[1], xys[2],
                dot(FPSF[1, :], xys[2:6])/h)

            fb = (xys) ->  f(xys[1], xys[6],
                dot(FPSF[1, :], xys[2:6])/h)

    elseif order == 3

            fc = (xys) ->  f(xys[1], xys[4],
                dot(FPSC[1, :], xys[2:6])/h,
                dot(FPSC[2, :], xys[2:6])/h^2)

            ff = (xys) ->  f(xys[1], xys[2],
                dot(FPSF[1, :], xys[2:6])/h,
                dot(FPSF[2, :], xys[2:6])/h^2)


            fb = (xys) ->  f(xys[1], xys[6],
                dot(-FPSF[1, 1:5], xys[6:2])/h,
                dot(-FPSF[2, 1:5], xys[6:2])/h^2)

    elseif order == 4

            fc = (xys) ->
            f(xys[1], xys[4],
                dot(FPSC[1, :], xys[2:6])/h,
                dot(FPSC[2, :], xys[2:6])/h^2,
                dot(FPSC[3, :], xys[2:6])/h^3)

            ff = (xys) ->  f(xys[1], xys[2],
                dot(FPSF[1, :], xys[2:6])/h,
                dot(FPSF[2, :], xys[2:6])/h^2,
                dot(FPSF[3, :], xys[2:6])/h^3)


            fb = (xys) ->  f(xys[1], xys[6],
                dot(-FPSF[1, 5:1], xys[6:2])/h,
                dot(-FPSF[2, 5:1], xys[6:2])/h^2,
                dot(-FPSF[3, 5:1], xys[6:2])/h^3)

    else
        error("This solver does not support 5th-order or higher order ODEs yet.")
    end

    copyM!(J, M, N)

    if order == 1
        gg = t -> ForwardDiff.gradient(fc, t, Chunk{2}())
        # J[2, 2] -= gg([x[1], y[2]])[2]
        # J[N, N] -= gg([x[N], y[N]])[2]
        for i = 3:N-1
            g = gg([x[i], y[i]])[2]
            J[i, i] -= g
        end
    else
        gg = t -> ForwardDiff.gradient(fc, t, Chunk{6}())
        gf = t -> ForwardDiff.gradient(ff, t, Chunk{6}())
        gb = t -> ForwardDiff.gradient(fb, t, Chunk{6}())
        J[2, 2:6]   -= gf([x[2], y[2], y[3], y[4], y[5], y[6]])[2:6]
        J[N, N-4:N] += reverse(gb([x[N], y[N-4], y[N-3], y[N-2], y[N-1], y[N]]))[1:5]

        for i = 3:N-1
            J[i, i-2:i+2] -= gg([x[i],y[i-2],y[i-1],y[i],y[i+1],y[i+2]])[2:6]
        end
    end
    J
end

function copyM!(J, M, N)
    for i = 1:2
        J[i, 1:i+4] = M[i, 1:i+4]
    end

    for i = N:N+1
        J[i, i-4:i] = M[i, i-4:i]
    end

    for i = 3:N-1
        J[i, i-2:i+2] = M[i, i-2:i+2]
    end
    J
end

# Five-Point Stencil Central
const FPSC =
[
[1/12   -2/3    0.  2/3 -1/12]
[-1/12   4/3  -5/2  4/3 -1/12]
[-1/2     1.    0.   -1.  1/2]
[1.      -4.    6.   -4.    1]
]
# Five-Point Stencil F&B
const FPSF =
[
[-25/12  4.      -3.      4/3    -1/4]
[35/12  -26/3   19/2    -14/3   11/12]
[-5/2    9.      -12.      7.    -3/2]
[1.     -4.       6.      -4.      1.]
]