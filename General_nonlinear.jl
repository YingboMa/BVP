using ForwardDiff
# y'' = F(x, y, y''),    x ∈ (a, b)
# a₀y(a) + b₀y'(a) = c₀
# a₁y(b) + b₁y'(b) = c₁
function BVP_nonlinear(f::Function, a::Number,  b::Number,
                       a₀::Number,  b₀::Number, c₀::Number,
                       a₁::Number,  b₁::Number, c₁::Number;
                       N=Int((b-a)*20),
                       tol = 1e-3)
    y=randn(N+1)
    # N = 50
    # a = 0.0
    # b = 10.0
    i = 1:N-1
    x = linspace(a, b, N)
    h = (b-a) / N

    # a₀ = 1.
    # b₀ = 0.
    # c₀ = 0.

    # a₁ = 1.
    # b₁ = 0.
    # c₁ = 1.

    A = Array(Float64, N+1)
    B = Array(Float64, N+1)
    C = Array(Float64, N+1)
    F = Array(Float64, N+1)

    h² = 1 / h^2
    h₂ = h*2

    A[1] = a₀ - 3*b₀/h₂
    C[1] = 2*b₀ / h
    B[1] = -b₀ / h₂
    F[1] = c₀

    A[N+1] = b₁ / h₂
    C[N+1] = -2*b₁ / h
    B[N+1] = a₁ + 3*b₁ / h₂
    F[N+1] = c₁

    for j in i
        A[i+1] = B[i+1] = h²
        C[i+1] = -2*h²
    end

    M = Tridiagonal([A[2:N], C[N+1];],
                    [A[1], C[2:N], B[N+1];],
                    [C[1], B[2:N];]
                    )
    y_n = y
    updateF!(F, f, x, y_n, N, h)

    for i in x
        J = compJ(f, x, y_n)
        G = M*y_n - F
        @show size(J\-G)
        @show size(y_n*transpose(y_n))
        #@show size(G)
        #@show size(J * y_n * transpose(y_n))
        y_n1 = (J * y_n) \ (J * y_n * transpose(y_n) - G)
        updateF!(F, f, x, y_n1)
        y_n = y_n1
    end
    return yn_1
end

function myBessel(x, y, yp)
    ( -x*yp - (x^2-α) ) / x^2
end

function updateF!(F::Vector, f::Function, x::LinSpace, y::Vector, N::Int, h::Number)
    for i = 2:N - 1
        F[i] = f(x[i], y[i], (y[i+1] - y[i-1])/(2h))
    end
    return F
end

function compJ(f::Function, x, y::Vector)
    N = length(y)
    J = Tridiagonal(zeros(N - 1), zeros(N), zeros(N - 1))
    ff(xys) = f(xys[1],xys[2],(xys[3]-xys[4])/(2h))
    for i = 2:N - 1
        g = ForwardDiff.gradient(ff)([x[i],y[i],y[i-1],y[i+1]])
        J.d[i] = g[2]
        J.dl[i] = g[3]
        J.du[i] = g[4]
    end
    J
end