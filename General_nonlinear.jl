using ForwardDiff
# y'' = F(x, y, y''),    x ∈ (a, b)
# a₀y(a) + b₀y'(a) = c₀
# a₁y(b) + b₁y'(b) = c₁
function BVP_nonlinear(f::Function, a::Number,  b::Number,
                       a₀::Number,  b₀::Number, c₀::Number,
                       a₁::Number,  b₁::Number, c₁::Number;
                       N=Int((b-a)*15),         y=0,
                       tol = 1e-3,  maxIt=Int(1e3))
    i = 1:N-1
    h = (b-a) / N
    x = [a:h:b;]
    if y==0
        slope = (c₁/a₁ - c₀/a₀) / (b-a)
        y = convert(Array{Float64}, slope*x+c₀/a₀)
    end

    A = Array(Float64, N+1)
    B = Array(Float64, N+1)
    C = Array(Float64, N+1)
    F = zeros(Float64, N+1)

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
        A[i+1] = h²
        B[i+1] = h²
        C[i+1] = -2*h²
    end

    M = Tridiagonal([A[2:N], C[N+1];],
                    [A[1], C[2:N], B[N+1];],
                    [C[1], B[2:N];]
                    )
    y_n = y
    y_n1 = similar(y_n)
    updateF!(F, f, x, y_n, N, h)
    JG = 1.
    k = 0

    while norm(JG) >= tol
        J = compJ(f, x, y_n, h, N+1, M)
        G = M*y_n - F
        JG = J \ -G
        y_n1 = JG + y_n
        @assert all(isfinite(y_n1))
        updateF!(F, f, x, y_n1, N, h)
        y_n = y_n1
            k+=1
        if k > maxIt
            warn("Max iterations reached")
            break
        end
    end
    return y_n1
end

function updateF!(F::Vector, f::Function, x::Vector, y::Vector, N::Int, h::Number)
    for i = 2:N
        F[i] = f(x[i], y[i], (y[i+1] - y[i-1])/(2h))
    end
    return F
end

function compJ(f::Function, x::Vector, y::Vector, h::Number, N::Number, M::Tridiagonal)
    J = Tridiagonal(zeros(N - 1), zeros(N), zeros(N - 1))
    ff(xys) = f(xys[1],xys[2],(xys[3]-xys[4])/(2h))
    J.d[1]  = M.d[1]
    J.du[1] = M.du[1]
    J.d[N]  = M.d[N]
    J.dl[N-1] = M.dl[N-1]
    for i = 2:N-1
        g = ForwardDiff.gradient(ff)([x[i],y[i],y[i+1],y[i-1]])
        J.d[i] = M.d[i] - g[2]
        J.dl[i-1] = M.dl[i-1] - g[3]
        J.du[i] = M.du[i] - g[4]
    end
    J
end