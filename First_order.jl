using ForwardDiff
# y' = F(x, y),    x ∈ (a, b)
# y(a) = c₀
# y(b) = c₁
function BVP_nonlinear(f::Function, a::Number,  b::Number,
                       c₀::Number,
                       c₁::Number;
                       N=Int((b-a)*15),         y=0,
                       tol = 1e-5,  maxIt=Int(1e3))
    N = (N-1)*2
    i = 1:N-1
    h = (b-a) / N
    ty = typeof(c₀)
    x = [a:h:b;]
    if y==0
        slope = (c₁ - c₀) / (b-a)
        y = convert(Array{ty}, slope*x+c₀)
    end

    F = zeros(ty, N+1)

    h₂ = 1 / h*2

    F[1] = c₀
    y[1] = c₀
    F[N+1] = c₁
    y[N+1] = c₁

    M = Tridiagonal([-h₂*ones(N-1), 0.;],
                    [1., zeros(N-1), 1.;],
                    [0., h₂*ones(N-1);]
                    )
    y_n = similar(y, N+1)
    y_n = copy(y)
    updateF!(F, f, x, y, N)
    JG = 1.
    k = 0
    J = Tridiagonal(zeros(N), zeros(N+1), zeros(N))
    K = 0
    while norm(JG) >= tol
        @show k+=1
        J = compJ(J, f, x, y, N, M)
        G = M*y_n - F
        JG = J \ G
        # getodd!(JG)
        y_n = -JG + y
        @show find(x->x<0, y_n)
        updateF!(F, f, x, y_n, N)
        y = y_n
        if k > maxIt
            warn("Max iterations reached")
            break
        end
    end
    return linspace(a, b, N+1)[1:2:end], y[1:2:end], M, F, J
end

# function getodd!(A)
#     A[2:2:end] = 0
#     A
# end

function updateF!(F::Vector, f::Function, x::Vector, y::Vector, N::Int)
    for i = 2:N
        F[i] = f(x[i], y[i])
    end
    return F
end

function compJ(J::Tridiagonal, f::Function, x::Vector, y::Vector, N::Number, M::Tridiagonal)
    J = copy(M)
    ff(xys) = f(xys[1],xys[2])
    gg = t -> ForwardDiff.gradient(ff, t, Chunk{2}())
    for i = 2:N
        g = gg([x[i],y[i]])
        J.d[i] -= g[2]
    end
    J
end