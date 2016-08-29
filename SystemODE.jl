using Optim
using ForwardDiff

function BVP_sys(f::Function, a::Number,  b::Number,
                 c₀::Vector,
                 c₁::Vector;
                 N=Int((b-a)*15),
                 tol = 1e-5,  maxIt=Int(1e3))
    Ntotal = N*2
    i = 1:Ntotal-1
    h = (b-a) / N
    ty = typeof(c₀[1])
    x = collect(a:h:b)

    F = zeros(ty, Ntotal+2)
    y = zeros(ty, Ntotal+2)

    h₂ = 1 / 2h

    F[1]   = c₀[1]
    y[1]   = c₀[1]
    F[N+1] = c₀[2]
    y[N+1] = c₀[2]

    F[N]      = c₁[1]
    y[N]      = c₁[1]
    F[Ntotal] = c₁[2]
    y[Ntotal] = c₁[2]

    A = h₂*ones(N-1)

    M = Tridiagonal([-A, 0., 0., -A, 0.;],
                    [1., zeros(N-1), 1., 1., zeros(N-1), 1.;],
                    [0., A,  0.,  0., A;]
                    )
    G(y::Vector) = Base.norm(M*y - updateF!(F, f, x, y, N))
    # J = deepcopy(M)
    sol = Optim.minimizer(
        optimize(G, y, LBFGS(),
        OptimizationOptions(autodiff = false, show_trace = true,
        show_every = 5, iterations = maxIt, g_tol = tol)
        ))
    return linspace(a, b, N+1), sol
    # y_new = similar(y)
    # updateF!(F, f, x, y, N)
    # JG = 1.
    # k = 0
    # J = deepcopy(M)
    # while norm(JG) >= tol
    #     compJ!(J, f, x, y, h, N)
    #     G = M*y_n - F
    #     JG = J \ -G
    #     y_n1 = JG + y_n
    #     @assert all(isfinite(y_new))
    #     updateF!(F, f, x, y_new, N, h)
    #     y = deepcopy(y_new)
    #         k+=1
    #     if k > maxIt
    #         warn("Max iterations reached")
    #         break
    #     end
    # end
end

function updateF!(F::Vector, f::Function, x::Vector, y::Vector, N::Int)
    for i = 2:N
        F[i] = f(x[i], [y[i], y[N+i]])[1]
    end

    for i = N+2:2N
        F[i] = f(x[i-N], [y[i-N], y[i]])[2]
    end
    return F
end

function compJ!(J::Tridiagonal, f::Function, x::Vector, y::Vector, h::Number, N::Number)
    ff(xys) = f(xys[1],xys[2],(xys[3]-xys[4])/(2h))
    gg = t -> ForwardDiff.gradient(ff, t, Chunk{4}())
    for i = 2:N
        g = gg([x[i],y[i],y[i+1],y[i-1]])
        J.d[i] -= g[2]
        J.dl[i-1] -= g[3]
        J.du[i] -= g[4]
    end
    for i = N+2:2N
        g = gg([x[i],y[i],y[i+1],y[i-1]])
        J.d[i] -= g[2]
        J.dl[i-1] -= g[3]
        J.du[i] -= g[4]
    end
    J
end