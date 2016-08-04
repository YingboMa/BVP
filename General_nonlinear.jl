using PlotlyJS
# y'' = F(x, y, y''),    x ∈ (a, b)
# a₀y(a) + b₀y'(a) = c₀
# a₁y(b) + b₁y'(b) = c₁

N = 50
a = 0.0
b = 10.0
i = 1:N-1
x = linspace(a, b, N)
h = (b-a) / N

function f(t, y, yp)
    xprime = -sin(y)
end

a₀ = 1.
b₀ = 0.
c₀ = 0.

a₁ = 1.
b₁ = 0.
c₁ = 1.

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
    # F[i] = f(x) = x[i], y[i], (y[i + 1] - y[i])/(2h)
end

M = Tridiagonal([A[2:N], C[N+1];],
                [A[1], C[2:N], B[N+1];],
                [C[1], B[2:N];]
                )

function updateF!(F::Vector, f::Function, y, x, h)
    N = length(y)
    for i = 2:N - 1
        F[i] = f(x[i], y[i], (y[i+1] - y[i-1])/(2h))
    end
    return F
end

function compJ(f::Function, x, y)
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