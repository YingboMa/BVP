using PlotlyJS
# y'' + p(x)*y' + q(x)*y = f(x)
# y'' + g/l*sin(x) = 0

N = 10000
a = 0.0
b = 10.0
i = 1:N-1
x = linspace(a,b,N)
h = (b-a) / N

g = 9.8 # m²/s
l = 1   # m

p(x) = 0
q(x) = 0
f(x) = - (g/l * sin(x))

a₀ = 1
b₀ = 0
a₁ = 1
b₁ = 0
c₀ = 0
c₁ = π

A = Array(Float64, N)
B = Array(Float64, N)
C = Array(Float64, N+1)
F = Array(Float64, N+1)

h² = 1 / h^2
h₂ = h*2

A[1] = a₀ - (3*b₀) / h₂
#B[1] = - b₀ / h₂    #
C[1] = 2*b₀ / h
F[1] = c₀

#A[N+1] = b₁ / h₂    # 
B[N] = a₁ + (3*b₁) / h₂
C[N+1] = - 2*b₁ / h
F[N+1] = c₁

for j in i
    p_i = p(i) / h₂
    A[i+1] = h² - p_i
    B[i] = h² + p_i       # Because of b₀ = b₁ = 0
    C[i+1] = -(2*h² - q(i))
    F[i+1] = f(x[i])
end

M = Tridiagonal([A[2:N], C[N+1];],
                [A[1], C[2:N], B[N];],
                [C[1], B[1:N-1];]
                )

Y = M\F

trace = scatter(;x=x, y=Y, mode="lines")
plot(trace)