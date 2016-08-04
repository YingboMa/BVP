using Gadfly
# y'' + p(x)*y' + q(x)*y = f(x)
#
# y'' + βx' + αx = 0

N = 50
a = 0.0
b = 10.0
i = 1:N-1
x = linspace(a, b, N)
h = (b-a) / N

g = 9.8 # m²/s
l = 1   # m

p(x) = 1.
q(x) = 1.
f(x) = 0 #- (g/l * sin(x))

a₀ = 1.
b₀ = 0.
c₀ = pi

a₁ = 0.
b₁ = 1.
c₁ = pi

@show a₀*b₁-a₁*b₀+a₀*a₁*(b-a)

A = Array(Float64, N+1)
B = Array(Float64, N+1)
C = Array(Float64, N+1)
F = Array(Float64, N+1)

h² = 1 / h^2
h₂ = h*2

A[1] = a₀ - (3*b₀) / h₂
B[1] = - b₀ / h₂
C[1] = 2*b₀ / h
F[1] = c₀

A[N+1] = b₁ / h₂
B[N+1] = a₁ + (3*b₁) / h₂
C[N+1] = - 2*b₁ / h
F[N+1] = c₁

for j in i
    p_i = p(x[i]) / h₂
    A[i+1] = h² - p_i
    B[i+1] = h² + p_i
    C[i+1] = -(2*h² - q(x[i]))
    F[i+1] = f(x[i])
end

if B[1] != 0
    ch = -B[1]/B[2]

    A[1] += ch*A[2]
    C[1] += ch*C[2]

    F[1] += ch*F[2]
end

if A[N+1] != 0
    ch = -A[N+1]/A[N]

    C[N+1] += C[N] * ch
    B[N+1] += C[N] * ch

    F[N+1] += F[N] * ch
end

M = Tridiagonal([A[2:N], C[N+1];],
                [A[1], C[2:N], B[N+1];],
                [C[1], B[2:N];]
                )

Y = M\F
Gadfly.plot(
    x=linspace(a, b, N+1), y=Y, Geom.point, Geom.line,
    Theme(panel_fill=colorant"black", default_color=colorant"orange"))
