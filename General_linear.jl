# y'' + p(x)*y' + q(x)*y = f(x)
#
# a₀y(a) + b₀y(a) = c₀ 
# a₁y(b) + b₁y'(b) = c₁
#
function BVP_linear(p::Function, q::Function, f::Function,
                    a::Number, b::Number,
                    a₀::Number,  b₀::Number,  c₀::Number,
                    a₁::Number,  b₁::Number,  c₁::Number;
                    N=Int((b-a)*20)
                    )

    if (a₀*b₁-a₁*b₀)+a₀*a₁*(b-a) == 0
        println("This system is not unique!")
    else
        i = 1:N-1
        x = linspace(a, b, N)
        h = (b-a) / N

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
        return linspace(a, b, N+1), Y
    end
end