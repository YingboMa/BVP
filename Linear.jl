using BandedMatrices
# y'' + p(x)*y' + q(x)*y = f(x)
#
# a₀y(a) + b₀y(a) = c₀ 
# a₁y(b) + b₁y'(b) = c₁
#
function BVP_linear(order::Int, p::Function, f::Function,
                    a::Number, b::Number,
                    a₀::Vector,   a₁::Vector;
                    N=Int((b-a)*50)
                    )

#===================================================================
The general idea of this solver is to solve a linear system
                        M * y = F
===================================================================#

    x = linspace(a, b, N)
    h = (b-a) / N
    T = typeof(a₀[end])
    A = Array(T, N+1, 6)

    A[1, :]     = zeros(6)
    A[end, :]   = zeros(6)
    A[1, 1]     = a₀[1]
    A[end, 5]   = a₁[1]
    A[1, end]   = a₀[end]
    A[end, end] = a₁[end]

#===================================================================
                    Construct matrix A
===================================================================#

    for j in 2:N
        p_i = p(x[j])
        # TODO: change d to 1:5 and try to access memory linearly
        for d in 1:5
            if j == 2
                A[j, d] = FPSF[order, d] / h^order
            elseif j == N
                A[j, d] = -FPSF[order, 6-d] / h^order
            else
                A[j, d] = FPSC[order, d] / h^order
            end

            for c in 1:order-1
                if j == 2
                    A[j, d] += FPSF[c, d] * p_i[c] / h^c
                elseif j== N
                    A[j, d] -= FPSF[c, 6-d] * p_i[c] / h^c
                else
                    A[j, d] += FPSC[c, d] * p_i[c] / h^c
                end
            end

        end
        A[j, 3] += p_i[end]
        A[j, 6] =  f(x[j])  # This is the matrix F
    end

    for d in 1:5
        for c in 1:order-1
            A[1,   d] += FPSF[c, d]   * a₀[c+1]  / h^c
            A[end, d] -= FPSF[c, 6-d] * a₁[c+1]  / h^c
        end
    end

#===================================================================
Use row operation to eliminate extra elements to form a banded matrix
===================================================================#

        if A[1, 5] != zero(T) && A[2, 5] != zero(T)
            ch = -A[1, 5]/A[2, 5]
            A[1, :] += ch * A[2, :]

        elseif A[2, 5] == zero(T)
            C = A[1, :]
            A[1, :] = A[2, :]
            A[2, :] = C
        end
        
        if A[2, 5] != zero(T) && A[3, 5] != zero(T)
            ch = -A[2, 5]/A[3, 5]
            A[2, :] += ch * A[3, :]

        elseif A[3, 5] == zero(T)
            C = A[2, :]
            A[2, :] = A[3, :]
            A[3, :] = C
        end

        if A[end, 1] != zero(T) && A[end-1, 1] != zero(T)
            ch = -A[end, 1]/A[end-1, 1]
            A[end, :] += A[end-1, :] * ch

        elseif A[end-1, 1] == zero(T)
            C = A[end-1, 1]
            A[end-1, 1] = A[end, 1]
            A[end, 1] = C
        end

        if A[end-1, 1] != zero(T) && A[end-2, 1] != zero(T)
            ch = -A[end-1, 1]/A[end-2, 1]
            A[end-1, :] += A[end-2, :] * ch

        elseif A[end-2, 1] == zero(T)
            C = A[end-2, 1]
            A[end-2, 1] = A[end-1, 1]
            A[end-1, 1] = C
        end

#===================================================================
            Construct matrix M and F from matrix A
===================================================================#

    # M = bzeros(N+1, N+1, 4, 4)

    # M[1, 1:4] = A[1, 1:4][:]
    # M[2, 2:6] = A[2, 1:5][:]
    # M[end-1, end-5:end-1]  = A[end-1, 1:5][:]
    # M[end, end-3:end]      = A[end, 2:5][:]

    # for i in 3:N-1
    #     M[i, i-2:i+2] = A[i, 1:5][:]
    # end

    M = bzeros(N+1, N+1, 3, 3)
    M[1, 1:4] = A[1, 1:4][:]
    M[2, 2:5] = A[2, 1:4][:]
    M[end-1, end-4:end-1]  = A[end-1, 2:5][:]
    M[end, end-3:end]      = A[end, 2:5][:]
    for i in 3:N-1
        M[i, i-2:i+2] = A[i, 1:5][:]
    end

#===================================================================
                    Solve the linear system
                            M * y = F
===================================================================#
    y = M \ A[:,6][:]
    return M, A, linspace(a, b, N+1), y
end

# Five-Point Stencil Central
FPSC =
[
[1/12   -2/3    0.  2/3 -1/12]
[-1/12   4/3  -5/2  4/3 -1/12]
[-1/2     1.    0.   -1.  1/2]
[1.      -4.    6.   -4.    1]
]
# Five-Point Stencil F&B
FPSF =
[
[-25/12  4.      -3.      4/3    -1/4]
[35/12  -26/3   19/2    -14/3   11/12]
[-5/2    9.      -12.      7.    -3/2]
[1.     -4.       6.      -4.      1.]
]