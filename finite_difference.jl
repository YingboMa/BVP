function forwardDifference(F::Function, x::Vector, order::Int; StepSize=1e-7)
    Σ = zeros(x)
    for j = 1:length(x)
        for i = 0:order
            Σ[j] += (-1)^i * Base.binomial(order, i) * F(zero(1), x + (order - i)*StepSize)[j]
        end
    end
    return Σ / StepSize^order
end

function backwardDifference(F::Function, x::Vector, order::Int; StepSize=1e-7)
    Σ = zeros(x)
    for j = 1:length(x)
        for i = 0:order
            Σ[j] += (-1)^i * Base.binomial(order, i) * F(zeros(1), x - i*StepSize)[j]
        end
    end
    return Σ / StepSize^order
end

function centralDifference(F::Function, x::Vector, order::Int; StepSize=1e-7)
    Σ = zeros(x)
    for j = 1:length(x)
        for i = 0:order
            Σ[j] += (-1)^i * Base.binomial(order, i) * F(zeros(1), x + (order/2 - i)*StepSize)[j]
        end
    end
    return Σ / StepSize^order
end