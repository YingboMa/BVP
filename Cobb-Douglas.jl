using PlotlyJS
using ODE
include("First_order.jl")

function Cobb(t, k)
    c = (α*k^(α-1)-ρ-δ)
    k^α - c - σ*k
end

α = 0.33
δ = 0.1
function init_guess()

end
tmp = BVP_nonlinear(Cobb, 0., 500., 1., 30.7, N=750, maxIt=100, y=);
sol = ode78(Cobb, 1., [0:100...])
layout = Layout(autosize=false, width=800, height=500, title="Cobb-Douglas")
trace1 = scatter(;x=tmp[1], y=tmp[2], name="BVP solver")
trace2 = scatter(;x=sol[1], y=sol[2], name="ODE78")
PlotlyJS.plot([trace1, trace2], layout)