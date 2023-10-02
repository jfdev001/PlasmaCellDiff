using Plots
# https://discourse.julialang.org/t/deactivate-plot-display-to-avoid-need-for-x-server/19359/6
Plots.default(show=false)

x = 0:0.1:3
y = exp.(x)

p = plot(x, y)
show(p, )
