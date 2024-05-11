### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 39fb9179-dd46-430e-8429-c9bb9650f7fd
begin
	import Pkg
	Pkg.activate(".")
	Pkg.add([
		"PlutoUI",
		"Random",
		"Turing",
		"Plots",
		"StatsPlots",
		"Distributions"
	])
end

# ╔═╡ 5404fef6-0ceb-11ef-0c2e-056033a575d4
using PlutoUI, Random, Turing, Plots, Distributions, StatsPlots

# ╔═╡ ac89180b-343a-4b91-92fa-7666b8da3fd3
begin
	w1 = (@bind s1 Slider(-1:.01:1));
	w2 = (@bind s2 Slider(-1:.01:1));
	w3 = (@bind s3 Slider(-1:.01:1));
end;

# ╔═╡ 8bd28075-76fb-464b-8a65-30f00aa4c58d
md"$(w1) w1 = $(s1)"

# ╔═╡ 3b1dac0c-af04-4292-a023-6da4e3ea3112
md"$(w2) w2 = $(s2)"

# ╔═╡ 57149f6a-a758-4eb2-a6db-db7fd2607195
md"$(w3) w3 = $(s3)"

# ╔═╡ f2c980d0-10f4-40ec-b7a1-214a8f0ffb1f
begin
	season_weight = @bind sw Slider(-1:0.01:1)
	trend_weight = @bind tw Slider(-1:0.01:1)
	constant = @bind c Slider(-1:.5:10)
end;

# ╔═╡ b84979a5-e22f-490d-a772-f9830a4fdba1
md"$(season_weight) sw = $(sw)"

# ╔═╡ 06da06b7-e983-4b18-89a0-b15095faf2ff
md"$(trend_weight) tw = $(tw)"

# ╔═╡ e1a40733-cc33-4b88-9b5d-ad384e7dd39e
md"$(constant) c = $(c)"

# ╔═╡ 811af22c-a9f7-45ff-84ec-4e610d5e2d01
let
	t = 0:1:155
	season =  sum([w .* sin.(f*2*π*t / 52.18) for (w, f) in [(s1, 1), (s2, 2), (s3, 3)]
	])
	noise = rand(Normal(0, .02), length(t))
	covar = .2 .* t ./ 52 .+ rand(Normal(0, .1), length(t)) .- .2 .* season
	obs_covar = exp.(covar)
	y = exp.(sw .* season .+ tw .* t ./ 52.0 .+ c .+ noise .+ .2 * covar)
	sc = scatter(t, y)
	pl = scatter(obs_covar, y)
	plot(sc, pl)
end

# ╔═╡ 180cea24-d49d-425b-872b-cf4861410627


# ╔═╡ a69d104c-84d2-4155-903c-24b36c33c7b1
@model function linear(X, y)
	obs, n = shape(X)
end

# ╔═╡ Cell order:
# ╠═39fb9179-dd46-430e-8429-c9bb9650f7fd
# ╠═5404fef6-0ceb-11ef-0c2e-056033a575d4
# ╟─ac89180b-343a-4b91-92fa-7666b8da3fd3
# ╟─8bd28075-76fb-464b-8a65-30f00aa4c58d
# ╟─3b1dac0c-af04-4292-a023-6da4e3ea3112
# ╟─57149f6a-a758-4eb2-a6db-db7fd2607195
# ╟─f2c980d0-10f4-40ec-b7a1-214a8f0ffb1f
# ╟─b84979a5-e22f-490d-a772-f9830a4fdba1
# ╟─06da06b7-e983-4b18-89a0-b15095faf2ff
# ╟─e1a40733-cc33-4b88-9b5d-ad384e7dd39e
# ╠═811af22c-a9f7-45ff-84ec-4e610d5e2d01
# ╠═180cea24-d49d-425b-872b-cf4861410627
# ╠═a69d104c-84d2-4155-903c-24b36c33c7b1
