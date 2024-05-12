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

# ╔═╡ ae9d3cc2-0ff0-11ef-0035-1d0e787b9280
begin
	import Pkg
	Pkg.activate(".")
	Pkg.add([
		"PlutoUI", "StatsModels",
		"GLM", "DataFrames", "StatsBase"
	])
end

# ╔═╡ 265d7845-a8bb-4cf2-9143-5fb173273ba2
using PlutoUI, Plots, StatsPlots, Distributions, Random, Turing, StatsModels, GLM, DataFrames, StatsBase

# ╔═╡ b6a4cd5f-3c52-4920-a686-bed33434c25c
sigma = @bind σ Slider(0.0:.1:1.0, show_value=true);

# ╔═╡ 21f4ee72-b2df-4a18-a181-7841c5790b0a
md"### Noise $(σ)"

# ╔═╡ 40578b4f-f335-47ba-8e70-8c4aa55eff64
sigma

# ╔═╡ ce121f90-eb39-40fc-a197-82266cbeaf3b
begin
	N_samples = 30
	if σ != 0
		ϵ = rand(Normal(0, σ), N_samples)
	else
		ϵ = zeros(N_samples)
	end
end;

# ╔═╡ fedd0a36-4de5-4f6f-99ae-37a93273968f
begin
	Random.seed!(42)
	x = -2.0:.01:2.0
	β = .4
	α = 1
	y = @. β*x + α
	test_points = rand(Uniform(-2, 2), N_samples);
	obs = @. β * test_points + α + ϵ;
end;

# ╔═╡ e9611cdf-deaf-4e99-898c-f142adf6a17f
data = DataFrame(x=test_points, y=obs); first(data, 5)

# ╔═╡ 050f629b-70ba-45f8-b103-595bcf9d0a6d
ols = lm(@formula(y~x), data); ols

# ╔═╡ 04362ea9-edac-4eb5-96d3-74e2ce926129
begin
	plot(x, y, label="Data Model\nα=$(α) β=$(β)", color=:black)
	α̂, β̂ = coef(ols)
	ŷ = @. β̂ * x + α̂
	plot!(x, ŷ, label="Model Estimate\nα̂=$(round(α̂,digits=2)) β̂=$(round(β̂,digits=2))", color=:blue, linestyle=:dash)
	scatter!(test_points, obs, label="Observed Data", color=:green)
	xlims!(-2, 2)
	ylims!(-1, 3)
	xlabel!("x")
	ylabel!("y")
	title!("OLS Example")
end

# ╔═╡ Cell order:
# ╟─21f4ee72-b2df-4a18-a181-7841c5790b0a
# ╟─40578b4f-f335-47ba-8e70-8c4aa55eff64
# ╟─04362ea9-edac-4eb5-96d3-74e2ce926129
# ╟─050f629b-70ba-45f8-b103-595bcf9d0a6d
# ╟─e9611cdf-deaf-4e99-898c-f142adf6a17f
# ╟─fedd0a36-4de5-4f6f-99ae-37a93273968f
# ╟─b6a4cd5f-3c52-4920-a686-bed33434c25c
# ╟─ce121f90-eb39-40fc-a197-82266cbeaf3b
# ╟─265d7845-a8bb-4cf2-9143-5fb173273ba2
# ╟─ae9d3cc2-0ff0-11ef-0035-1d0e787b9280
