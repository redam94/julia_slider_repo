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

# ╔═╡ 42f55973-0cc6-4837-9b01-ac72337a1fe1
begin
	import Pkg;
	Pkg.activate(".")
	Pkg.add([
		"Colors", "LambertW", "SpecialFunctions", "NLsolve",
	])
end

# ╔═╡ 6c37b344-1b74-4a78-a7f3-394e62ca8a8a
using Colors, LambertW, SpecialFunctions, NLsolve;

# ╔═╡ d6458a2c-0405-4529-80d1-a01b7193997e
using Plots, PlutoUI, Random, Distributions;

# ╔═╡ 7d14f706-cd3e-11ec-07c6-75ef1918e381
md"
# S Shaped Curves
## How they come about and their usefulness.
Physical processes are often limited to available quantities, i.e., you can't sell more dog food than you have, or to more people than those that exist. This fundemental limitation introduces a non-linearity that would need to be captured by a model to avoid unrealistic conclusions. This is where s shaped curves come into play.

There are 3 main properties a good choice of s curve should have. 
* They should start at ``0`` (No spend in a media channel ⟹ no effect)
* They should be monotonically increasing (An increase in spend ⟹ an increase in effect)
* They should saturate (The amount of effect has an upper limit)
!!! note
	There are a few other properties that should also be considered but not at the expense of the 3 above: 
	* Ease of estimating a unique set of parameters (Parameters that can not be uniquely determined we will produce poor estimates for contirbutions)
	* Ease of implimentation (If we can't readily compute the s curve we can't use it in the model in the first place)
	* Have some basis in the underlying physical process. (This is mainly to limit the space of functions without introducing too much bias)
I will be looking at three main categories of S curves below.
* Current double exponential functions.
* Functions that arise from physical sciences.
* Cummulative distribution functions.
---
# 1. Double Exponential S Curves
"

# ╔═╡ 4a82e1e6-1fbe-4aad-bc1e-96ce4423c1cb
md"""
The Double Exponential S curves have 2 forms:
* The standard S Curve: ``S(x| \beta, \alpha) = \beta^{\alpha^x}``
* The S Origin: ``S_{origin}(x| \beta, \alpha) = \beta^{\alpha^x}-\beta``
!!! note
	* ``\alpha`` acts like it stretches the x axis.
	* ``\beta`` translates the x axis.
"""

# ╔═╡ b7874e3d-fde1-4e10-a538-2547d7a1faea
md"""
## The standard S Curve ``\beta^{\alpha^x}``
The standard s curve is relatively easy to compute only relying on exponentiation. It saturates at 1. It has 2 shape parameters ``\alpha`` and ``\beta``. The derivatives are also relatively simple to calculate and allow closed form solutions to the optimal point, saturation begins, and the breakthrough point. 
!!! warning "Beware!!!"
	* This choice of curve does not start at zero
	* Double exponential functions rarely come about in physical process
## The S Origin Curve ``\beta^{\alpha^x} - \beta``
This curve solves the starting point issue with the standard s curve. This curve starts at 0. Like the standard S curve this is also easy to calculate. Also has 2 shape parameters. 
!!! note 
	* Unlike the standard S curve this curve saturates at ``1-\beta``.
"""

# ╔═╡ 0ec2fa7b-0f42-406e-a5c9-fa4ff95d37a5
md"""
# 2. Curves from biology!
These function are rooted in biology and physical chemistry. They tend to focus on the equalibrium solutions to chemical reactions or population growth. 
## The Hill Function ``\theta(x|K,n) = \frac{x^n}{x^n+K^n}``
This promising S shaped curve function comes from a simplified model of ligand binding. As with ligand binding, an advertisement has to "react" with a person in order to have an effect on revenue. Also there is a finite quantity of people those advertisements can "react" with.Here we can think of advertisements as ligands and people as the macromolecules. Then the hill function can be thought of as the fraction of a population affected by the advertising channel!
!!! note
	* The hill function starts at 0 for n>0!
	* K is the value at which 50% of the response is captured
	* n can be loosely thought of as how multiple views of the channel affect a persons conversion rate
	* Realistically n>1
"""

# ╔═╡ d6eba7dd-c7d3-4e77-bd94-45175bba93f5
md"""
* ``K``: $(@bind K Slider(1:100, default=20))
* ``n``: $(@bind n Slider(1.1:.1:5, default=2))
"""

# ╔═╡ 60cb5e1b-a3bc-4dfa-90b6-f93810c71190
md"""
#### ``K``: $(K), ``n``: $(n)
"""

# ╔═╡ cf54bcc5-3e94-40b3-8261-4a00c8f12e1e
md"""
## The Logistic Function: ``L(x| x_0, k) = \frac{1}{1+e^{-k(x-x_0)}}``
This function typically appears in population growth models. Growth begins geometrically and falls off as x gets larger.
"""

# ╔═╡ fab95351-5115-44ed-aead-2401eb1a864c
md"""
!!! warning "Beware!!!"
	* This function does not start at zero
	* This function is very sensitive to choice of k
"""

# ╔═╡ 390061cc-ac32-4c0d-8b15-cc24b077d4b9
md"""
* ``x_0``: $(@bind x₀ Slider(1:100, default=20))
* ``k``: $(@bind k Slider(0.001:.001:1, default=.1))
"""

# ╔═╡ 78d5266c-3013-43c6-9d4a-93bcb55b792b
md"""
#### ``x_0``: $(x₀), ``k``: $(k)
"""

# ╔═╡ 994ad38d-122e-4c4c-bd12-091698e92184
md"""
# 3. Cumulative Probability Distributions
We are looking at population effects. If we imagine the s curve as being the proportion of a suseptable population affected by an impression level at of at least ``x``. The captured effect would be proportional to ``P(X\le x)``.
"""

# ╔═╡ 710cb03c-06bc-4536-8463-6cea0c7ff9a4
md"""
## Weibull CDF: ``Weibull(x|b, k) = 1-e^{-bx^k}``
"""

# ╔═╡ 319002bd-7ad5-4cb6-b406-217122c57464
md"""
* ``b``: $(@bind b Slider(.1:.1:5, default=3))
* ``k``: $(@bind k_w Slider(1:.1:5, default=2))
"""

# ╔═╡ d794899f-c3eb-47d9-bb95-57eb13f68257
md"""
#### ``b``: $(b), ``k``: $(k_w)
"""

# ╔═╡ cc1c3404-922c-4a28-9093-8ac5ee8a3034
md"""
## Gamma CDF: ``Gamma(x; \alpha, \beta) = \frac{1}{\Gamma(\alpha)}\gamma(\alpha, \beta x)``
"""

# ╔═╡ 3ddb51d9-e2c1-4079-b5ce-7e210b57b963
md"""
* ``\alpha``: $(@bind α_g Slider(.1:.1:10, default=3))
* ``\beta``: $(@bind β_g Slider(-1:.01:1.5, default=1))
"""

# ╔═╡ 1347ab68-5ad5-43cd-8c8e-7726f5adf7ef
md"""
### ``\alpha =`` $(α_g), ``\beta =`` $(10.0^-β_g)
"""

# ╔═╡ 4561a7bd-6a11-452f-8505-16676b387285
function ddxgammarev(F, v)
	b_g = 10.0^(-β_g)
	x = v[1]
	F[1] = (α_g-1)*gamma(α_g)*gamma_inc(α_g, b_g*abs(x))[1] - gamma(α_g+1) * gamma_inc(α_g+1, b_g*abs(x))[1]
end

# ╔═╡ a9e9f1ff-3ade-4d75-814c-0ab6ee28de0c
function sCurveDoubleExponential(x; α=.95, β=10.0^-3)
	return β^(α^x)
end

# ╔═╡ 78523353-e1b0-4fb9-91a7-dc757cacc59f
function sCurveDoubleExponentialOrigin(x; α=.95, β=10.0^-3)
	return β^(α^x) - β
end

# ╔═╡ 30b8b00d-e484-4a2f-8a14-bb6c6757e682
md"""
* ``\alpha``: $(@bind α Slider(.80:.01:.99, default=.9))
* ``\beta``: $(@bind β Slider(1:.1:13, default=3))
* Function: $(@bind eef Select([sCurveDoubleExponential => "S_Curve", sCurveDoubleExponentialOrigin => "S_Origin"]))
"""

# ╔═╡ 5d7cfb58-4692-407c-991c-3540b0f6cf71
md"""
#### ``\alpha``: $(α), ``\beta``: $(10.0^(-β))
"""

# ╔═╡ 6804bd6a-bf2c-4432-95f1-c1d765bfa036
function hill(x; K=20, n=2)
	return x^n/(x^n+K^n)
end

# ╔═╡ c44e4a5a-86d8-47ca-b1fa-af9737389da1
function logistic(x; x₀=1, k=2)
	return 1/(1+ℯ^(-k*(x-x₀)))
end

# ╔═╡ 9d2d096d-c009-4520-a6f2-9ab2854c97c2
function Weibull(x; b=1.0, k=1.0)
	return 1-exp(-b*(x^k))
end

# ╔═╡ da261ae9-8b1e-43dd-bfa6-8bd7f961e8ab
colors = Dict("DarkBlue"=>RGB(0, 0, 139/255), "LightBlue"=>colorant"#008B8B", "Red"=>colorant"#8B0000", "Gold"=>colorant"#8B8B00"); 

# ╔═╡ f0ae0619-b095-490b-8ece-38e1df3e2765
let
	f(x) = eef(x, α=α, β=10.0^(-β))
	true_optimal = -100
	if β>=1.8
		true_optimal = lambertw((-1/(β*log(10))),-1)/log(α)
	end
	if Symbol(eef)==:sCurveDoubleExponentialOrigin && β>1
		true_optimal = lambertw((1-10.0^(-β))/(-β*log(10)), -1)/log(α)
	end
	optimal_point = -log(β*log(10))/log(α)
	optimal_value = f(optimal_point)
	sb_point = log((3-sqrt(5))/(2*β*log(10)))/log(α)
	sb_value = f(sb_point)
	x = 0:.01:120
	name = "S_Curve"
	if Symbol(eef) == :sCurveDoubleExponentialOrigin
		name = "S_Origin"
	end
	p = plot(x, f.(x), ylims=(0, 1), xlims=(0, 120), label=name, linecolor=RGB(0, 0, 139/255))
	plot!([optimal_point, optimal_point], [0, 1], label="Optimal begins/ optimal mROAS $(round(optimal_point, digits=2))", linestyle=:dash, linecolor=colors["LightBlue"])
	scatter!([optimal_point], [optimal_value], label="", markercolor=colors["LightBlue"])
	plot!([sb_point, sb_point], [0, 1], linestyle=:dash, linecolor=colors["Gold"], label="Saturation Begins $(round.(sb_point, digits=2))")
	scatter!([sb_point], [sb_value], label="")
	if true_optimal>0
		plot!([true_optimal, true_optimal], [0, 1], label="Optimal ROI $(round(true_optimal, digits=2))", linestyle=:dash, linecolor=:gray)
		scatter!([true_optimal], [f(true_optimal)], label="", markercolor=:gray)
	end
	p
end

# ╔═╡ 05f72657-dfbe-442b-a3ef-a1d7f8c7e78a
let
	f(x) = hill(x, K=K, n=n)
	true_optimal = 0
	optimal_point = 0
	sb_point = 0
	if n>1
		optimal_point = ((n-1)/(n+1))^(1/n)*K
		true_optimal = (n-1)^(1/n)*K
		sb_point = (((3.0)^.5*abs(n^2-n)+2*n^2-2)/(n^2+3*n+2))^(1.0/n)*K
	end
	
	x = 0:.01:120
	name = "Hill Function"
	p = plot(x, f.(x), ylims=(0, 1), xlims=(0, 120), label=name, linecolor=RGB(0, 0, 139/255))
	plot!([optimal_point, optimal_point], [0, 1], label="Optimal Point $(round(optimal_point, digits=2))", linestyle=:dash, linecolor=colors["LightBlue"])
	scatter!([optimal_point], [f(optimal_point)], markercolor=colors["LightBlue"], label="")
	plot!([sb_point, sb_point], [0, 1], label="Saturation begins $(round(sb_point, digits=2))", linestyle=:dash, linecolor=colors["Gold"])
	scatter!([sb_point], [f(sb_point)], markercolor=colors["Gold"], label="")
	plot!([true_optimal, true_optimal], [0, 1], label="Optimal ROI $(round(true_optimal, digits=2))", linestyle=:dash, linecolor=:gray)
	scatter!([true_optimal], [f(true_optimal)], markercolor=:gray, label="")
	p
end

# ╔═╡ 7c2a0863-93fb-49c4-9a73-1cb3625f7df6
let
	f(x) = logistic(x, x₀=x₀, k=k)
	true_optimal = 0
	optimal_point = 0
	sb_point = 0
	
	optimal_point = x₀
	true_optimal = -lambertw(-k/(exp(k*x₀)), -1)/k
	sb_point = -log(10/(5+sqrt(5))-1)/k+x₀
	
	
	x = 0:.01:120
	name = "Hill Function"
	p = plot(x, f.(x), ylims=(0, 1), xlims=(0, 120), label=name, linecolor=RGB(0, 0, 139/255))
	plot!([optimal_point, optimal_point], [0, 1], label="Optimal Point $(round(optimal_point, digits=2))", linestyle=:dash, linecolor=colors["LightBlue"])
	scatter!([optimal_point], [f(optimal_point)], markercolor=colors["LightBlue"], label="")
	plot!([sb_point, sb_point], [0, 1], label="Saturation begins $(round(sb_point, digits=2))", linestyle=:dash, linecolor=colors["Gold"])
	scatter!([sb_point], [f(sb_point)], markercolor=colors["Gold"], label="")
	plot!([true_optimal, true_optimal], [0, 1], label="Optimal ROI $(round(true_optimal, digits=2))", linestyle=:dash, linecolor=:gray)
	scatter!([true_optimal], [f(true_optimal)], markercolor=:gray, label="")
	p
end

# ╔═╡ 6f88f527-b6bc-44e2-b318-5436ea60a076
let
	b_w = 10.0^(-b)
	f(x) = Weibull(x, b=10.0^(-b), k=k_w)
	true_optimal = 0
	optimal_point = 0
	sb_point = 0
	if k_w>1
		optimal_point = 10.0^(b/k_w)*((k_w-1)/k_w)^(1/k_w)
		sb_point=2^(-1/k_w)*((b_w*k_w*sqrt(5*k_w^2-6*k_w+1)-b_w*(3-3*k_w)*k_w)/(b_w^2*k_w^2))^(1/k_w)
		true_optimal = ((-k_w*lambertw(-1/(exp(1.0/k_w)*k_w), -1)-1)/(b_w*k_w))^(1/k_w)
	end
	
	
	x = 0:.1:100
	name = "Wiebull CDF"
	p = plot(x, f.(x), ylims=(0, 1), xlims=(0, 100), label=name, linecolor=RGB(0, 0, 139/255))
	plot!([optimal_point, optimal_point], [0, 1], label="Optimal Point $(round(optimal_point, digits=2))", linestyle=:dash, linecolor=colors["LightBlue"])
	scatter!([optimal_point], [f(optimal_point)], markercolor=colors["LightBlue"], label="")
	plot!([sb_point, sb_point], [0, 1], label="Saturation begins $(round(sb_point, digits=2))", linestyle=:dash, linecolor=colors["Gold"])
	scatter!([sb_point], [f(sb_point)], markercolor=colors["Gold"], label="")
	plot!([true_optimal, true_optimal], [0, 1], label="Optimal ROI $(round(true_optimal, digits=2))", linestyle=:dash, linecolor=:gray)
	scatter!([true_optimal], [f(true_optimal)], markercolor=:gray, label="")
	p
end

# ╔═╡ 32a73eb8-4fa4-4d74-9018-291c22144b06
let
	b_g = 10.0^(-β_g)
	gamma = Gamma(α_g, 1/b_g)
	f(x) = cdf(gamma, x)
	true_optimal = 0
	optimal_point = 0
	sb_point = 0
	if α_g>=1
		optimal_point = (α_g-1)/b_g
		sb_point = ((α_g-1)+sqrt(α_g-1))/(b_g)
		true_optimal= nlsolve(ddxgammarev, [sb_point;]).zero[1]
	end
	#	optimal_point = 10.0^(b/k_w)*((k_w-1)/k_w)^(1/k_w)
	#	sb_point=2^(-1/k_w)*((b_w*k_w*sqrt(5*k_w^2-6*k_w+1)-b_w*(3-3*k_w)*k_w)/(b_w^2*k_w^2))^(1/k_w)
	#	true_optimal = ((-k_w*lambertw(-1/(exp(1.0/k_w)*k_w), -1)-1)/(b_w*k_w))^(1/k_w)
	#end
	
	
	x = 0:.1:120
	name = "Gamma CDF"
	p = plot(x, f.(x), ylims=(0, 1), xlims=(0, 120), label=name, linecolor=RGB(0, 0, 139/255))
	plot!([optimal_point, optimal_point], [0, 1], label="Optimal Point $(round(optimal_point, digits=2))", linestyle=:dash, linecolor=colors["LightBlue"])
	scatter!([optimal_point], [f(optimal_point)], markercolor=colors["LightBlue"], label="")
	plot!([sb_point, sb_point], [0, 1], label="Saturation begins $(round(sb_point, digits=2))", linestyle=:dash, linecolor=colors["Gold"])
	scatter!([sb_point], [f(sb_point)], markercolor=colors["Gold"], label="")
	plot!([true_optimal, true_optimal], [0, 1], label="Optimal ROI $(round(true_optimal, digits=2))", linestyle=:dash, linecolor=:gray)
	scatter!([true_optimal], [f(true_optimal)], markercolor=:gray, label="")
	p
end

# ╔═╡ Cell order:
# ╟─7d14f706-cd3e-11ec-07c6-75ef1918e381
# ╟─4a82e1e6-1fbe-4aad-bc1e-96ce4423c1cb
# ╟─30b8b00d-e484-4a2f-8a14-bb6c6757e682
# ╟─5d7cfb58-4692-407c-991c-3540b0f6cf71
# ╟─f0ae0619-b095-490b-8ece-38e1df3e2765
# ╟─b7874e3d-fde1-4e10-a538-2547d7a1faea
# ╟─0ec2fa7b-0f42-406e-a5c9-fa4ff95d37a5
# ╟─d6eba7dd-c7d3-4e77-bd94-45175bba93f5
# ╟─60cb5e1b-a3bc-4dfa-90b6-f93810c71190
# ╟─05f72657-dfbe-442b-a3ef-a1d7f8c7e78a
# ╟─cf54bcc5-3e94-40b3-8261-4a00c8f12e1e
# ╟─fab95351-5115-44ed-aead-2401eb1a864c
# ╟─390061cc-ac32-4c0d-8b15-cc24b077d4b9
# ╟─78d5266c-3013-43c6-9d4a-93bcb55b792b
# ╟─7c2a0863-93fb-49c4-9a73-1cb3625f7df6
# ╟─994ad38d-122e-4c4c-bd12-091698e92184
# ╟─710cb03c-06bc-4536-8463-6cea0c7ff9a4
# ╟─319002bd-7ad5-4cb6-b406-217122c57464
# ╟─d794899f-c3eb-47d9-bb95-57eb13f68257
# ╟─6f88f527-b6bc-44e2-b318-5436ea60a076
# ╟─cc1c3404-922c-4a28-9093-8ac5ee8a3034
# ╟─3ddb51d9-e2c1-4079-b5ce-7e210b57b963
# ╟─1347ab68-5ad5-43cd-8c8e-7726f5adf7ef
# ╟─32a73eb8-4fa4-4d74-9018-291c22144b06
# ╟─4561a7bd-6a11-452f-8505-16676b387285
# ╟─a9e9f1ff-3ade-4d75-814c-0ab6ee28de0c
# ╟─78523353-e1b0-4fb9-91a7-dc757cacc59f
# ╟─6804bd6a-bf2c-4432-95f1-c1d765bfa036
# ╟─c44e4a5a-86d8-47ca-b1fa-af9737389da1
# ╟─9d2d096d-c009-4520-a6f2-9ab2854c97c2
# ╟─da261ae9-8b1e-43dd-bfa6-8bd7f961e8ab
# ╠═6c37b344-1b74-4a78-a7f3-394e62ca8a8a
# ╠═d6458a2c-0405-4529-80d1-a01b7193997e
# ╟─42f55973-0cc6-4837-9b01-ac72337a1fe1
