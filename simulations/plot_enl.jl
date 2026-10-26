### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 9bdeb39e-cd95-11f1-a2c0-1f9fecca428e
using Pkg

# ╔═╡ f7de87a8-f6d4-4bb8-875d-d2768a96db9e
Pkg.activate(".")

# ╔═╡ d42003c7-124a-4a81-b5ca-a465e29057e8
using Revise

# ╔═╡ 85e64ebb-6e8e-4bf2-80a9-4d592cf10db0
using Polyhedron

# ╔═╡ bd2b7321-834c-456b-85ab-b1be5f4f08fd
S = [1/3 1/3;
	-1/2 1/2
	1/2 -1/2
	-1/3 -1/3;
	-1/4 0;
	1/4 0]

# ╔═╡ 05ff9c37-ac6a-41e8-b1ec-41a7cdaab39e
plot_poly(S, ones(6))

# ╔═╡ 71eb9bee-ee52-4817-a3fc-b3bbf67e6e59
X = [1 0;
	-1 0;
	0 1;
	0 -1]

# ╔═╡ f872df35-814a-41f5-8a43-f5bc5600a89a
plot_poly(X, ones(4))

# ╔═╡ 8389120d-4023-4390-b612-a9736b468947
ip_int, ip_ext = enl_simples(X, S)

# ╔═╡ 56b841cd-d16c-4faa-95f8-4d1d02f1e755
plot_poly(ip_int, ip_ext, title="Representação enl(.) 2D")

# ╔═╡ 0ca5494b-0049-48d7-ae9d-c70baa174383
1/12

# ╔═╡ Cell order:
# ╠═9bdeb39e-cd95-11f1-a2c0-1f9fecca428e
# ╠═f7de87a8-f6d4-4bb8-875d-d2768a96db9e
# ╠═d42003c7-124a-4a81-b5ca-a465e29057e8
# ╠═85e64ebb-6e8e-4bf2-80a9-4d592cf10db0
# ╠═bd2b7321-834c-456b-85ab-b1be5f4f08fd
# ╠═05ff9c37-ac6a-41e8-b1ec-41a7cdaab39e
# ╠═71eb9bee-ee52-4817-a3fc-b3bbf67e6e59
# ╠═f872df35-814a-41f5-8a43-f5bc5600a89a
# ╠═8389120d-4023-4390-b612-a9736b468947
# ╠═56b841cd-d16c-4faa-95f8-4d1d02f1e755
# ╠═0ca5494b-0049-48d7-ae9d-c70baa174383
