### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 54a3ea38-6b27-11f1-1ec0-bfdb353c2801
using Pkg

# ╔═╡ 948bafad-c0b3-451d-896c-36830fa3fed4
Pkg.activate()

# ╔═╡ 3fb6f754-bec1-4820-b6db-29a5ff605edb
using Revise

# ╔═╡ dfb4041a-5908-4db1-b628-339c88697df7
using Polyhedron

# ╔═╡ 3e262867-828d-4a9d-a50e-ba2b95af8b62
begin
        A = [1 1; 0 1]
        B = [2; 1;;]
        C = [1 0]
        X = [0.8 0; 0 1; -1 0; 0 -1]
        U = [1.2; -1.5;;]
end

# ╔═╡ dc080723-e220-48e2-9365-df47ef603fad
ex_capitulo = is_pinvariant(A, B, C, U, X)

# ╔═╡ 388eb623-4dff-4b39-a0f1-a47bd4e06b99
ex_capitulo["lambda"]

# ╔═╡ a7d83b1b-2cf3-4706-8245-fb670b514dec
plt = plot_poly(X, ones(size(X, 1)))

# ╔═╡ 081e9d9c-a614-42e4-9939-badc67332d82
x0 = [[-1.0, 1.0]]

# ╔═╡ 6a871c50-9595-4a57-9f4e-cde3a476c31c
A

# ╔═╡ 3601a766-c79e-421f-9830-4db0c0e69437
A*x0[1]

# ╔═╡ fbc792d2-f65b-4641-89f1-a54eb4dfb28b
passos = 15

# ╔═╡ 9f32c252-075e-4e65-b6ba-1cf8b581b68c
pipe_trajectory(A, B, ex_capitulo["F"], X, passos=passos)

# ╔═╡ ed1ce78d-b2a2-49f1-90e2-1d302f15e8ea
begin
	A_ = [1.2 0.2; -0.4 0.6]
	Ad_ = [-0.3 -0.2; 0.4 0.3]
	
	X_ = [-1 -1; 2 1]
end

# ╔═╡ f4220962-484d-441e-a1c2-76e14ccb8342
ex_artigo = is_pinvariant(A_, Ad_, X, d=1; symetric=true)

# ╔═╡ 6b6ecbad-6342-43e0-8d5e-d3bf57b415c6


# ╔═╡ 3b9ba8f1-0b51-42e3-b9b6-c4a04d7d1b62
is_pinvariant(A, B, C, U, X, SOF=true)

# ╔═╡ Cell order:
# ╟─54a3ea38-6b27-11f1-1ec0-bfdb353c2801
# ╟─948bafad-c0b3-451d-896c-36830fa3fed4
# ╟─3fb6f754-bec1-4820-b6db-29a5ff605edb
# ╠═dfb4041a-5908-4db1-b628-339c88697df7
# ╠═3e262867-828d-4a9d-a50e-ba2b95af8b62
# ╠═dc080723-e220-48e2-9365-df47ef603fad
# ╠═388eb623-4dff-4b39-a0f1-a47bd4e06b99
# ╠═a7d83b1b-2cf3-4706-8245-fb670b514dec
# ╠═081e9d9c-a614-42e4-9939-badc67332d82
# ╠═6a871c50-9595-4a57-9f4e-cde3a476c31c
# ╠═3601a766-c79e-421f-9830-4db0c0e69437
# ╠═fbc792d2-f65b-4641-89f1-a54eb4dfb28b
# ╠═9f32c252-075e-4e65-b6ba-1cf8b581b68c
# ╠═ed1ce78d-b2a2-49f1-90e2-1d302f15e8ea
# ╠═f4220962-484d-441e-a1c2-76e14ccb8342
# ╠═6b6ecbad-6342-43e0-8d5e-d3bf57b415c6
# ╟─3b9ba8f1-0b51-42e3-b9b6-c4a04d7d1b62
