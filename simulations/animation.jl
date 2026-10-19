### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ e432788a-03ca-4eb8-8901-0392b8e55796
using Pkg

# ╔═╡ b3135adf-ec4d-4826-825c-f5821faacc17
Pkg.activate(".")

# ╔═╡ 8cac86a0-23c7-4e42-8ced-9801145f0c17
using Revise

# ╔═╡ 6c2c31dd-e12e-4635-bb6c-8b35eddad057
using MAT

# ╔═╡ b4762019-eb77-4316-8574-b3cefbf88e09
using Polyhedron

# ╔═╡ 2e769adc-c80c-4a50-a18b-fa37f54fa83f
dados = matread("results_logs/artigo_cba_new_cond.mat")

# ╔═╡ da1beccc-ae17-44a1-ae0c-c5faa4ca1ac6
keys(dados)

# ╔═╡ 0cee0caf-b81e-49f3-a8b0-fb4d5b33dba1
begin
	F = dados["F"]
	traj_sul = dados["traj_sul"]
	traj_norte = dados["traj_norte"]
	traj_leste = dados["traj_leste"]
	traj_oeste = dados["traj_oeste"]

	traj_sul_zero = dados["traj_sul_zero"]
	traj_norte_zero = dados["traj_norte_zero"]
	traj_leste_zero = dados["traj_leste_zero"]
	traj_oeste_zero = dados["traj_oeste_zero"]
end

# ╔═╡ 1ec5cbb1-309d-4c13-8a62-2859aa2c00f3
traj_sul

# ╔═╡ b04a0d8a-605a-4029-8c15-9c34da43a3be
begin
	gif_trajectories(F, traj_sul, dest="traj_sul")
	gif_trajectories(F, traj_norte, dest="traj_norte")
	gif_trajectories(F, traj_leste, dest="traj_leste")
	gif_trajectories(F, traj_oeste, dest="traj_oeste")
	gif_trajectories(F, traj_sul_zero, dest="traj_sul_zero")
	gif_trajectories(F, traj_norte_zero, dest="traj_norte_zero")
	gif_trajectories(F, traj_leste_zero, dest="traj_leste_zero")
	gif_trajectories(F, traj_oeste_zero, dest="traj_oeste_zero")
end

# ╔═╡ Cell order:
# ╠═e432788a-03ca-4eb8-8901-0392b8e55796
# ╠═b3135adf-ec4d-4826-825c-f5821faacc17
# ╠═6c2c31dd-e12e-4635-bb6c-8b35eddad057
# ╠═8cac86a0-23c7-4e42-8ced-9801145f0c17
# ╠═b4762019-eb77-4316-8574-b3cefbf88e09
# ╠═2e769adc-c80c-4a50-a18b-fa37f54fa83f
# ╠═da1beccc-ae17-44a1-ae0c-c5faa4ca1ac6
# ╠═0cee0caf-b81e-49f3-a8b0-fb4d5b33dba1
# ╠═1ec5cbb1-309d-4c13-8a62-2859aa2c00f3
# ╠═b04a0d8a-605a-4029-8c15-9c34da43a3be
