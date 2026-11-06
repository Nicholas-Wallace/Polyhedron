### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 8791174e-da1d-11f1-9cd9-41cdb38b9c55
using Pkg

# ╔═╡ 8981e562-141b-4f13-8c85-25acee1ae09a
Pkg.activate(".")

# ╔═╡ e94a447e-351a-463b-97aa-c72551e14f2e
using Revise

# ╔═╡ 2e71cca9-3dce-4f54-99a6-72183b06ec53
using Polyhedron

# ╔═╡ 9449376a-daff-4bab-bcaa-fccdb3f4fa05
begin
	A_2 = [1.2 0.2; -0.4 0.6]
	Ad = [-0.3 -0.2; 0.4 0.2]
	C = [1 1]
end

# ╔═╡ 1f96048b-4eb1-40f0-93f7-8d6b0aa348fd
sys_2 = SystemDelayed(A_2, Ad, C, 1, SSF)

# ╔═╡ e5d2bf50-ae04-45ab-a36c-37df9e50b2dc
ip_2 = finding_L_pinvariant(sys_2, [10 10; 10 -10; -10 10; -10 -10], [-1 -1; 2 1; 1 1; -2 -1;])

# ╔═╡ cb738856-db4b-4956-9c55-080e85887b5f
U = [10 10; 10 -10; -10 10; -10 -10]

# ╔═╡ 65f33442-d868-4d9e-aec8-ce84a0a90299
F = float([-1 -1; 2 1; 1 1; -2 -1;])

# ╔═╡ 560daede-3c3b-472d-9316-f04337e1f33b
plot_poly(ip_2)

# ╔═╡ 780d7457-603b-4a86-9a44-113665bc3b99
traj = trajectory(sys_2, [2 0; -3 1], 100)

# ╔═╡ f83e256d-47c2-4c65-8e3c-dfee858f08e7
plot_trajectory(ip_2, traj[:, 2:end])

# ╔═╡ 622fd96e-1273-407b-8ec7-e8cb13d10fe6
ip = is_pinvariant(sys_2, F)

# ╔═╡ e8477d78-2a17-480d-81e2-c785981ac1d0


# ╔═╡ Cell order:
# ╠═8791174e-da1d-11f1-9cd9-41cdb38b9c55
# ╠═8981e562-141b-4f13-8c85-25acee1ae09a
# ╠═2e71cca9-3dce-4f54-99a6-72183b06ec53
# ╠═e94a447e-351a-463b-97aa-c72551e14f2e
# ╠═9449376a-daff-4bab-bcaa-fccdb3f4fa05
# ╠═1f96048b-4eb1-40f0-93f7-8d6b0aa348fd
# ╠═e5d2bf50-ae04-45ab-a36c-37df9e50b2dc
# ╠═cb738856-db4b-4956-9c55-080e85887b5f
# ╠═65f33442-d868-4d9e-aec8-ce84a0a90299
# ╠═560daede-3c3b-472d-9316-f04337e1f33b
# ╠═780d7457-603b-4a86-9a44-113665bc3b99
# ╠═f83e256d-47c2-4c65-8e3c-dfee858f08e7
# ╠═622fd96e-1273-407b-8ec7-e8cb13d10fe6
# ╠═e8477d78-2a17-480d-81e2-c785981ac1d0
