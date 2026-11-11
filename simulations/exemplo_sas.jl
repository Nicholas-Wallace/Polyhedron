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
	A = [1.2 0.2; -0.4 0.6]
	Ad = [-0.3 -0.2; 0.4 0.2]
	C = [1 1]
end

# ╔═╡ 1f96048b-4eb1-40f0-93f7-8d6b0aa348fd
sys = SystemDelayed(A, Ad, C, 1, SSF)

# ╔═╡ ba0062eb-7f96-48e8-a396-90faaf9a412b
sys_d10 = SystemDelayed(A, Ad, C, 10, SSF)

# ╔═╡ e46f047f-f5f0-48c1-baa3-83aff72c9738
ip = pinvariant_poly(sys, [-1 -1; 2 1; 1 1; -2 -1], enl=state_constraint_enl!)

# ╔═╡ 34156e26-460f-40c2-b8ee-349c6c361fa8
plot_poly(ip)

# ╔═╡ fd9c2ad1-09f8-4487-80e6-c5ad319bd167
ip_2 = pinvariant_poly(sys, [-1 -1; 2 1; 1 1; -2 -1], t=4, F_init=ip.F)

# ╔═╡ 1cf268c5-b83b-4d6c-86b0-26a0e50b0128
plot_poly(ip_2)

# ╔═╡ 39f68b78-3f30-4409-a123-227fc95ec4f2
is_pinvariant(sys_d10, [-1 -1; 2 1; 1 1; -2 -1])

# ╔═╡ 36ef58cd-7bba-4482-af9e-2eac8fb87d76
ip_init_cond = pinvariant_poly(sys, [-1 -1; 2 1; 1 1; -2 -1], t=4, F_init=ip_2.F)

# ╔═╡ c34c0124-b5c7-4e6d-8a2f-3336a21bb022
plot_poly(ip_init_cond)

# ╔═╡ 0fae5baf-ec1a-4d18-aec2-2536d46fae9a
begin
    A_t = [0.9666 0;
        0.0328 0.9666] 
    B_t = [0.1980;
        0.0033;;]
    C_t = [0 1]
    X = [1/15 0
          0 1/15;
         -1/15 0;
          0 -1/15;]
end

# ╔═╡ 8bf6ff87-d489-4c1a-beec-ab03c023c6c6
sys_t = SLIT(A_t, B_t, C_t, SSF)

# ╔═╡ 2944508a-299e-4650-aefd-6ab08773fecf
ip_t = pinvariant_poly(sys_t, X, enl=state_constraint_enl!)

# ╔═╡ a0df670c-5877-4d0f-a10a-0bbf88f3517b
plot_poly(pinvariant_poly(sys_t, X, enl=state_constraint_enl!, F_init=ip_t.F))

# ╔═╡ bfe6553f-c5d6-4a0a-8c94-c00715a86491
sys_t_d = SystemDelayed(A_t, B_t, C_t, 1, SSF)

# ╔═╡ 6d1a8a56-bac9-4563-a904-d7e1277e839e
model = pinvariant_poly(sys_t_d, X, return_model=true)

# ╔═╡ ee069dbd-c5a0-4373-880e-c30eb5555ec5
constraint_by_name(model, "contractiave_pi_H")

# ╔═╡ aab98b53-7247-4886-892a-3c6b47aa922a
plot_poly(X, ones(4))

# ╔═╡ Cell order:
# ╠═8791174e-da1d-11f1-9cd9-41cdb38b9c55
# ╠═8981e562-141b-4f13-8c85-25acee1ae09a
# ╠═2e71cca9-3dce-4f54-99a6-72183b06ec53
# ╠═e94a447e-351a-463b-97aa-c72551e14f2e
# ╠═9449376a-daff-4bab-bcaa-fccdb3f4fa05
# ╠═1f96048b-4eb1-40f0-93f7-8d6b0aa348fd
# ╠═ba0062eb-7f96-48e8-a396-90faaf9a412b
# ╠═e46f047f-f5f0-48c1-baa3-83aff72c9738
# ╠═34156e26-460f-40c2-b8ee-349c6c361fa8
# ╠═fd9c2ad1-09f8-4487-80e6-c5ad319bd167
# ╠═1cf268c5-b83b-4d6c-86b0-26a0e50b0128
# ╠═39f68b78-3f30-4409-a123-227fc95ec4f2
# ╠═36ef58cd-7bba-4482-af9e-2eac8fb87d76
# ╠═c34c0124-b5c7-4e6d-8a2f-3336a21bb022
# ╠═0fae5baf-ec1a-4d18-aec2-2536d46fae9a
# ╠═8bf6ff87-d489-4c1a-beec-ab03c023c6c6
# ╠═2944508a-299e-4650-aefd-6ab08773fecf
# ╠═a0df670c-5877-4d0f-a10a-0bbf88f3517b
# ╠═bfe6553f-c5d6-4a0a-8c94-c00715a86491
# ╠═6d1a8a56-bac9-4563-a904-d7e1277e839e
# ╠═ee069dbd-c5a0-4373-880e-c30eb5555ec5
# ╠═aab98b53-7247-4886-892a-3c6b47aa922a
