### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 556a343a-c71e-11f1-bf7d-5518879398df
using Pkg

# ╔═╡ bdaf5791-83b8-4f74-93d9-7d634f12906f
Pkg.activate(".")

# ╔═╡ b1c3a51b-4d94-4b11-a1b2-14413010fd3a
using Revise

# ╔═╡ a50459be-747a-4bee-9780-bbbe87985a7a
using Polyhedron

# ╔═╡ 5a8ba776-4499-443c-96be-94bad9ed69da
begin
	A = [1 1;
     0 1]
    
    B = [2;
         1;;]
    C = [1 0]
    
    X = [0.8 0;
          0  1;
         -1  0;
          0  -1]
    
    U = [1.2;
         -1.5;;]
end

# ╔═╡ d8a582db-06e9-4ef1-8f3a-254335ca1344
sys = SLIT(A, B, C, SSF)

# ╔═╡ 952b7f5f-2c7c-454d-bc62-ae30549f48aa
ip_1 = is_pinvariant(sys, U, X)

# ╔═╡ a3118aa7-5661-45c4-9abe-2cf9ff68a1c6
traj = trajectory(sys, [1.25; 1;;], 10)

# ╔═╡ 4f0556d2-153d-4f51-8910-08b19a1185ca
sys_sof = SLIT(A, B, C, SOF)

# ╔═╡ 1fd24b26-761c-4775-bfed-0e2be072988c
ip = finding_L_pinvariant(sys, U, X)

# ╔═╡ 2a3a23b1-10dc-421c-99f4-38d06226b919
plot_poly(ip)

# ╔═╡ f9c8a1c1-8172-4ce2-946a-d7837d9b86e6
plot_trajectory(ip, traj)

# ╔═╡ 6ea49cb5-9e62-4390-94a5-bcbd30d236e8
plot_poly(ip_1, ip)

# ╔═╡ Cell order:
# ╠═556a343a-c71e-11f1-bf7d-5518879398df
# ╠═bdaf5791-83b8-4f74-93d9-7d634f12906f
# ╠═b1c3a51b-4d94-4b11-a1b2-14413010fd3a
# ╠═a50459be-747a-4bee-9780-bbbe87985a7a
# ╠═5a8ba776-4499-443c-96be-94bad9ed69da
# ╠═d8a582db-06e9-4ef1-8f3a-254335ca1344
# ╠═952b7f5f-2c7c-454d-bc62-ae30549f48aa
# ╠═2a3a23b1-10dc-421c-99f4-38d06226b919
# ╠═a3118aa7-5661-45c4-9abe-2cf9ff68a1c6
# ╠═f9c8a1c1-8172-4ce2-946a-d7837d9b86e6
# ╠═4f0556d2-153d-4f51-8910-08b19a1185ca
# ╠═1fd24b26-761c-4775-bfed-0e2be072988c
# ╠═6ea49cb5-9e62-4390-94a5-bcbd30d236e8
