### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ d7d2ec55-8de5-40f5-be93-18e1ffc1ef31
using Pkg

# ╔═╡ 45e5da98-e38c-4795-a6b6-b7f547b98b9d
Pkg.activate()

# ╔═╡ c11c50fb-fd54-4ace-a3b9-b53d94267287
Pkg.add("MAT")

# ╔═╡ fb34bda2-7b79-4de6-8bdc-e2ff6d102dc9
using Revise

# ╔═╡ a12754f5-6d7a-46ed-88ff-c3fe56080557
using PlutoUI

# ╔═╡ 6328bb51-e5c9-485c-b2c6-07b954faaf37
using Polyhedron

# ╔═╡ 36bf5979-f314-4f08-a192-a9b0c3eab65f
using LinearAlgebra

# ╔═╡ 2094eec6-7d71-4338-886c-722573924746
using MAT

# ╔═╡ f9a4c39e-7ef6-11f1-ba82-979f074de1d9
# ╠═╡ disabled = true
#=╠═╡
using Pkg
  ╠═╡ =#

# ╔═╡ b818e60c-86c9-4dcd-b86d-01bf740a3255
md"""
Definindo o sistema de interesse
"""

# ╔═╡ d8793666-4b32-498b-a60a-c0660d242f18
begin 
    A = [0.9666 0;
        0.0328 0.9666] 
    B = [0.1980;
        0.0033;;]
    
    C = [0 1]
    
    t = 0.5
    
    delta = 0.8
    
    A_exp = hcat(A, [0 0; 0 0])
    A_exp = vcat(A_exp, [0 -t 1 t; 0 0 0 delta])
    
    B_exp = vcat(B, [0;0])
    
    #E_exp = vcat(zeros(3), [1-delta])
    E_exp = [0; 0; 0; 1-delta;;]
    
    Sx = [1/15 0
          0 1/15;
         -1/15 0;
          0 -1/15;]
    Sv = [1/300;
          -1/300;]
    
    Sw = [1/50;
          -1/50;]
    
    Sx = vcat(Sx, zeros(4, 2))
    Sv = vcat(vcat(zeros(4, 1), Sv) , zeros(2, 1))
    Sw = vcat(zeros(6, 1), Sw)
    
    S = hcat(hcat(Sx, Sv), Sw)
    
    V = [1/4;
         -1/4;]
    
    # Fr para a referência
    R = [1/15;
          -1/15;;]
    d = 2
end

# ╔═╡ 2801e34e-df43-41c7-8096-dd1f39812b84
md"""
Definindo as matrixes que foram encontradas com o problema de otimização
"""

# ╔═╡ b6ceb5ce-ed29-4434-8b66-a7c07b0bec5e
begin
	F = [-0.0007192032411083585	2.260487767141468e-5	8.491472830817068e-5	-0.0007847902975979626;
	-0.004166202736310583	0.00022095520069825207	0.00038915477017592966	-0.0021081038352975465;
	-0.0011205093461843726	0.0017385509839190623	2.3163241772739214e-5	-0.0002117698787387944;
	-0.001352507917059543	-0.0030521771248577354	0.00029474513366686047	-0.0012470459289986766;
	0.044375168486375724	0.01230231671880651	-0.00474412596218661	0.008309627845662223;
	-0.002428916102814619	0.002423802629940089	0.00013777582714907246	-0.0005739661859077883;
	-0.024678537974451756	-0.009208254337837691	-0.0006051678451496056	0.0496092563242573;
	-0.0005447875238555245	0.005721888308299788	-0.0003045989847743259	0.0012997051333361728;
	-0.06404817009395036	-0.09767737304371042	0.007484084223159192	-0.0048393669518952045;
	-0.17726742031550285	-0.22264665841005715	0.017063371425514312	-6.973579432502995e-6;
	-0.03577003993250923	0.0351916188542822	0.0034433557522892664	-1.4114389559207376e-6;
	0.011050716932604152	-0.020241225507051104	0.00021435167121497276	0.011245719163074715;
	0.04296187859448264	-0.028498979540241495	-0.003682254790122036	-0.003455937816841924;
	0.03067193873052531	-0.021901707125223286	-0.0028610756948431303	0.0010429068572820124;
	-0.05909990451402138	-0.02959626916114162	0.004441994771655258	-0.06512980609329033;
	0.055367342628872177	0.07336737524228158	-0.005839908059192106	0.017440164449567967;
	-0.01845445862511624	-0.08410998782159876	0.0062816852165281075	-0.006315311147534877;
	0.0863080917260766	0.07282300820682695	-0.004334892682033207	-0.13956650420947583;
	-0.09647422646425953	-0.10976298685520446	0.008722697334035962	-0.004199651787151327;
	-0.003380865361867741	-0.01285441615312553	0.000857011052713737	-0.0022248922016910363;
	0.0025848565671949962	0.011696203139332421	-0.0007663021202077218	0.0020608212712371703;
	0.03063850754677806	-0.02327874159031448	-0.0029341592171167113	0.0015554841074633514;
	-0.05367366051405434	0.031423010971826246	-0.002458173643197253	-2.085730277022526e-6;
	-0.00555688527469321	-0.015785627877193908	0.0010873056513416974	-0.0027329838751536826]
	F, w = elimred(F, ones(size(F, 1)))
	F = vcat(F, -F)
end

# ╔═╡ 053b9d0b-e1e7-4bd9-813f-bde4fae8046a


# ╔═╡ e8b79a18-bd27-4773-8f3f-30372a91f9ec
T = poly_projection(F[:,3:4])

# ╔═╡ dd24d7d9-46e9-4bd6-8212-f455b11ceb05
plot_poly(T*F[:,1:2], T*ones(size(F, 1)))

# ╔═╡ 5839333a-3a5a-49f1-95a5-bcd65df98611
# ╠═╡ disabled = true
#=╠═╡
traj_ini = [[2.7388208921740618,	-14.040137368808697,	-300.04220168579127,	-15.011439973468073],
[2.1540512279393313,	-13.489585006201114,	-300.52785298812097	,-10.52915197877446],
[2.686181927381699,	-12.958312053209994,	-299.0476364744076,	-6.943321583019569]]
  ╠═╡ =#

# ╔═╡ fa66c678-3e1e-4016-861e-eb0179f793ab
G = [-0.09043589303665445	-0.08656210844107828	0.005122875252774146	0.2380944797007607]


# ╔═╡ b0060dbb-ea29-498b-b17b-e42a6ba3078d
md"""
Agora vou testar com os vertices do Poliedro de condicoes iniciais admissiveis
"""

# ╔═╡ 67f685f5-816b-4d4f-a0ec-0c12336fccd7
nx = size(A_exp, 1)

# ╔═╡ e5235f94-2b34-4e57-8a98-0e74b9eca789
begin
    n = size(A_exp, 1)
    
    Z = zeros(n, n)
    Zf = zeros(size(F, 1), n)
    
    Fe = [F Zf Zf;
          Zf F Zf;
          Zf Zf F]
    
    Ee = [E_exp; [0;0;0;0]; [0;0;0;0]]
    
    Ae = [A_exp  Z  B_exp*G;
            I(n)    Z  Z;
           Z      I(n)  Z]
end

# ╔═╡ 5ba9d134-301e-4a52-b39f-950df23687d5
begin
	N = size(Fe, 1)
	ones_N = ones(N)
	
	v1 = calcular_v(Fe, Ee, R)
	v2 = calcular_v(Fe, Ae * Ee, R)
	
	l = ones_N - v1
	m = ones_N - v1 - v2
	
	init_cond_w = vcat(ones_N, vcat(l, m))
	init_cond_F = [Fe; Fe*Ae; Fe*Ae*Ae]
	init_cond_F, init_cond_w = elimred(init_cond_F, init_cond_w)
end

# ╔═╡ 4ac93b77-7567-430f-82fe-d2e2ed75a893
# ╠═╡ disabled = true
#=╠═╡
init_cond_F, init_cond_w = admissable_initCond(A_exp, B_exp*G, F, d, w; fixed_d=true)
  ╠═╡ =#

# ╔═╡ a97a764b-0b36-426e-99b9-69b7fbe1f543
md"""
Vou tentar multiplas trajetórias para a mesma referência mas com condições iniciais diferentes
"""

# ╔═╡ 6a50db38-84c1-4eb9-8ef9-c6dee43ce933
md"""
Considerando o vetor de condicoes iniciais nessa ordem:

[ x[k], x[k-1], ... , x[k-d] ] 
"""

# ╔═╡ 87aa2fea-d2ed-4e47-b8b0-fdf30062caa8
v = get_extVert_tuple(init_cond_F, init_cond_w, 4)

# ╔═╡ 68b3b178-67ae-4d13-b37c-9b7b17e12220
v_sorted = y = sort(v; by = x -> x[1][1])

# ╔═╡ 2621005a-c65d-4962-9243-18d0daf7d2a3
md"""
Todos os vertices:
"""

# ╔═╡ 6d9edc1a-63fd-4c2e-a69b-6e8fd007d828
plot_vertices(F, p[1] for p in v_sorted)

# ╔═╡ e8eb024b-a10c-44be-8f6a-1349270e3cc4
# leste 6, leste 89, norte [362, 364], sul 37, oeste 447, leste 39

# ╔═╡ 11ffd751-41d1-44b4-a8ed-1f533882a1ca
md"""
fui procurando por todos os vertices no poliedro
"""

# ╔═╡ 09598248-1558-4a1b-9673-43847b6d9de1
@bind i Slider(1:size(v_sorted, 1), show_value=true)

# ╔═╡ dbe8ea2f-6f3c-496c-9c52-1a464ce92a55
plot_vertices(F, [v_sorted[i][1]])

# ╔═╡ 0c9e7151-4906-4e5e-a274-7f3d2a36dd07
md"""
Vamos pegar 4 pontos de partida
"""

# ╔═╡ 0372b717-c18c-4e30-b42e-839c0811e2c4
begin
	plt_sul, traj_sul = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[94108]), passos=300)
	plt_sul, traj_sul = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_sul, passos=300)
end

# ╔═╡ bf396529-f045-4b00-aac7-9fdc1f02047f
#=╠═╡
plt_sul
  ╠═╡ =#

# ╔═╡ fe099f74-a943-4fa9-8b10-db54714855a9
begin
	plt_norte, traj_norte = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[31824]), passos=300)
	plt_norte, traj_norte = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_norte, passos=300)
end

# ╔═╡ 3694f734-e34a-4e6f-b6dc-abc3c73bee6c
#=╠═╡
plt_norte
  ╠═╡ =#

# ╔═╡ 8d7b357e-451a-4098-985d-89c884092fed
begin
	_, traj_leste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[107852]), passos=300)
	plt_leste, traj_leste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_leste, passos=300)
end

# ╔═╡ 4bdd0dc6-6f5b-4be2-b41e-aa8f73839bf2
#=╠═╡
plt_leste
  ╠═╡ =#

# ╔═╡ 18c35444-7fbb-495f-823c-f7e3d3c252be
begin
	_, traj_oeste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[9536]), passos=300)
	plt_oeste, traj_oeste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_oeste, passos=300)
end

# ╔═╡ ae317296-62b3-40ac-99cb-fc6cba3cb3f7
#=╠═╡
plt_oeste
  ╠═╡ =#

# ╔═╡ 67aeadc8-a22b-44f5-a45c-7782c59342ab
md"""
Todas trajetórias observadas no mesmo plot
"""

# ╔═╡ 0b61751d-646d-40f1-99d3-d61aabeea274
plot_trajectories(F, [traj_norte, traj_sul, traj_leste, traj_oeste])

# ╔═╡ 2ee46796-0afc-4b05-854c-abb5b32c2f57
md"""
Vizualizando demais estados
"""

# ╔═╡ f490a3ea-180b-41a8-a850-2d7052494742
plot_expanded_state_trajectory([traj_norte, traj_sul, traj_leste, traj_oeste])

# ╔═╡ b0674ad3-758c-4aad-91b2-dc0147691697
plot_v_state_trajectory([traj_norte, traj_sul, traj_leste, traj_oeste])

# ╔═╡ f29a1973-7a18-45b4-8946-bdf36afb4592
md"""
Passando as matrizes para um arquivo .mat
"""

# ╔═╡ 6bf29487-e6af-41c0-9899-42593ce22649
begin
	K = [-0.005895724889024788	-0.004852451763849023	0.00023072302830382636	0.02523718111724989;
	-0.0028750606375070767	0.0006813146818040654	-0.0002771916753272453	0.0007857547976951455;
	0.09175481815071201	0.14138426919747357	-0.011897449486050867	-0.2275783948857918;
	0.003230800225057837	0.003708901096689136	-0.0003192354887556665	-0.006480028934478674]
	u_norte = [G*x for x in traj_norte]
	u_norte = [elem[1] for elem in u_norte]
	u_sul = [G*x for x in traj_sul]
	u_sul = [elem[1] for elem in u_sul]
	u_leste = [G*x for x in traj_leste]
	u_leste = [elem[1] for elem in u_leste]
	u_oeste = [G*x for x in traj_oeste]
	u_oeste = [elem[1] for elem in u_oeste]
end

# ╔═╡ eb96fe20-7b7a-4004-9e27-73eb4293e873
plot_along_k([u_leste, u_oeste, u_norte, u_sul])

# ╔═╡ b7e3dcbf-b17a-4ba3-8e43-7f9e1e4327db
maximum(abs.(u_leste))

# ╔═╡ 0d119224-7aee-44c3-a2c0-2a01f11d3a46
maximum(abs.(u_oeste))

# ╔═╡ 5fc556df-9203-4a8b-ab5d-e6c683f574bd
maximum(abs.(u_norte))

# ╔═╡ 30129ad8-a6d1-45a2-821c-48a5fdb4e3ea
maximum(abs.(u_sul))

# ╔═╡ dc0468c1-24bd-48b0-ac1e-aa0c4825fa86
begin
	file = matopen("results_logs/artigo_cba_new_cond.mat", "w")
	
	write(file, "traj_norte", traj_norte)
	write(file, "traj_sul", traj_sul)
	write(file, "traj_leste", traj_leste)
	write(file, "traj_oeste", traj_oeste)

	write(file, "F", F)
	write(file, "G", G)

	write(file, "u_norte", u_norte)
	write(file, "u_sul", u_sul)
	write(file, "u_leste", u_leste)
	write(file, "u_oeste", u_oeste)
	write(file, "K", K)

	write(file, "A", A)
	write(file, "B", B)
	write(file, "C", C)
	
	close(file)
end

# ╔═╡ Cell order:
# ╠═f9a4c39e-7ef6-11f1-ba82-979f074de1d9
# ╠═45e5da98-e38c-4795-a6b6-b7f547b98b9d
# ╠═fb34bda2-7b79-4de6-8bdc-e2ff6d102dc9
# ╠═a12754f5-6d7a-46ed-88ff-c3fe56080557
# ╠═6328bb51-e5c9-485c-b2c6-07b954faaf37
# ╟─b818e60c-86c9-4dcd-b86d-01bf740a3255
# ╠═d8793666-4b32-498b-a60a-c0660d242f18
# ╟─2801e34e-df43-41c7-8096-dd1f39812b84
# ╠═b6ceb5ce-ed29-4434-8b66-a7c07b0bec5e
# ╠═053b9d0b-e1e7-4bd9-813f-bde4fae8046a
# ╟─e8b79a18-bd27-4773-8f3f-30372a91f9ec
# ╠═dd24d7d9-46e9-4bd6-8212-f455b11ceb05
# ╟─5839333a-3a5a-49f1-95a5-bcd65df98611
# ╠═fa66c678-3e1e-4016-861e-eb0179f793ab
# ╟─b0060dbb-ea29-498b-b17b-e42a6ba3078d
# ╟─67f685f5-816b-4d4f-a0ec-0c12336fccd7
# ╠═36bf5979-f314-4f08-a192-a9b0c3eab65f
# ╠═e5235f94-2b34-4e57-8a98-0e74b9eca789
# ╠═5ba9d134-301e-4a52-b39f-950df23687d5
# ╠═4ac93b77-7567-430f-82fe-d2e2ed75a893
# ╟─a97a764b-0b36-426e-99b9-69b7fbe1f543
# ╟─6a50db38-84c1-4eb9-8ef9-c6dee43ce933
# ╠═87aa2fea-d2ed-4e47-b8b0-fdf30062caa8
# ╠═68b3b178-67ae-4d13-b37c-9b7b17e12220
# ╟─2621005a-c65d-4962-9243-18d0daf7d2a3
# ╠═6d9edc1a-63fd-4c2e-a69b-6e8fd007d828
# ╠═e8eb024b-a10c-44be-8f6a-1349270e3cc4
# ╟─11ffd751-41d1-44b4-a8ed-1f533882a1ca
# ╠═09598248-1558-4a1b-9673-43847b6d9de1
# ╠═dbe8ea2f-6f3c-496c-9c52-1a464ce92a55
# ╟─0c9e7151-4906-4e5e-a274-7f3d2a36dd07
# ╠═0372b717-c18c-4e30-b42e-839c0811e2c4
# ╟─bf396529-f045-4b00-aac7-9fdc1f02047f
# ╠═fe099f74-a943-4fa9-8b10-db54714855a9
# ╟─3694f734-e34a-4e6f-b6dc-abc3c73bee6c
# ╠═8d7b357e-451a-4098-985d-89c884092fed
# ╟─4bdd0dc6-6f5b-4be2-b41e-aa8f73839bf2
# ╠═18c35444-7fbb-495f-823c-f7e3d3c252be
# ╟─ae317296-62b3-40ac-99cb-fc6cba3cb3f7
# ╟─67aeadc8-a22b-44f5-a45c-7782c59342ab
# ╠═0b61751d-646d-40f1-99d3-d61aabeea274
# ╟─2ee46796-0afc-4b05-854c-abb5b32c2f57
# ╠═f490a3ea-180b-41a8-a850-2d7052494742
# ╠═b0674ad3-758c-4aad-91b2-dc0147691697
# ╠═f29a1973-7a18-45b4-8946-bdf36afb4592
# ╠═d7d2ec55-8de5-40f5-be93-18e1ffc1ef31
# ╟─c11c50fb-fd54-4ace-a3b9-b53d94267287
# ╠═2094eec6-7d71-4338-886c-722573924746
# ╠═6bf29487-e6af-41c0-9899-42593ce22649
# ╠═eb96fe20-7b7a-4004-9e27-73eb4293e873
# ╠═b7e3dcbf-b17a-4ba3-8e43-7f9e1e4327db
# ╠═0d119224-7aee-44c3-a2c0-2a01f11d3a46
# ╠═5fc556df-9203-4a8b-ab5d-e6c683f574bd
# ╠═30129ad8-a6d1-45a2-821c-48a5fdb4e3ea
# ╠═dc0468c1-24bd-48b0-ac1e-aa0c4825fa86
