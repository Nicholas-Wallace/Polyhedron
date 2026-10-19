### A Pluto.jl notebook ###
# v0.20.28

using Markdown
using InteractiveUtils

# ╔═╡ 45e5da98-e38c-4795-a6b6-b7f547b98b9d
#=╠═╡
Pkg.activate()
  ╠═╡ =#

# ╔═╡ c11c50fb-fd54-4ace-a3b9-b53d94267287
#=╠═╡
Pkg.add("MAT")
  ╠═╡ =#

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
    F = [
        5.94566678455784e-5	5.9246039809246274e-5	-2.9920033957354245e-6	-1.2846373703718837e-5;
        5.689263130924019e-5	5.352477082043005e-5	-3.122205451871746e-6	-3.053238781691081e-5;
        2.929895314019125e-5	2.56367722917431e-5	-1.6542706868975309e-6	-2.4643730195993224e-5;
        3.2481023226809313e-5	2.655548020378451e-5	-1.3784556324959852e-6	-2.618151623485802e-5;
        -8.486684880810956e-7	-1.2717360016303154e-6	5.929427276088602e-9	-2.2986527656584684e-6;
        0.0027406397181953595	-0.16795860853742575	0.005052715637588139	1.4852663923278467e-8;
        -5.075608773582398e-6	-8.700810065480087e-6	2.200267422292824e-7	-1.4673582117634128e-5;
        2.7925164813274055e-9	1.0079270985576511e-9	-0.0033332727090779695	-5.253504821672977e-9;
        -2.4683248695573825e-6	-3.313395524007632e-6	1.3155243534754887e-7	-2.8556084117376687e-6;
        6.710326582287261e-6	8.449495095276045e-6	-3.745929850227469e-7	5.118155787338942e-6;
        -3.127050720942536e-7	5.006815995900341e-6	7.404972305778283e-8	2.32467556142376e-5;
        -2.052323698117072e-6	-7.018153482904333e-6	2.0621088614066433e-7	-1.8140842873361475e-5;
        3.621040111617345e-7	4.393497476330485e-6	1.2182910475349818e-8	1.7216871554962623e-5;
        -1.9337774467431847e-6	-7.325388017017013e-6	1.629239198154696e-7	-2.065910992115748e-5;
        1.3365192456435227e-5	1.2445313259199236e-5	-6.306907133964221e-7	-5.715344554166864e-6;
        3.1979215018482246e-6	4.3167693511362844e-6	-3.955265124985359e-7	-6.220450990952028e-7;
        -0.37945274628167014	-0.40744733858562043	0.018672667452048625	6.125281875385242e-8;
        -4.288821560241404e-8	-8.820147699314112e-10	-9.557858734682634e-7	0.06666646756831834;
        2.7190630376599423e-5	2.246886618450714e-5	-1.2863608897510603e-6	-2.3502806277301522e-5;
        -1.385135949534879e-7	5.712231594514103e-6	-6.58524630645277e-8	2.284366014823067e-5;
        -2.2505145763523486e-7	7.596563489273779e-7	2.220874560444653e-8	4.332787245059746e-6;
        0.013839048411046693	0.0148192200259226	-0.0007230260965961522	-9.524675801923828e-5;
        1.978895574569648e-6	-3.0223692511159758e-6	1.1233734491599425e-8	-1.8951026901221304e-5;
        -1.2902914331425894e-7	-1.884844258075871e-6	3.6404547034731274e-8	-6.659713987262626e-6;
    ]
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
#=╠═╡
v = get_extVert_tuple(init_cond_F, init_cond_w, 4)
  ╠═╡ =#

# ╔═╡ 68b3b178-67ae-4d13-b37c-9b7b17e12220
#=╠═╡
v_sorted = y = sort(v; by = x -> x[1][1])
  ╠═╡ =#

# ╔═╡ 2621005a-c65d-4962-9243-18d0daf7d2a3
md"""
Todos os vertices:
"""

# ╔═╡ 6d9edc1a-63fd-4c2e-a69b-6e8fd007d828
#=╠═╡
plot_vertices(F, p[1] for p in v_sorted)
  ╠═╡ =#

# ╔═╡ e8eb024b-a10c-44be-8f6a-1349270e3cc4
# leste 6, leste 89, norte [362, 364], sul 37, oeste 447, leste 39

# ╔═╡ 11ffd751-41d1-44b4-a8ed-1f533882a1ca
md"""
fui procurando por todos os vertices no poliedro
"""

# ╔═╡ 09598248-1558-4a1b-9673-43847b6d9de1
#=╠═╡
@bind i Slider(1:size(v_sorted, 1), show_value=true)
  ╠═╡ =#

# ╔═╡ dbe8ea2f-6f3c-496c-9c52-1a464ce92a55
#=╠═╡
plot_vertices(F, [v_sorted[i][1]])
  ╠═╡ =#

# ╔═╡ 0c9e7151-4906-4e5e-a274-7f3d2a36dd07
md"""
Vamos pegar 4 pontos de partida
"""

# ╔═╡ 0372b717-c18c-4e30-b42e-839c0811e2c4
#=╠═╡
begin
	plt_sul, traj_sul = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[94108]), passos=300)
	plt_sul, traj_sul = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_sul, passos=300)
end
  ╠═╡ =#

# ╔═╡ bf396529-f045-4b00-aac7-9fdc1f02047f
#=╠═╡
plt_sul
  ╠═╡ =#

# ╔═╡ fe099f74-a943-4fa9-8b10-db54714855a9
#=╠═╡
begin
	plt_norte, traj_norte = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[31824]), passos=300)
	plt_norte, traj_norte = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_norte, passos=300)
end
  ╠═╡ =#

# ╔═╡ 3694f734-e34a-4e6f-b6dc-abc3c73bee6c
#=╠═╡
plt_norte
  ╠═╡ =#

# ╔═╡ 8d7b357e-451a-4098-985d-89c884092fed
#=╠═╡
begin
	_, traj_leste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[107852]), passos=300)
	plt_leste, traj_leste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_leste, passos=300)
end
  ╠═╡ =#

# ╔═╡ 4bdd0dc6-6f5b-4be2-b41e-aa8f73839bf2
#=╠═╡
plt_leste
  ╠═╡ =#

# ╔═╡ 18c35444-7fbb-495f-823c-f7e3d3c252be
#=╠═╡
begin
	_, traj_oeste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, 8.0,reverse(v_sorted[9536]), passos=300)
	plt_oeste, traj_oeste = pipe_trajectory(F, A_exp, B_exp, E_exp, G, d, -8.0,traj_oeste, passos=300)
end
  ╠═╡ =#

# ╔═╡ ae317296-62b3-40ac-99cb-fc6cba3cb3f7
#=╠═╡
plt_oeste
  ╠═╡ =#

# ╔═╡ 67aeadc8-a22b-44f5-a45c-7782c59342ab
md"""
Todas trajetórias observadas no mesmo plot
"""

# ╔═╡ 0b61751d-646d-40f1-99d3-d61aabeea274
#=╠═╡
plot_trajectories(F, [traj_norte, traj_sul, traj_leste, traj_oeste])
  ╠═╡ =#

# ╔═╡ 2ee46796-0afc-4b05-854c-abb5b32c2f57
md"""
Vizualizando demais estados
"""

# ╔═╡ f490a3ea-180b-41a8-a850-2d7052494742
#=╠═╡
plot_expanded_state_trajectory([traj_norte, traj_sul, traj_leste, traj_oeste])
  ╠═╡ =#

# ╔═╡ b0674ad3-758c-4aad-91b2-dc0147691697
#=╠═╡
plot_v_state_trajectory([traj_norte, traj_sul, traj_leste, traj_oeste])
  ╠═╡ =#

# ╔═╡ f29a1973-7a18-45b4-8946-bdf36afb4592
md"""
Passando as matrizes para um arquivo .mat
"""

# ╔═╡ 6bf29487-e6af-41c0-9899-42593ce22649
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ eb96fe20-7b7a-4004-9e27-73eb4293e873
#=╠═╡
plot_along_k([u_leste, u_oeste, u_norte, u_sul])
  ╠═╡ =#

# ╔═╡ b7e3dcbf-b17a-4ba3-8e43-7f9e1e4327db
#=╠═╡
maximum(abs.(u_leste))
  ╠═╡ =#

# ╔═╡ 0d119224-7aee-44c3-a2c0-2a01f11d3a46
#=╠═╡
maximum(abs.(u_oeste))
  ╠═╡ =#

# ╔═╡ 5fc556df-9203-4a8b-ab5d-e6c683f574bd
#=╠═╡
maximum(abs.(u_norte))
  ╠═╡ =#

# ╔═╡ 30129ad8-a6d1-45a2-821c-48a5fdb4e3ea
#=╠═╡
maximum(abs.(u_sul))
  ╠═╡ =#

# ╔═╡ dc0468c1-24bd-48b0-ac1e-aa0c4825fa86
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ d7d2ec55-8de5-40f5-be93-18e1ffc1ef31
#=╠═╡
using Pkg
  ╠═╡ =#

# ╔═╡ 4ac93b77-7567-430f-82fe-d2e2ed75a893
# ╠═╡ disabled = true
#=╠═╡
init_cond_F, init_cond_w = admissable_initCond(A_exp, B_exp*G, F, d, w; fixed_d=true)
  ╠═╡ =#

# ╔═╡ f9a4c39e-7ef6-11f1-ba82-979f074de1d9
# ╠═╡ disabled = true
#=╠═╡
using Pkg
  ╠═╡ =#

# ╔═╡ 5ba9d134-301e-4a52-b39f-950df23687d5
#=╠═╡
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
  ╠═╡ =#

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
