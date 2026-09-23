struct LTISystem
  A::Matrix{Float64}
  B::Matrix{Float64}
  C::Matrix{Float64}
end

struct SystemDelayed
  A::Matrix{Float64}
  Ad::Matrix{Float64}
  d::Int
end
