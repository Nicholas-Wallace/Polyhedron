function plot_trajectory(ip::InvariantPolyhedron, traj::Matrix{Float64})
  plt = plot_poly(ip)

  plot_trajectory!(plt, traj)

  return plt
end


function plot_trajectory!(plt, traj::Matrix{Float64})
  plot!(plt, traj[1, :], traj[2, :],
        seriestype = :path,          # connects points with lines
        linewidth = 2,
        marker = :circle,
        markersize = 5,
        markercolor = :red,
        label = "trajectory",
        xlabel = "x1",
        ylabel = "x2")
	
  scatter!(plt, [traj[1, 1]], [traj[2, 1]], 
        marker = :x, 
        markersize = 7, 
        markercolor = :green, 
        label = "Start")

  # Add a red marker at the end
  scatter!(plt, [traj[1, end]], [traj[2, end]], 
        marker = :x, 
        markersize = 7, 
        markercolor = :black, 
        label = "End")
end
