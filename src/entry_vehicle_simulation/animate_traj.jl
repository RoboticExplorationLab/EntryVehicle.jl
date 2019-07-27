#File contains function for Visualization of entry Trajectory - Mars Environment

function animate_traj(t_sim, Z)
vis = Visualizer()
open(vis)
delete!(vis)
#Plot Mars in MeshCat
image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "Mars.png"))
texture = Texture(image=image)
material = MeshLambertMaterial(map=texture)
planet = HyperSphere(Point(0.,0,0), 10.0)
geometry = planet
setobject!(vis["planet"], geometry, material)
settransform!(vis["planet"], LinearMap(AngleAxis(pi/2, 1, 0, 0))) #rotate Planet

#Plot Spacecraft
image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "Rex.png"))
texture= Texture(image = image)
material = MeshLambertMaterial(map=texture)
cap = load(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "orion_100_smaller.obj"), GLUVMesh)
setobject!(vis["vehicle"], cap, material)
#settransform!(vis["vehicle"], LinearMap(AngleAxis(pi/2, 1.0, 0, 0)))
#settransform!(vis["vehicle"], LinearMap(Quat(Q...)))
    #Material
red_material = MeshPhongMaterial(color=RGBA(1, 0, 0, 1.0))
green_material = MeshPhongMaterial(color=RGBA(0, 1, 0, 1.0))

    #Points Trajectory
sphere_small = HyperSphere(Point(0.0,0.0,0.0), 0.005)

    #Plot Trajectory
traj = vis["traj"]
vehicle = vis["vehicle"]

N = length(t_sim)


    #Building Animation
anim = MeshCat.Animation()
for i = 1:N
    MeshCat.atframe(anim,vis,i) do frame
        settransform!(frame["vehicle"], compose(Translation(Z[1:3, i].*10...),LinearMap(Quat(qmult(Z[4:7, i], QQ)...))))
    end
    setobject!(vis["traj"]["t$i"],sphere_small, green_material)
    settransform!(vis["traj"]["t$i"], Translation(Z[1, i]*10, Z[2, i]*10, Z[3, i]*10))

end

MeshCat.setanimation!(vis,anim)
settransform!(vis["/Cameras/default"], Translation(10, 0, 0))

end
