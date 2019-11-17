
function animate_footprint(pos_end)
vis = Visualizer()
open(vis)

#Plot Mars in MeshCat
image = PngImage(joinpath(MeshCat.VIEWER_ROOT, "..", "data", "8k_mars.png"))
texture = Texture(image=image)
material = MeshLambertMaterial(map=texture)
planet = HyperSphere(Point(0.,0,0), 10.0)
geometry = planet
setobject!(vis["planet"], geometry, material)
settransform!(vis["planet"], LinearMap(AngleAxis(pi/2, 1, 0, 0))) #rotate Planet

#Points (Footprint)
red_material = MeshPhongMaterial(color=RGBA(1, 0, 0, 1.0))
green_material = MeshPhongMaterial(color=RGBA(0, 1, 0, 1.0))


for i in 1:1:36
    if i != 10
        point = HyperSphere(Point(pos_end[i, 1]*10, pos_end[i, 2]*10, pos_end[i, 3]*10), 0.005)
        setobject!(vis["point_$i"], point, red_material)
    end
end

#Departure point
pt_start = HyperSphere(Point(1.036878595663077*10, 0.0, 0.0), 0.01)
setobject!(vis["pt_start"], pt_start, green_material)
end
