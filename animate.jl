_path = "/home/user/Documents/Julia_Code_Remy/animation_no_damping"

using FileIO
using ImageMagick
using Plots
# using Images, ImageView

anim = @animate for _img in readdir(_path)[1:5:end]
    im_ = FileIO.load(joinpath(_path,_img))
    im_crop = im_[200:end-50,200:end-50]

    plot(im_crop,axis=false)
end
gif(anim, joinpath(_path,"entry2.gif"), fps = 25)
