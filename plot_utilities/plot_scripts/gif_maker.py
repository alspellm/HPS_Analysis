import glob
import os
from PIL import Image

def make_gif(frame_folder,name):
    frames = [Image.open(image) for image in sorted(glob.glob('%s.png'%(frame_folder)),key=lambda name: int(os.path.basename(name).replace(".png",'').split("_")[5]))]
    frame_one = frames[0]
    frame_one.save("%s.gif"%(name), format="GIF", append_images=frames,
               save_all=True, duration=1000, loop=0)

if __name__ == "__main__":
    make_gif("/home/alic/HPS/projects/simps/zbi_opt/zalpha_opt/test_plots/tritrig*",'background_models')
