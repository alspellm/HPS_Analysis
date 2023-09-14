from array import array
import os,sys,glob
import argparse

parser = argparse.ArgumentParser(description="Flags ")
parser.add_argument('--images', '-i', type=str, dest="images",
                  help="Specify directory where images are stored", metavar="images", default=".")
parser.add_argument('--title', '-t', type=str, dest="title",
                  help="Specify HTML Page title", metavar="title", default="Plots")

options = parser.parse_args()

image_dir = options.images
title = options.title

def makeHTML(image_dir,title):
    os.chdir(image_dir)
    plots = glob.glob('*.png') + glob.glob('*.pdf')
    plots.sort()
    outfilebase = os.path.split(image_dir)[1]
    f = open(outfilebase+'.html',"w+")
    f.write("<!DOCTYPE html\n")
    f.write(" PUBLIC \"-//W3C//DTD HTML 3.2//EN\">\n")
    f.write("<html>\n")
    f.write("<head><title>"+ title +" </title></head>\n")
    f.write("<body bgcolor=\"EEEEEE\">\n")
    f.write("<table border=\"0\" cellspacing=\"5\" width=\"100%\">\n")
    for i in range(0,len(plots)):
        pname = ""
        offset = 1
        if i==0 or i%2==0: f.write("<tr>\n")
        f.write("<td width=\"10%\"><a target=\"_blank\" href=\"" + plots[i] + "\"><img src=\"" + plots[i] + "\" alt=\"" + plots[i] + "\" title=\"" + pname + "\" width=\"85%\" ></a></td>\n")
        if i==offset:
            f.write("</tr>\n")
        elif (i>offset and (i-offset)%2==0) or i==len(plots):
            f.write("</tr>\n")

    f.write("</table>\n")
    f.write("</body>\n")
    f.write("</html>")
    f.close()


#makeHTML('/home/alic/HPS/projects/simp_analysis/scripts/plots/trackhit/trackhit_ana','Run 7800 vs Tritrig_Beam TrackHitAna')
#makeHTML('/home/alic/HPS/projects/simp_analysis/scripts/plots/trackhit/nhit_plots','Run 7800 vs Tritrig_Beam Momentum')
#makeHTML('/home/alic/HPS/projects/simp_analysis/scripts/plots/trackhit/debug_plots','Run 7800 vs Tritrig_Beam Momentum')
print(image_dir)
print(title)
makeHTML('%s'%(image_dir),'%s'%(title))
