# Definitions for my custom legends
# The handlers got changed by hand depending on the different plots
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.colors
import seaborn
import colorsys

cmap = matplotlib.colors.ListedColormap(seaborn.color_palette('muted', 5))
cmap_colors = []

for i in range(cmap.N):
    rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
    cmap_colors.append(matplotlib.colors.rgb2hex(rgb))

print cmap_colors
print len(cmap_colors)


def adjust_lightness(color, amount=0.5):
    try:
        c = matplotlib.colors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*matplotlib.colors.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


colors = {'0': cmap_colors[0],#'#E24A33',  # red
          '1': cmap_colors[1],  # blue
          '2': cmap_colors[2],#'#988ED5',  # purple
          '3': cmap_colors[3],#'#777777',   # gray
          '4': cmap_colors[4],#'#FBC15E',  # yellow
          }  # '#8EBA42'}   # green
# '#FFB5B8'  # pink
#          '5': adjust_lightness(cmap_colors[8], amount=0.95),#'#FBC15E',  # yellow
