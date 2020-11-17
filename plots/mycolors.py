# Definitions for my custom color palettes

import matplotlib.colors
import seaborn
import colorsys


def adjust_lightness(color, amount=0.5):
    # lower numbers darker, higher number brighter
    try:
        c = matplotlib.colors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*matplotlib.colors.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


cmap = matplotlib.colors.ListedColormap(seaborn.color_palette('muted', 10))
cmap_colors = []

for i in range(cmap.N):
    rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
    cmap_colors.append(matplotlib.colors.rgb2hex(rgb))

# Basic colors, not linked to simulation runs
# Just more pleasant to look at (to me) than basic matplotlib colors
colors = {'blue': cmap_colors[0],  # blue
          'orange': cmap_colors[1],  # orange
          'green': adjust_lightness(cmap_colors[2], 0.6),  # green
          'red': cmap_colors[3],  # red
          'purple': cmap_colors[4],  # purple
          'brown': cmap_colors[5],  # brown
          'pink': cmap_colors[6],  # pink
          'gray': cmap_colors[7],  # gray
          'yellow': cmap_colors[8],  # yellow
          'turquoise': adjust_lightness(cmap_colors[9], 0.7),  # turquoise
          }

cmap = matplotlib.colors.ListedColormap(seaborn.color_palette('deep', 10))
cmap_colors = []

for i in range(cmap.N):
    rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
    cmap_colors.append(matplotlib.colors.rgb2hex(rgb))

# Colors for each simulation run, to keep things consistent accross plots
runcolors = {0: cmap_colors[0],  # blue
             1: cmap_colors[1],  # orange
             2: cmap_colors[2],  # green
             3: cmap_colors[3],  # red
             4: cmap_colors[4],  # purple
             5: cmap_colors[5],  # brown
             6: cmap_colors[6],  # pink
             7: cmap_colors[7],  # gray
             8: cmap_colors[8],  # yellow
             9: cmap_colors[9],  # turquoise
          }