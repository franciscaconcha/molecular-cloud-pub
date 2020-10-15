# Definitions for my custom legends and other stuff
# The handlers got changed by hand depending on the different plots
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.colors
import seaborn
import colorsys

# Making this to make sure label, colors, and line styles are consistent through all the plots
labels = {'R01': 'R = 0.1 pc',
          'R03': 'R = 0.3 pc',
          'R05': 'R = 0.5 pc',
          'R1': 'R = 1 pc',
          'R25': 'R = 2.5 pc',
          'R5': 'R = 5 pc',
          'N1E4R25': 'R = 2.5 pc',
          'N1E4R5': 'R = 5 pc'}

Ns = {0: 6369,
      1: 6273,
      2: 6161,
      3: 6258,
      4: 5717,
      5: 6288,
      6: 5218,
      7: 6081,
      8: 6221,
      9: 5817}

cmap = matplotlib.colors.ListedColormap(seaborn.color_palette('cubehelix', 13))
cmap_colors = []

for i in range(cmap.N):
    rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
    cmap_colors.append(matplotlib.colors.rgb2hex(rgb))

print len(cmap_colors)


def adjust_lightness(color, amount=0.5):
    try:
        c = matplotlib.colors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*matplotlib.colors.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])


colors = {'R01': adjust_lightness(cmap_colors[2], amount=1.5),#'#E24A33',  # red
          'R03': adjust_lightness(cmap_colors[5], amount=0.95),#'#348ABD',  # blue
          'R05': adjust_lightness(cmap_colors[6], amount=1.5),#'#988ED5',  # purple
          'R1': adjust_lightness(cmap_colors[7], amount=1.1),#'#777777',   # gray
          'R25': adjust_lightness(cmap_colors[8], amount=0.95),#'#FBC15E',  # yellow
          'R5': adjust_lightness(cmap_colors[9], amount=0.4),#'#8EBA42'}   # green
          'data1': matplotlib.colors.to_rgba(adjust_lightness(cmap_colors[10], amount=0.4)),  # '#777777',   # gray
          'N1E4R25': adjust_lightness(cmap_colors[12], amount=0.4),  # '#FBC15E',  # yellow
          'data2': matplotlib.colors.to_rgba(adjust_lightness(cmap_colors[11], amount=0.8))}  # '#8EBA42'}   # green
# '#FFB5B8'  # pink
lines = {'R01': '-',
         'R03': '-',
         'R05': '-',
         'R1': '-',
         'R25': '-',
         'R5': '-'}

right_limits = {'R01': 1E7,
                'R03': 1E6,
                'R05': 1E5,
                'R1': 1E4,
                'R25': 1E3,
                'R5': 1E2}

folders = ['N1E3_R01', 'N1E3_R03', 'N1E3_R05', 'N1E3_R1', 'N1E3_R25', 'N1E3_R5']


class EmptyObject(object):
    pass


class DottedObject(object):
    pass


class DashedObject(object):
    pass


class SolidObject(object):
    pass


class SolidObjectLight(object):
    pass


class SolidShadedObject(object):
    pass


class DashedShadedObject(object):
    pass


class DottedShadedObject(object):
    pass


class R01Object(object):
    pass


class R03Object(object):
    pass


class R05Object(object):
    pass


class R1Object(object):
    pass


class R25Object(object):
    pass


class R5Object(object):
    pass


class R01ShadedObject(object):
    pass


class R03ShadedObject(object):
    pass


class R05ShadedObject(object):
    pass


class R1ShadedObject(object):
    pass


class R25ShadedObject(object):
    pass


class R5ShadedObject(object):
    pass


class EmptyObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, x0 + 1],
                           [0.1 * height, 0.1 * height],
                           lw=3, ls="-",
                           color='black')
        #handlebox.add_artist(l1)
        return [l1]


class DottedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls=":",
                           color='black')  # Have to change color by hand for different plots
        handlebox.add_artist(l1)
        return [l1]


class DashedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color='black')  # Have to change color by hand for different plots
        handlebox.add_artist(l1)
        return [l1]


class SolidObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 3],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="-",
                           color='black')  # Have to change color by hand for different plots
        handlebox.add_artist(l1)
        return [l1]


class SolidObjectLightHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 3],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-", alpha=0.4,
                           color='black')
        handlebox.add_artist(l1)
        return [l1]


class SolidShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0+1, y0 + width + 8],
                           [0.4 * height, 0.4 * height],
                           lw=3, ls="-",
                           color='black')
        l2 = patches.Rectangle(
            (x0 - 1, y0 + width - 40),  # (x,y)
            y0 + width + 11,  # width
            1.4 * height,  # height
            fill='black',
            facecolor='black',
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        return [l1, l2]


class DashedShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color='black')  # Have to change color by hand for different plots
        l2 = patches.Rectangle(
            (x0 - 1, y0 + width - 33),  # (x,y)
            y0 + width + 11,  # width
            1.4 * height,  # height
            fill='black',
            facecolor='black',
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        return [l1, l2]


class DottedShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls=":",
                           color='black')  # Have to change color by hand for different plots
        l2 = patches.Rectangle(
            (x0 - 1, y0 + width - 33),  # (x,y)
            y0 + width + 11,  # width
            1.4 * height,  # height
            fill='black',
            facecolor='black',
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        return [l1, l2]


class R01ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R01'])  # Have to change color by hand for different plots
        """l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['R01'])"""
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        #return [l1, l2]
        return [l1]


class R03ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 8],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R03'])  # Have to change color by hand for different plots
        """l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['R03'])"""
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        #return [l1, l2]
        return [l1]


class R05ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 8],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R05'])  # Have to change color by hand for different plots
        """l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['R03'])"""
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        #return [l1, l2]
        return [l1]


class R1ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 8],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R1'])  # Have to change color by hand for different plots
        """l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['R03'])"""
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        #return [l1, l2]
        return [l1]


class R25ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 8],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R25'])  # Have to change color by hand for different plots
        """l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['R25'])"""
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        #return [l1, l2]
        return [l1]


class R5ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R5'])  # Have to change color by hand for different plots
        """l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['R5'])"""
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        #return [l1, l2]
        return [l1]


"""class R01ShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0+1, y0 + width + 8],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R01'])  # Have to change color by hand for different plots
        l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.2 * height, 0.2 * height],
                           lw=3, ls="--",
                           color=colors['R01'])
        l3 = patches.Rectangle(
            (x0, y0 + width - 34),  # (x,y)
            y0 + width + 10,  # width
            1.4 * height,  # height
            fill=colors['R01'],
            facecolor=colors['R01'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        handlebox.add_artist(l3)
        #return [l1, l2, l3]
        return [l1, l3]


class R03ShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls=":",
                           color=colors['R03'])  # Have to change color by hand for different plots
        l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.2 * height, 0.2 * height],
                           lw=3, ls="--",
                           color=colors['R03'])
        l3 = patches.Rectangle(
            (x0, y0 + width - 36),  # (x,y)
            y0 + width + 10,  # width
            1.4 * height,  # height
            fill=colors['R03'],
            facecolor=colors['R03'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        handlebox.add_artist(l3)
        return [l1, l2, l3]
        #return [l2, l3]


class R05ShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls=":",
                           color=colors['R05'])  # Have to change color by hand for different plots
        l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.2 * height, 0.2 * height],
                           lw=3, ls="--",
                           color=colors['R05'])
        l3 = patches.Rectangle(
            (x0, y0 + width - 36),  # (x,y)
            y0 + width + 10,  # width
            1.4 * height,  # height
            fill=colors['R05'],
            facecolor=colors['R05'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        handlebox.add_artist(l3)
        return [l1, l2, l3]
        #return [l1, l3]


class R1ShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls=":",
                           color=colors['R1'])  # Have to change color by hand for different plots
        l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.2 * height, 0.2 * height],
                           lw=3, ls="--",
                           color=colors['R1'])
        l3 = patches.Rectangle(
            (x0, y0 + width - 36),  # (x,y)
            y0 + width + 10,  # width
            1.4 * height,  # height
            fill=colors['R1'],
            facecolor=colors['R1'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        handlebox.add_artist(l3)
        return [l1, l2, l3]
        #return [l1, l3]


class R25ShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0, y0 + width + 10],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls=":",
                           color=colors['R25'])  # Have to change color by hand for different plots
        l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.2 * height, 0.2 * height],
                           lw=3, ls="--",
                           color=colors['R25'])
        l3 = patches.Rectangle(
            (x0, y0 + width - 36),  # (x,y)
            y0 + width + 10,  # width
            1.4 * height,  # height
            fill=colors['R25'],
            facecolor=colors['R25'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        handlebox.add_artist(l3)
        return [l1, l2, l3]
        #return [l2, l3]


class R5ShadedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0+2, y0 + width + 8],
                           [0.5 * height, 0.5 * height],
                           lw=3, ls="-",
                           color=colors['R5'])  # Have to change color by hand for different plots
        l2 = mlines.Line2D([x0, y0 + width + 10],
                           [0.2 * height, 0.2 * height],
                           lw=3, ls="--",
                           color=colors['R5'])
        l3 = patches.Rectangle(
            (x0, y0 + width - 34),  # (x,y)
            y0 + width + 10,  # width
            1.4 * height,  # height
            fill=colors['R5'],
            facecolor=colors['R5'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        #handlebox.add_artist(l2)
        handlebox.add_artist(l3)
        #return [l1, l2, l3]
        return [l1, l3]"""