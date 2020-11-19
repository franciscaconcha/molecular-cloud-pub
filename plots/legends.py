# Definitions for my custom legends and other stuff
# The handlers got changed by hand depending on the different plots
import matplotlib.lines as mlines
import matplotlib.patches as patches

from mycolors import *

# Number of stars created in each run
Ns = {0: 6289,
      1: 6315,
      2: 6162,
      3: 6283,
      4: 5691,
      5: 6137}

# t_end for star formation
t_end = {0: 4.246,
         1: 4.288,
         2: 10.006,
         3: 5.928,
         4: 8.870,
         5: 7.735}

# t_ini for evolution without star formation
t_ini = {0: 4.250,
         1: 4.290,
         2: 10.010,
         3: 5.930,
         4: 8.875,
         5: 7.740}

class SolidObject(object):
    pass


class DashedObject(object):
    pass


class GasSolidDashedObject(object):
    pass


class DustSolidDashedObject(object):
    pass


class SolidObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0 + 1, y0 + width + 8],
                           [0.3 * height, 0.3 * height],
                           lw=3, ls="-",
                           color='k')
        handlebox.add_artist(l1)
        return [l1]


class DashedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0 + 1, y0 + width + 8],
                           [0.3 * height, 0.3 * height],
                           lw=3, ls="--",
                           color='k')
        handlebox.add_artist(l1)
        return [l1]


class GasSolidDashedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0 + 1, y0 + width + 8],
                           [0.3 * height, 0.3 * height],
                           lw=3, ls="-",
                           color=colors['turquoise'])
        l2 = mlines.Line2D([x0 + 1, y0 + width + 8],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['turquoise'])
        r1 = patches.Rectangle(
            (x0 - 1, y0 + width - 48),  # (x,y)
            y0 + width + 11,  # width
            1.5 * height,  # height
            fill=colors['turquoise'],
            facecolor=colors['turquoise'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        handlebox.add_artist(r1)
        return [l1, l2, r1]


class DustSolidDashedObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        l1 = mlines.Line2D([x0+1, y0 + width + 8],
                           [0.3 * height, 0.3 * height],
                           lw=3, ls="-",
                           color=colors['navy'])
        l2 = mlines.Line2D([x0+1, y0 + width + 8],
                           [0.6 * height, 0.6 * height],
                           lw=3, ls="--",
                           color=colors['navy'])
        r1 = patches.Rectangle(
            (x0 - 1, y0 + width - 47),  # (x,y)
            y0 + width + 11,  # width
            1.5 * height,  # height
            fill=colors['navy'],
            facecolor=colors['navy'],
            # edgecolor="black",
            alpha=0.2,
            # hatch="/",
        )
        handlebox.add_artist(l1)
        handlebox.add_artist(l2)
        handlebox.add_artist(r1)
        return [l1, l2, r1]