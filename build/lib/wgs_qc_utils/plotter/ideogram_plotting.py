from matplotlib.collections import BrokenBarHCollection
import numpy as np


def plot(ideogram, axis):
    xranges = ideogram[['start', 'width']].values
    colors = ideogram["color"].values
    collection = BrokenBarHCollection(xranges, (0, 1), facecolors=colors)
    axis.add_collection(collection)
    axis.set_xlim(0, ideogram.start.max())
    axis.get_yaxis().set_visible(False)
    axis.set_xticks(np.arange(0, ideogram.start.max(), 25))

    return axis
