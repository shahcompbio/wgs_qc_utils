
class EmptyReader():
    pass

def empty_plot(axis):
    axis = axis.text(0.1, 0.1, 'No Data', 
        horizontalalignment='center', 
        verticalalignment='center', 
        transform=axis.transAxes
    )
    return axis