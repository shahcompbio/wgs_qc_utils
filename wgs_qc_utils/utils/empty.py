def empty_plot(axis, label, snv_cn=False):

    x = axis.get_xlim()[1]/2
    y = axis.get_ylim()[1]/2
    text = "{}:  No data to display".format(label)
    size = 20
    if snv_cn:
        text = "{}: \n No data \n to display".format(label)
        size = 10
    axis.text(x,y, text, 
        horizontalalignment='center', 
        verticalalignment='center', 
        transform=axis.transAxes,
        fontsize=size
    )
    axis.xaxis.set_visible(False)
    axis.yaxis.set_visible(False)
    axis.spines["bottom"].set_color("red")
    axis.spines["left"].set_color("red")
    axis.spines["right"].set_color("red")
    axis.spines["top"].set_color("red")

    return axis