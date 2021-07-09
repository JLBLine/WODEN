import aplpy
import matplotlib.pyplot as plt


##Thanks to the legend here: https://github.com/aplpy/aplpy/issues/423
##that solved the stupid 4 vs 2 axis issue in aplpy
def fix_aplpy_fits(aplpy_obj, dropaxis=2):
    """This removes the degenerated dimensions in APLpy 2.X...
    The input must be the object returned by aplpy.FITSFigure().
    `dropaxis` is the index where to start dropping the axis (by default it assumes the 3rd,4th place).
    """
    temp_wcs = aplpy_obj._wcs.dropaxis(dropaxis)
    temp_wcs = temp_wcs.dropaxis(dropaxis)
    aplpy_obj._wcs = temp_wcs

def do_plot(ax,label,vmin,vmax, colour_bar_label=False, hide_dec_labels=False):
    """Do an aplpy plot with some nice looking things"""
    ##Means ax.recenter will work
    fix_aplpy_fits(ax)

    ax.show_colorscale(cmap='gnuplot',vmin=vmin,vmax=vmax)  #cmap='Blues_r',
    ax.add_colorbar()
    ax.colorbar.show()

    if colour_bar_label:
        ax.colorbar.set_axis_label_text('Jy/beam')

    ax.tick_labels.set_xformat('hh:mm')
    ax.tick_labels.set_yformat('dd:mm')

    ax.ticks.set_xspacing(1.0)

    if hide_dec_labels:
        ax.axis_labels.hide_y()
        ax.tick_labels.hide_y()

    ax.recenter(60.0,-40.0,radius=2.5)

    ax.add_grid()
    ax.grid.set_alpha(0.3)
    ax.grid.show()

    # ax.axis_labels.hide()
    ax.set_title(label)
    # ax.axis_labels.set_ypad(-1)
    # ax.set_nan_color('grey')

    ax.ticks.hide_y()
    ax.ticks.hide_x()

    # ax.ticks.set_color('white')

    ax.ticks.set_minor_frequency(3)

fig = plt.figure(figsize=(12,4))

plot_width = 0.3
plot_height = 0.9
plot_bottom = 0.05

edge_pad = 0.05

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../test_installation/grid_component_models/images/grid_point-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)
ax2 = aplpy.FITSFigure('../../../test_installation/grid_component_models/images/grid_gauss-image.fits',
                        subplot=[plot_width + edge_pad,plot_bottom,plot_width, plot_height],figure=fig)
ax3 = aplpy.FITSFigure('../../../test_installation/grid_component_models/images/grid_shapelet-image.fits',
                        subplot=[(plot_width + edge_pad)*2,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'Point Sources', 0.0, 0.7)
do_plot(ax2,'Gaussian Sources',0.0, 0.28,hide_dec_labels=True)
do_plot(ax3,'Shapelet Sources',0.0, 0.28, hide_dec_labels=True, colour_bar_label=True)

fig.savefig('grid_component_plots.png',bbox_inches='tight', dpi=300)
plt.close()
