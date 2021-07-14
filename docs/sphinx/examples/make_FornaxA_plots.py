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

def do_plot(ax,label,vmin,vmax, colour_bar_label=False, hide_dec_labels=False,
                                radius=10):
    """Do an aplpy plot with some nice looking things"""
    ##Means ax.recenter will work
    fix_aplpy_fits(ax)

    ax.show_colorscale(cmap='viridis',vmin=vmin,vmax=vmax)  #cmap='Blues_r',
    ax.add_colorbar()
    ax.colorbar.show()

    if colour_bar_label:
        ax.colorbar.set_axis_label_text('Jy/beam')

    ax.tick_labels.set_xformat('hh:mm')
    ax.tick_labels.set_yformat('dd:mm')

    if hide_dec_labels:
        ax.axis_labels.hide_y()
        ax.tick_labels.hide_y()

    ax.recenter(50.67, -37.2, radius=radius)

    ax.set_title(label)

    ax.ticks.hide_y()
    ax.ticks.hide_x()

fig = plt.figure(figsize=(6,4))

plot_width = 0.9
plot_height = 0.9
plot_bottom = 0.05

edge_pad = 0.08

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../examples/FornaxA/images/FornaxA_msclean-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'MSClean Fornax A MWA phase II extended', -0.02, 1.6, radius=0.5, colour_bar_label=True)

fig.savefig('FornaxA_msclean-image.png',bbox_inches='tight', dpi=100)
plt.close()




fig = plt.figure(figsize=(6,4))

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../examples/FornaxA/images/FornaxA_shapelets-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'Shapelet Fornax A MWA phase II extended', -0.02, 1.6, radius=0.5, colour_bar_label=True)

fig.savefig('FornaxA_shapelets-image.png',bbox_inches='tight', dpi=100)
plt.close()
