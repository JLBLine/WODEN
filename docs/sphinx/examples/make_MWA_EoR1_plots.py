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
            ra=60.0, dec=-27.0, radius=13.0):
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

    # ax.ticks.set_xspacing(1.0)

    if hide_dec_labels:
        ax.axis_labels.hide_y()
        ax.tick_labels.hide_y()

    ax.recenter(ra,dec,radius=radius)

    # ax.add_grid()
    # ax.grid.set_alpha(0.3)
    # ax.grid.show()

    # ax.axis_labels.hide()
    ax.set_title(label)
    # ax.axis_labels.set_ypad(-1)
    # ax.set_nan_color('grey')

    ax.ticks.hide_y()
    ax.ticks.hide_x()

    # ax.ticks.set_color('white')

    # ax.ticks.set_minor_frequency(3)

fig = plt.figure(figsize=(8,8))

plot_width = 0.95
plot_height = 0.95
plot_bottom = 0.05

edge_pad = 0.05

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../examples/MWA_EoR1/images/MWA_EoR1_smaller-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'MWA EoR smaller', 0.0, 0.2, colour_bar_label=True)

fig.savefig('MWA_EoR1_plot_smaller.svg',bbox_inches='tight', dpi=300)
plt.close()

fig = plt.figure(figsize=(12,8))

plot_width = 0.4
plot_height = 0.95
plot_bottom = 0.05

edge_pad = 0.15

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../examples/MWA_EoR1/images/MWA_EoR1_large-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)
ax2 = aplpy.FITSFigure('../../../examples/MWA_EoR1/images/MWA_EoR1_large-image.fits',
                        subplot=[plot_width + edge_pad,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'MWA EoR Larger', 0.0, 0.1, radius=40)

do_plot(ax2,'MWA EoR Larger Zoom', 0.0, 0.1, radius=15, dec=0)

fig.savefig('MWA_EoR1_plot_larger.svg',bbox_inches='tight', dpi=300)
plt.close()







# ax2 = aplpy.FITSFigure('../../../examples/MWA_EoR1/images/MWA_EoR1_large-image.fits',
#                         subplot=[plot_width + edge_pad,plot_bottom,plot_width, plot_height],figure=fig)
# ax3 = aplpy.FITSFigure('../../../examples/MWA_EoR1/images/MWA_EoR1_large-image.fits',
#                         subplot=[(plot_width + edge_pad)*2,plot_bottom,plot_width, plot_height],figure=fig)
