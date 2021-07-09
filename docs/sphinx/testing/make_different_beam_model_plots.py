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
                                radius=10, grid=True):
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



    if hide_dec_labels:
        ax.axis_labels.hide_y()
        ax.tick_labels.hide_y()

    ax.recenter(60.0,-40.0,radius=radius)

    if grid:

        ax.ticks.set_xspacing(5.0)
        ax.ticks.set_yspacing(5.0)

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

    # ax.ticks.set_minor_frequency(3)

fig = plt.figure(figsize=(12,4))

plot_width = 0.3
plot_height = 0.9
plot_bottom = 0.05

edge_pad = 0.05

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../test_installation/different_beam_models/images/multi-comp_grid_None-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)
ax2 = aplpy.FITSFigure('../../../test_installation/different_beam_models/images/multi-comp_grid_Gaussian-image.fits',
                        subplot=[plot_width + edge_pad,plot_bottom,plot_width, plot_height],figure=fig)
ax3 = aplpy.FITSFigure('../../../test_installation/different_beam_models/images/multi-comp_grid_MWA_FEE-image.fits',
                        subplot=[(plot_width + edge_pad)*2,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'No Beam', 0.0, 0.7)
do_plot(ax2,'Gaussian Beam \nFWHM 20.0 deg @ 150MHz',0.0, 0.7,hide_dec_labels=True)
do_plot(ax3,'MWA FEE Beam',0.0, 0.7, hide_dec_labels=True, colour_bar_label=True)

fig.savefig('different_beam_plots.png',bbox_inches='tight', dpi=300)
plt.close()





fig = plt.figure(figsize=(8,4))

plot_width = 0.45
plot_height = 0.9
plot_bottom = 0.05

edge_pad = 0.08

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../test_installation/different_beam_models/images/multi-comp_grid_MWA_FEE-psf.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)
ax2 = aplpy.FITSFigure('../../../test_installation/different_beam_models/images/multi-comp_grid_EDA2-psf.fits',
                        subplot=[plot_width + edge_pad,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'multi-comp_grid_MWA_FEE-psf.fits', 0.0, 0.7, grid=False, radius=2)
do_plot(ax2,'multi-comp_grid_EDA2-psf.fits',0.0, 0.7,hide_dec_labels=True,
             colour_bar_label=True, grid=False, radius=2)

fig.savefig('MWA-vs-EDA2_psf.png',bbox_inches='tight', dpi=300)
plt.close()


fig = plt.figure(figsize=(6,4))

plot_width = 0.9
plot_height = 0.9
plot_bottom = 0.05

edge_pad = 0.08

##Make some aplpy figures
##Manually set the subplots to maximise space usage
ax1 = aplpy.FITSFigure('../../../test_installation/different_beam_models/images/multi-comp_grid_EDA2-image.fits',
                        subplot=[0.01,plot_bottom,plot_width, plot_height],figure=fig)

do_plot(ax1,'EDA2 beam and array layout', -2.0, 16)

fig.savefig('EDA2_layout_image.png',bbox_inches='tight', dpi=300)
plt.close()
