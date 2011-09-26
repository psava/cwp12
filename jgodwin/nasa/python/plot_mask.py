from enthought.mayavi import mlab
import numpy
import sys,m8r_numpy as rsf

par = rsf.Par()
slices = par.bool('slices',False)

input = rsf.Input()
n1,n2,n3 = input.n1,input.n2,input.n3
data = input.read()
input.close()

mlab.figure()
plot_data = mlab.pipeline.scalar_field(data)
if slices:
    mlab.pipeline.image_plane_widget(plot_data,
        plane_orientation='x_axes',
        slice_index=10)
    mlab.pipeline.image_plane_widget(plot_data,
        plane_orientation='y_axes',
        slice_index=10)
    mlab.pipeline.image_plane_widget(plot_data,
        plane_orientation='z_axes',
        slice_index=10)
else:
    mlab.pipeline.volume(plot_data,vmin=0,vmax=0.9)
mlab.show()
