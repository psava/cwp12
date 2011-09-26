from enthought.mayavi import mlab
import numpy
import sys,m8r_numpy as rsf

par = rsf.Par()

n4 = par.int('n4',1)
input = rsf.Input()
input2 = rsf.Input('vel')
n1,n2,n3 = input.n1,input.n2,input.n3

for i in range(n4):
    data = input.read((n3,n2,n1))

vel  = input2.read((n3,n2,n1))


mlab.figure()
wfld_data = mlab.pipeline.scalar_field(data)
vel_data  = mlab.pipeline.scalar_field(vel)

#print plot_data, dir(plot_data)

#data2 = input.read((n3,n2,n1))
#data2 = input.read((n3,n2,n1))
#data2 = input.read((n3,n2,n1))
#plot_data.mlab_source.scalars = data2

mlab.pipeline.image_plane_widget(wfld_data,
    plane_orientation='x_axes',
    slice_index=10,colormap='gray')

mlab.pipeline.image_plane_widget(wfld_data,
    plane_orientation='y_axes',
    slice_index=10,colormap='gray')

mlab.pipeline.image_plane_widget(wfld_data,
    plane_orientation='z_axes',
    slice_index=10,colormap='gray')

#mlab.pipeline.volume(vel_data,vmin=0,vmax=0.7)
#mlab.contour3d(vel)
#@mlab.show
#@mlab.animate
#def animate():
#    f = mlab.gcf()
#    while 1:
#        #f.scene.camera.azimuth(10)
#        plot_data.mlab_source.scalars = input.read((n3,n2,n1))
#        #f.scene.render()
#        yield
#
##mlab.show()
#animate()
input.close()
input2.close()
mlab.show()

