volume.f90: 
compute volumetric statistics or statistics in xy, xz, yz directions,
including iso-surface of swirl, vorticity, density and pressure gradients
outputs: 
         volume_xxxxxxxx.dat (if ivolume = 1)
         plane_xy_xxxxxxxx_kxxxx.dat (if iplane_xy = 1)
         plane_xz_xxxxxxxx_jxxxx.dat (if iplane_xz = 1)
         plane_yz_xxxxxxxx_ixxxx.dat (if iplane_yz = 1)
         (xxxxxxxx is flowfield index and xxxx is location i/j/k index)
usage: ./volume < volume.inp

