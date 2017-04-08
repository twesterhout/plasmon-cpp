#!/usr/bin/env python2

import sys
import time
import argparse

import numpy
import matplotlib
matplotlib.use('Agg')

sys.path.append('/home/twesterh/Bachelor/project/plasmon-cpp/python-tipsi')
import tipsi 



# Define some config parameters
# =============================

# lattice_constant = 0.24612E-9  # [m]
# onsite_potential = 0.0         #
# hopping_value    = 2.8         # [eV]



def parse_options(argv):
    parser = argparse.ArgumentParser('Generate system')
    parser.add_argument( '--type'
                       , dest='sample_type'
                       , choices=[ 'triangle:zigzag', 'triangle:armchair'
                                 , 'sierpinski:carpet'
                                 , 'square'
                                 , 'periodic' 
                                 ]
                       , required=True )
    parser.add_argument( '--lattice-constant'
                       , dest='lattice_constant'
                       , type=float
                       , required=True )
    parser.add_argument( '--hopping-value'
                       , dest='hopping_value'
                       , type=float
                       , required=True )
    parser.add_argument( '--width'
                       , type=int
                       , dest='width'
                       , required=False )
    parser.add_argument( '--start'
                       , type=int
                       , dest='start_width'
                       , required=False )
    parser.add_argument( '--depth'
                       , type=int
                       , dest='depth'
                       , required=False )

    print argv
    options = parser.parse_args(argv)
    print options
    take_needed = \
        { 'triangle:zigzag'   : lambda x: (x.width,)
        , 'triangle:armchair' : lambda x: (x.width,)
        , 'periodic'          : lambda x: (x.width,)
        , 'sierpinski:carpet' : lambda x: (x.start_width, x.depth)
        , 'square'            : lambda x: (x.width,)
        }
    return (options.sample_type,) \
         + take_needed[options.sample_type](options) \
         + (options.lattice_constant, options.hopping_value)


def save_hamiltonian(sample, hamiltonian_filename):
    print "[*] Saving H to file..."
    H = sample.hamiltonian()
    N, _ = H.shape
    with open(hamiltonian_filename, 'w') as f:
        for i in range(N):
            for x in numpy.nditer(H[i]):
                f.write('({0.real},{0.imag})\t'.format(x))
            f.write('\n')
    print "[+] Done."


def save_coordinates(sample, coordinates_filename):
    print "[*] Saving (x,y,z)'s to file..."
    with open(coordinates_filename, 'w') as f:
        for pos in sample._r:
            f.write('{}\t{}\t{}\n'.format(pos[0], pos[1], pos[2]))
    print "[+] Done."


def sierpinski_carpet( start_width
                     , iteration
                     , lattice_constant
                     , hopping_value ):
    W = start_width * 3**iteration
    H = W

    # First we choose a system size (width, height) in unit cells.
    sample = tipsi.square_sheet( W, H
                               , pbc=False
                               , latconst=lattice_constant)
    # Add fractal holes.
    deletesites = []
    for i in xrange(iteration):
        scale = W / (3.**i)

        def in_hole(tag):
            x, y, _, _ = tag
            x = x % scale
            y = y % scale
            return x >= scale / 3. and x < 2. * scale / 3. \
               and y >= scale / 3. and y < 2. * scale / 3.

        for tag in sample.sites:
            if in_hole(tag):
                deletesites.append(tag)

    for tag in deletesites:
        sample.delete(tag)

    # Add hoppings and finalize
    sample.finalize_sites()
    sample.neighbor_hopping(-hopping_value)
    sample.plot()
    sample = sample.finalize()

    # Save results
    save_hamiltonian(sample, "Hamiltonian"
                             + "." + str(start_width) 
                             + "." + str(iteration)
                             + ".dat")
    save_coordinates(sample, "Coordinates"
                             + "." + str(start_width) 
                             + "." + str(iteration)
                             + ".dat")

def square( width
          , lattice_constant
          , hopping_value ):
    # First we choose a system size (width, height) in unit cells.
    sample = tipsi.square_sheet( width, width
                               , pbc=False
                               , latconst=lattice_constant)
    # Add hoppings and finalize
    sample.finalize_sites()
    sample.neighbor_hopping(-hopping_value)
    sample.plot()
    sample = sample.finalize()

    # Save results
    save_hamiltonian(sample, "Hamiltonian"
                             + "." + str(width) 
                             + ".dat")
    save_coordinates(sample, "Coordinates"
                             + "." + str(width)
                             + ".dat")


def triangle_zigzag( width
                   , lattice_constant
                   , hopping_value ):
    sample = tipsi.sample(tipsi.honeycomb_2d_lattice(lattice_constant))
    
    # add sites
    for y in xrange(width+2):
        for x in xrange(width+2-y):
            if (y==0):
                if ((x!=0) and (x!=width+1)):
                    sample.set((x,y,0,1),onsite_potential)
            elif (y==width+1):
                sample.set((x,y,0,0),onsite_potential)
            else:
                sample.set((x,y,0,0),onsite_potential)  
                sample.set((x,y,0,1),onsite_potential)

    # Add hoppings and finalize
    sample.finalize_sites()
    sample.neighbor_hopping(-hopping_value)
    sample = sample.finalize()

    # Save results
    save_hamiltonian(sample, "Hamiltonian"
                             + "." + str(width) 
                             + ".dat")
    save_coordinates(sample, "Coordinates"
                             + "." + str(width)
                             + ".dat")


def triangle_armchair( width
                     , lattice_constant
                     , hopping_value ):
    sample = tipsi.sample(tipsi.honeycomb_2d_lattice(lattice_constant))
    
    # Add sites
    for y in xrange(-width,width+1):
        x_min = 0
        x_max = int(ceil(width*1.5))
        if (width%2==0):
            x_max = x_max-int(floor(abs(y)/2.))
        elif (width%2==1):
            x_max = x_max-int(ceil(abs(y)/2.))
        if (y!=0):
            x_min = abs(y)-1
        if (y<0):
            x_min = x_min-y
            x_max = x_max-y
        for x in xrange(x_min,x_max):
            if ((y!=0) and x==x_min):
                if (y<0):
                    sample.set((x,y,0,1),onsite_potential)
                else:
                    sample.set((x,y,0,0),onsite_potential)
            else:
                sample.set((x,y,0,0),onsite_potential)
                sample.set((x,y,0,1),onsite_potential)
    
    # Add hoppings and finalize
    sample.finalize_sites()
    sample.neighbor_hopping(-hopping_value)
    sample = sample.finalize()

    # Save results
    save_hamiltonian(sample, "Hamiltonian"
                             + "." + str(width) 
                             + ".dat")
    save_coordinates(sample, "Coordinates" +
                             + "." + str(width)
                             + ".dat")



def periodic( width
            , lattice_constant
            , hopping_value ):
    sample = tipsi.honeycomb_sheet( width, width
                                  , pbc=True
                                  , latconst=lattice_constant )
    sample.finalize_sites()
    sample.neighbor_hopping(-hopping_value)
    sample = sample.finalize()

    # Save results
    save_hamiltonian(sample, "Hamiltonian"
                             + "." + str(width) 
                             + ".dat")
    save_coordinates(sample, "Coordinates" +
                             + "." + str(width)
                             + ".dat")


def main():
    options = parse_options(sys.argv[1:])
    constructor = { 'triangle:zigzag'   : triangle_zigzag
                  , 'triangle:armchair' : triangle_armchair
                  , 'periodic'          : periodic
                  , 'sierpinski:carpet' : sierpinski_carpet
                  , 'square'            : square
                  }

    sample_type = options[0]
    arguments   = options[1:]

    print "[*] Building sample..."
    constructor[sample_type](*arguments)
    print "[+] Done."


if __name__ == '__main__':
    main()
