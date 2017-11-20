#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 13:56:12 2017

@author: johannes
"""

import numpy as np

def _read_cube_header(f):
    # Read the title
    title = f.readline().strip()
    # skip the second line
    f.readline()

    def read_grid_line(line):
        """Read a grid line from the cube file"""
        words = line.split()
        return (
            int(words[0]),
            np.array([float(words[1]), float(words[2]), float(words[3])], float)
            # all coordinates in a cube file are in atomic units
        )

    # number of atoms and origin of the grid
    natom, origin = read_grid_line(f.readline())
    print origin
    # numer of grid points in A direction and step vector A, and so on
    shape0, axis0 = read_grid_line(f.readline())
    shape1, axis1 = read_grid_line(f.readline())
    shape2, axis2 = read_grid_line(f.readline())
    shape = np.array([shape0, shape1, shape2], int)
    axes = np.array([axis0, axis1, axis2])

#    cell = Cell(axes*shape.reshape(-1,1))
#    ugrid = UniformGrid(origin, axes, shape, np.ones(3, int))

    def read_coordinate_line(line):
        """Read an atom number and coordinate from the cube file"""
        words = line.split()
        return (
            int(words[0]), float(words[1]),
            np.array([float(words[2]), float(words[3]), float(words[4])], float)
            # all coordinates in a cube file are in atomic units
        )

    numbers = np.zeros(natom, int)
    pseudo_numbers = np.zeros(natom, float)
    coordinates = np.zeros((natom, 3), float)
    for i in range(natom):
        numbers[i], pseudo_numbers[i], coordinates[i] = read_coordinate_line(f.readline())
        # If the pseudo_number field is zero, we assume that no effective core
        # potentials were used.
        if pseudo_numbers[i] == 0.0:
            pseudo_numbers[i] = numbers[i]

    return title, coordinates, numbers, origin, shape, axes, pseudo_numbers


def _read_cube_data(f, shape):
    data = np.zeros(tuple(shape), float)
    tmp = data.ravel()
    counter = 0
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        words = line.split()
        for word in words:
            tmp[counter] = float(word)
            counter += 1
    return data


def load_cube(filename):
    '''Load data from a cube file
       **Arguments:**
       filename
            The name of the cube file
       **Returns** a dictionary with ``title``, ``coordinates``, ``numbers``,
       ``cube_data``, ``grid``, ``pseudo_numbers``.
    '''
    with open(filename) as f:
        title, coordinates, numbers, origin, shape, axes, pseudo_numbers = _read_cube_header(f)
        data = _read_cube_data(f, shape)
        return {
            'title': title,
            'coordinates': coordinates,
            'numbers': numbers,
            'origin': origin,
            'shape': shape,
            'axes': axes,
            'cube_data': data,
            'pseudo_numbers': pseudo_numbers,
        }


def _write_cube_header(f, title, coordinates, numbers, origin, shape, axes, pseudo_numbers):
    print >> f, title
    print >> f, 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
    natom = len(numbers)
    x, y, z = origin
    print >> f, '%5i % 11.6f % 11.6f % 11.6f' % (natom, x, y, z)
    rvecs = axes
    for i in range(3):
        x, y, z = rvecs[i]
        print >> f, '%5i % 11.6f % 11.6f % 11.6f' % (shape[i], x, y, z)
    for i in range(natom):
        q = pseudo_numbers[i]
        x, y, z = coordinates[i]
        print >> f, '%5i % 11.6f % 11.6f % 11.6f % 11.6f' % (numbers[i], q, x, y, z)


def _write_cube_data(f, cube_data):
    counter = 0
    for value in cube_data.flat:
        f.write(' % 12.5E' % value)
        if counter%6 == 5:
            f.write('\n')
        counter += 1


def dump_cube(filename, data):
    '''Write a IOData to a .cube file.
       **Arguments:**
       filename
            The name of the file to be written. This usually the extension
            ".cube".
       data
            An IOData instance. Must contain ``coordinates``, ``numbers``,
            ``grid``, ``cube_data``. May contain ``title``, ``pseudo_numbers``.
    '''
    with open(filename, 'w') as f:
        _write_cube_header(f, data["title"], data["coordinates"], data["numbers"],
                           data["origin"], data["shape"], data["axes"], data["pseudo_numbers"])
        _write_cube_data(f, data["cube_data"])




