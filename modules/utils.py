#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function


def random_points_in_circle(n, x, y, rad):
  """
  get n random points in a circle.
  """

  from numpy import zeros, logical_not
  from numpy import column_stack
  from numpy import cos
  from numpy import sin
  from numpy import array
  from numpy import pi
  from numpy import reshape
  from numpy.random import random

  rnd = random(size=(n,3))
  t = 2.*pi*rnd[:,0]
  u = rnd[:,1:].sum(axis=1)
  r = zeros(n,'float')
  mask = u>1.
  xmask = logical_not(mask)
  r[mask] = 2.-u[mask]
  r[xmask] = u[xmask]
  xyp = reshape(rad*r,(n,1))*column_stack( (cos(t),sin(t)) )
  dartsxy  = xyp + array([x,y])
  return dartsxy


def darts(n, x, y, rad, dst, old_darts=None):
  """
  get at most n random, uniformly distributed, points in a circle.
  centered at (x,y), with radius rr. points are no closer to each other
  than dst.
  """

  from numpy import array
  from numpy import zeros
  from numpy import row_stack 
  from scipy.spatial import cKDTree as kdt

  dartsxy = random_points_in_circle(n, x, y, rad)
  jj = zeros(n,'bool')

  tree = kdt(dartsxy)
  for j,near in enumerate(tree.query_ball_point(dartsxy, dst)):
    if len(near)<2:
      jj[j] = True

  res = dartsxy[jj,:]

  if old_darts is not None:
    return row_stack([old_darts, res])

  return res


def spatial_sort(paths, init_rad=0.01):

  from numpy import row_stack
  from numpy import array
  from numpy import zeros
  from numpy.linalg import norm
  from scipy.spatial import cKDTree as kdt

  num = len(paths)

  res = []

  unsorted = set(range(2*num))

  xs = zeros((2*num,2), 'float')
  x_path = zeros(2*num, 'int')

  for i, path in enumerate(paths):
    xs[i,:] = path[0,:]
    xs[num+i,:] = path[-1,:]

    x_path[i] = i
    x_path[num+i] = i

  tree = kdt(xs)

  count = 0
  pos = array([0,0],'float')

  while count<num:

    rad = init_rad
    while True:

      near = tree.query_ball_point(pos, rad)
      cands = list(set(near).intersection(unsorted))
      if not cands:
        rad *= 2.0
        continue

      dst = norm(pos - xs[cands,:], axis=1)
      cp = dst.argmin()
      uns = cands[cp]
      break

    path_ind = x_path[uns]
    path = paths[path_ind]

    if uns>=num:
      res.append(path[::-1])
      pos = paths[path_ind][0,:]
      unsorted.remove(uns)
      unsorted.remove(uns-num)

    else:
      res.append(path)
      pos = paths[path_ind][-1,:]
      unsorted.remove(uns)
      unsorted.remove(uns+num)

    count += 1

  return res


def export_svg(fn, paths, size):

  from cairo import SVGSurface, Context
  from numpy import array

  one = 1.0/size
  s = SVGSurface(fn, size, size)
  c = Context(s)

  c.set_line_width(0.1)

  paths = spatial_sort(paths)

  for path in paths: 
    path *= size

    c.new_path()
    c.move_to(*path[0,:])
    for p in path[1:]:
      c.line_to(*p)
    c.stroke()

