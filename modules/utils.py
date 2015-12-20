#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function


def random_points_in_circle(n,xx,yy,rr):
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
  xyp = reshape(rr*r,(n,1))*column_stack( (cos(t),sin(t)) )
  dartsxy  = xyp + array([xx,yy])
  return dartsxy

def darts(n, xx, yy, rr, dst):
  """
  get at most n random, uniformly distributed, points in a circle.
  centered at (xx,yy), with radius rr. points are no closer to each other
  than dst.
  """

  from numpy import array
  from scipy.spatial import cKDTree as kdt

  ## remove new nodes that are too close to other
  ## new nodes

  dartsxy = random_points_in_circle(n, xx, yy, rr)
  tree = kdt(dartsxy)
  near = tree.query_ball_point(dartsxy, dst)
  jj = []
  for j,n in enumerate(near):
    if len(n)<2:
      jj.append(j)

  res = dartsxy[jj,:]
  return res

def export_svg(fn, paths, size):

  from cairo import SVGSurface, Context
  from numpy import array

  one = 1.0/size
  s = SVGSurface(fn, size, size)
  c = Context(s)

  c.set_line_width(0.1)

  for path in paths: 
    path *= size

    c.new_path()
    c.move_to(*path[0,:])
    for p in path[1:]:
      c.line_to(*p)
    c.stroke()

