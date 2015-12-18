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
  from numpy.random import random

  t = 2.*pi*random(n)
  u = random(n)+random(n)
  r = zeros(n,'float')
  mask = u>1.
  xmask = logical_not(mask)
  r[mask] = 2.-u[mask]
  r[xmask] = u[xmask]
  xyp = column_stack( (rr*r*cos(t),rr*r*sin(t)) )
  dartsxy  = xyp + array([xx,yy])

  return dartsxy

def darts(n, xx, yy, rr, dst):
  """
  get at most n random, uniformly distributed, points in a circle.
  centered at (xx,yy), with radius rr. points are no closer to each other
  than dst.
  """

  from scipy.spatial import distance
  from numpy import array
  cdist = distance.cdist

  dartsxy = random_points_in_circle(n, xx, yy, rr)

  jj = []

  ## remove new nodes that are too close to other
  ## new nodes
  dists = cdist(dartsxy,dartsxy,'euclidean')
  for j in xrange(n-1):
    if all( dists[j,j+1:] > dst ):
      jj.append(j)

  res = dartsxy[array(jj,'int'),:]

  return res

