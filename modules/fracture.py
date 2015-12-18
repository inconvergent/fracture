# -*- coding: utf-8 -*-

from numpy import pi
from numpy import array

TWOPI = pi*2
HPI = pi*0.5

class Fracture(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      init_dst=0.0, 
      crack_dot=0.95,
      crack_dst=0.05
      ):

    self.init_num = init_num
    self.init_rad = init_rad
    self.init_dst = init_dst
    self.crack_dot = crack_dot
    self.crack_dst = crack_dst

    self.fractures = []
    self.old_fractures = []

    self.hit = set()

    self.__make_sources()

    self.make_fracture(x=array([0.5,0.5]), dx=array([0.0,1.0]))
    self.make_fracture(x=array([0.5,0.5]), dx=array([0.0,-1.0]))

    self.make_fracture(x=array([0.5,0.5]), dx=array([1.0,0.0]))
    self.make_fracture(x=array([0.5,0.5]), dx=array([-1.0,0.0]))

  def __make_sources(self):

    from scipy.spatial import cKDTree as kdt
    from utils import darts

    sources = darts(self.init_num, 0.5, 0.5, self.init_rad, self.init_dst)
    tree = kdt(sources)

    self.sources = sources
    self.tree = tree

  def make_fracture(self, x, dx):

    _,p = self.tree.query(x,1)
    self.fractures.append([(p, dx)])

  def make_fracture_from_old(self):

    from numpy.random import randint
    from numpy.random import random
    from numpy import cos
    from numpy import sin

    cands = array(list(self.hit))
    i = cands[randint(len(cands))]
    x = self.sources[i,:]
    a = random()*TWOPI
    dx = array([cos(a), sin(a)])
    self.make_fracture(x=x, dx=dx)

  def fracture(self):

    from numpy import arctan2
    from numpy import logical_and
    from numpy import arange
    from numpy import reshape
    from numpy import abs
    from numpy.linalg import norm

    query = self.tree.query_ball_point
    sources = self.sources
    crack_dst = self.crack_dst
    dl = self.crack_dot

    keep = set()

    for i in range(len(self.fractures)):

      try:

        p,dx = self.fractures[i][-1]
        dx = dx.reshape((1,2))
        px = sources[p,:]

        near = query(px, crack_dst)
        diff = sources[near,:] - px
        nrm = norm(diff,axis=1).reshape((-1,1))

        zeromask = nrm<=0
        nrm[zeromask] = 1000.0

        diff /= nrm

        dot = (dx * diff).sum(axis=1)
        mask = dot>dl

        masked_arange = arange(len(mask))[mask]
        masked_sources = sources[mask,:]
        masked_diff = diff[mask]
        masked_nrm = nrm[mask]

        mp = masked_nrm.argmin()
        m = masked_nrm[mp]
        ind = masked_arange[mp]

        if m<=0.0 or m>1.0:
          raise ValueError

      except ValueError:

        pass

      else:

        h = near[ind]

        if h in self.hit:
          continue

        self.fractures[i].append((h, masked_diff[mp,:]))
        print(mp, m, near[ind])
        keep.add(i)
        self.hit.add(h)

    fracs = []
    for i,c in enumerate(self.fractures):
      if i in keep:
        fracs.append(c)
      else:
        self.old_fractures.append(c)

    self.fractures = fracs

    return len(self.fractures)>0

