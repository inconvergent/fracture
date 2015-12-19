# -*- coding: utf-8 -*-

from numpy import pi
from numpy import array
from numpy import arange
from numpy import reshape
from numpy.random import random
from numpy.random import randint
from numpy.linalg import norm
from numpy import cos
from numpy import sin
from numpy import arctan2

TWOPI = pi*2
HPI = pi*0.5

class Fracture(object):

  def __init__(self, fractures, start, dx):
    
    self.fractures = fractures
    self.start = start
    self.inds = [start]
    self.dxs = [dx]
    self.alive = True

  def step(self):

    fractures = self.fractures
    sources= fractures.sources
    query = fractures.tree.query_ball_point
    frac_dst = fractures.frac_dst
    dt = fractures.frac_dot
    hit = fractures.hit

    try:

      p = self.inds[-1]
      dx = self.dxs[-1].reshape((1,2))
      px = sources[p,:]

      near = query(px, frac_dst)
      diff = sources[near,:] - px
      nrm = norm(diff,axis=1).reshape((-1,1))

      zeromask = nrm<=0
      nrm[zeromask] = 1000.0

      diff /= nrm

      dot = (dx * diff).sum(axis=1)
      mask = dot>dt

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

      self.dxs.append(masked_diff[mp,:])
      self.inds.append(h)

      if not h in hit:
        hit[h] = masked_diff[mp,:]
      else:
        self.alive = False

class Fractures(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      init_dst=0.0, 
      frac_dot=0.95,
      frac_dst=0.05
      ):

    self.init_num = init_num
    self.init_rad = init_rad
    self.init_dst = init_dst
    self.frac_dot = frac_dot
    self.frac_dst = frac_dst

    self.alive_fractures = []
    self.dead_fractures = []

    self.hit = {}

    self.__make_sources()

  def blow(self,n=5,x=array([0.5,0.5])):

    for a in random(size=n)*TWOPI:

      dx = array([cos(a), sin(a)])
      self.make_fracture(x=x, dx=dx)

  def __make_sources(self):

    from scipy.spatial import cKDTree as kdt
    from utils import darts

    sources = darts(self.init_num, 0.5, 0.5, self.init_rad, self.init_dst)
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree

  def make_fracture(self, x, dx):

    _,p = self.tree.query(x,1)
    self.alive_fractures.append(Fracture(self,p,dx))

  def make_fracture_from_old(self):

    cands = array(self.hit.keys())
    i = cands[randint(len(cands))]

    dx = self.hit[i]
    a = arctan2(dx[1], dx[0]) + (-1)**randint(2) * HPI
    dx = array([cos(a), sin(a)])

    x = self.sources[i,:]
    self.make_fracture(x=x, dx=dx)

  def step(self):

    fracs = []
    for f in self.alive_fractures:
      f.step()
      if f.alive:
        fracs.append(f)
      else:
        self.dead_fractures.append(f)

    self.alive_fractures = fracs
    return len(fracs)>0

  def get_fracture_paths(self):

    from numpy import row_stack

    paths = []

    for f in self.alive_fractures + self.dead_fractures:
      path = row_stack([self.sources[p,:] for p in f.inds])
      paths.append(path)

    return paths

