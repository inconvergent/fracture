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
    self.query = fractures.tree.query_ball_point
    self.sources = fractures.sources
    self.start = start
    self.inds = [start]
    self.dxs = [dx]
    self.alive = True

  def __find_near_fractures(self,h):

    x = self.sources[h,:]
    near = self.query(x, self.fractures.frac_dst)
    # print(len(near))

  def step(self):

    fractures = self.fractures
    sources = self.sources
    frac_dst = fractures.frac_dst
    dt = fractures.frac_dot
    hit = fractures.hit

    p = self.inds[-1]
    dx = self.dxs[-1].reshape((1,2))
    px = sources[p,:]

    near = self.query(px, frac_dst)
    diff = sources[near,:] - px
    nrm = norm(diff,axis=1).reshape((-1,1))

    nrm[nrm<=0] = 1e10
    diff /= nrm

    dot = (dx * diff).sum(axis=1)
    mask = dot>dt

    nonz = mask.nonzero()[0]
    masked_sources = sources[mask,:]
    masked_diff = diff[mask]
    masked_nrm = nrm[mask]

    if len(masked_nrm)<1:
      self.alive = False
      return

    mp = masked_nrm.argmin()
    m = masked_nrm[mp]
    ind = nonz[mp]

    if m<=0.0 or m>1.0:
      self.alive = False
      return

    h = near[ind]
    dx = masked_diff[mp,:]

    self.__find_near_fractures(h)

    self.dxs.append(dx)
    self.inds.append(h)

    if not h in hit:
      hit[h] = dx
    else:
      self.alive = False
      return

class Fractures(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      source_dst=0.0, 
      frac_dot=0.95,
      frac_dst=0.05
      ):

    self.init_num = init_num
    self.init_rad = init_rad
    self.source_dst = source_dst 
    self.frac_dot = frac_dot
    self.frac_dst = frac_dst

    self.alive_fractures = []
    self.dead_fractures = []

    self.hit = {}

    self.__make_sources()

  def blow(self,n=5,x=array([0.5,0.5])):

    for a in random(size=n)*TWOPI:

      dx = array([cos(a), sin(a)])
      self.__make_fracture(x=x, dx=dx)

  def __make_sources(self):

    from scipy.spatial import cKDTree as kdt
    from utils import darts

    sources = darts(self.init_num, 0.5, 0.5, self.init_rad, self.source_dst)
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree

  def __make_fracture(self, x, dx):

    _,p = self.tree.query(x,1)
    self.alive_fractures.append(Fracture(self,p,dx))

  def make_random_fracture(self):

    # a = arctan2(dx[1], dx[0]) + (-1)**randint(2) * HPI
    cands = array(self.hit.keys())
    i = cands[randint(len(cands))]

    dx = self.hit[i]
    a = arctan2(dx[1], dx[0])
    x = self.sources[i,:]

    a1 = a - HPI
    dx1 = array([cos(a1), sin(a1)])
    self.__make_fracture(x=x, dx=dx1)

    a2 = a + HPI
    dx2 = array([cos(a2), sin(a2)])
    self.__make_fracture(x=x, dx=dx2)

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

  def print_stats(self):

    alive = len(self.alive_fractures)
    dead = len(self.dead_fractures)
    print('a: {:d} d: {:d}\n'.format(alive, dead))

