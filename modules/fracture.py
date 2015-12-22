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
from collections import defaultdict

TWOPI = pi*2
HPI = pi*0.5

class Fracture(object):

  def __init__(self, fractures, start, dx):

    self.i = 0
    self.fractures = fractures
    self.tree = fractures.tree

    self.hit = set([start])
    self.start = start
    self.inds = [start]
    self.dxs = [dx]
    self.alive = True

  def __find_near_fractures(self, c, h):

    from numpy import column_stack
    from numpy import logical_not
    from operator import itemgetter

    sources = self.fractures.sources

    cx = sources[c,:]
    hx = sources[h,:]
    mid = 0.5*(hx+cx)
    u = array(
      self.tree.query_ball_point(mid,self.fractures.frac_dst),
      'int'
    )

    uc = norm(cx-sources[u,:], axis=1)
    uh = norm(hx-sources[u,:], axis=1)

    ch = norm(cx-hx)
    mm = column_stack([uc,uh]).max(axis=1)
    mask = ch<mm

    a = set(u[logical_not(mask)])
    relative_neigh_sources = a.difference([c,h])
    relative_neigh_sources_hit = relative_neigh_sources.intersection(
      self.fractures.hit
    )

    if relative_neigh_sources_hit:
      rnh = [ (r, norm(sources[r,:]-mid)) for r in relative_neigh_sources_hit]
      rnh.sort(key=itemgetter(1))
      return rnh[0][0]

    else:
      return -1

  def step(self, add_sources=False):

    self.i += 1

    fractures = self.fractures
    sources = fractures.sources
    frac_dst = fractures.frac_dst
    dt = fractures.frac_dot
    hit = fractures.hit

    c = self.inds[-1]
    dx = self.dxs[-1].reshape((1,2))
    cx = sources[c,:]

    near = self.tree.query_ball_point(cx, frac_dst)
    diff = sources[near,:] - cx
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
      return False

    mp = masked_nrm.argmin()
    m = masked_nrm[mp]
    ind = nonz[mp]

    if m<=0.0 or m>1.0:
      self.alive = False
      return False

    h = near[ind]
    dx = masked_diff[mp,:]

    collide = self.__find_near_fractures(c, h)
    if collide>-1:
      self.alive = False
      h = collide

    if not h in hit:
      hit[h] = dx
    else:
      self.alive = False

    self.dxs.append(dx)
    self.inds.append(h)

    if add_sources:
      diff = sources[h,:] - sources[c,:]
      nrm = norm(diff)
      source_dst = self.fractures.source_dst*5
      d = arange(source_dst,nrm,source_dst)
      new_sources = sources[c,:]+diff/nrm*reshape(d, (-1,1))
      if len(new_sources)>0:
        self.fractures.add_hit_sources(new_sources, dx)
        # print(new_sources)

    return self.alive

class Fractures(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      source_dst=0.0, 
      frac_dot=0.95,
      frac_dst=0.05
    ):

    self.i = 0
    self.init_num = init_num
    self.init_rad = init_rad
    self.source_dst = source_dst 
    self.frac_dot = frac_dot
    self.frac_dst = frac_dst

    self.alive_fractures = []
    self.dead_fractures = []

    self.hit = {}

    self.__make_sources()

  def blow(self,n=5,x=array([0.5,0.5]), add_sources=False):

    for a in random(size=n)*TWOPI:
      dx = array([cos(a), sin(a)])
      self.__make_fracture(x, dx, add_sources=add_sources)

  def __make_sources(self):

    from scipy.spatial import cKDTree as kdt
    from utils import darts

    sources = darts(self.init_num, 0.5, 0.5, self.init_rad, self.source_dst)
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree

  def more_sources(self,n):

    from scipy.spatial import cKDTree as kdt
    from utils import darts
    
    sources = darts(
      n, 
      0.5, 
      0.5, 
      self.init_rad, 
      self.source_dst, 
      self.sources
    )
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree

  def add_hit_sources(self, src, dx):

    from numpy import row_stack
    from scipy.spatial import cKDTree as kdt

    sources = self.sources

    n1 = len(sources)
    sources = row_stack([sources, src])
    n2 = len(sources)
    tree = kdt(sources)

    for n in range(n1,n2):
      self.hit[n] = dx

    # self.hit.update(range(n1,n2))
    self.sources = sources
    self.tree = tree

  def __make_fracture(self, x, dx, add_sources=False):

    _,p = self.tree.query(x,1)
    f = Fracture(self,p,dx)
    res = f.step(add_sources)
    if res:
      self.alive_fractures.append(f)
    return res

  def make_random_alive_fracture(self, i, angle=0.3, add_sources=False):

    if not self.alive_fractures:
      return False

    dx = self.alive_fractures[i].dxs[-1]
    a = arctan2(dx[1], dx[0])
    x = self.sources[self.alive_fractures[i].inds[-1],:]

    a1 = a + (-1)**randint(2)*HPI + (0.5-random()) * angle
    return self.__make_fracture(x, array([cos(a1), sin(a1)]), add_sources)

  def make_random_fracture(self, i, angle=0, add_sources=False):

    ## i ∈ self.sources
    dx = self.hit[i]
    a = arctan2(dx[1], dx[0])
    x = self.sources[i,:]

    # a1 = a + (-1)**randint(2)*HPI
    a1 = a + (-1)**randint(2)*HPI + (0.5-random()) * angle
    dx1 = array([cos(a1), sin(a1)])

    return self.__make_fracture(x, dx1, add_sources)


  def step(self, add_sources=False):

    self.i += 1

    fracs = []
    for f in self.alive_fractures:
      f.step(add_sources)
      if f.alive:
        fracs.append(f)
      else:
        if len(f.inds)>1:
          self.dead_fractures.append(f)
        else:
          print('discarding path')

    self.alive_fractures = fracs
    return len(fracs)>0

  def get_fracture_paths(self):

    from numpy import row_stack

    paths = []

    for f in self.alive_fractures + self.dead_fractures:
      if len(f.inds)<2:
        continue
      path = row_stack([self.sources[p,:] for p in f.inds])
      paths.append(path)

    return paths

  def print_stats(self):

    alive = len(self.alive_fractures)
    dead = len(self.dead_fractures)
    print('a: {:d} d: {:d} s: {:d}\n'
      .format(alive, dead, len(self.sources))
    )

