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
    u = array(
      self.tree.query_ball_point(0.5*(hx+cx), 
      self.fractures.frac_dst),'int'
    )

    uc = norm(cx-sources[u,:], axis=1)
    uh = norm(hx-sources[u,:], axis=1)

    ch = norm(cx-hx)
    mm = column_stack([uc,uh]).max(axis=1)
    mask = ch<mm

    a = set(u[logical_not(mask)])
    b = set([c,h])
    relative_neigh_sources = a.difference(b)
    relative_neigh_sources_hit = relative_neigh_sources.intersection(
      self.fractures.hit
    )

    if relative_neigh_sources_hit:
      rnh = [ (r, norm(sources[r,:]-cx)) for r in relative_neigh_sources_hit]
      rnh.sort(key=itemgetter(1))
      return rnh[0][0]

    else:
      return -1

  def step(self):

    self.i += 1

    fractures = self.fractures
    sources = fractures.sources
    frac_dst = fractures.frac_dst
    dt = fractures.frac_dot
    hit = fractures.hit
    count = fractures.count

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
      self.dxs.append(dx)
      self.inds.append(collide)
      return False

    self.dxs.append(dx)
    self.inds.append(h)

    if h in self.hit:
      print('warning')
    self.hit.add(h)

    if not h in hit:
      hit[h] = dx
      count[h] += 1
    else:
      self.alive = False
      return False

    return True 

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
    self.count = defaultdict(int)

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
    self.count[p] += 1
    f = Fracture(self,p,dx)
    res = f.step()
    if res:
      self.alive_fractures.append(f)
    return res

  def make_random_alive_fracture(self):

    if not self.alive_fractures:
      return False

    cands = arange(len(self.alive_fractures))
    i = cands[randint(len(cands))]

    dx = self.alive_fractures[i].dxs[-1]
    a = arctan2(dx[1], dx[0])
    x = self.sources[self.alive_fractures[i].inds[-1],:]

    # a1 = a + (0.5-random()) * pi
    if random()<0.5:
      a1 = a - HPI + (0.5-random()) * 0.3
    else:
      a1 = a + HPI + (0.5-random()) * 0.3
    dx1 = array([cos(a1), sin(a1)])

    return self.__make_fracture(x=x, dx=dx1)

  def make_random_fracture(self):

    cands = array(self.hit.keys())
    i = cands[randint(len(cands))]

    dx = self.hit[i]
    a = arctan2(dx[1], dx[0])
    x = self.sources[i,:]

    if random()<0.5:
      a1 = a - HPI
    else:
      a1 = a + HPI
    dx1 = array([cos(a1), sin(a1)])

    return self.__make_fracture(x=x, dx=dx1)


  def step(self):

    self.i += 1

    fracs = []
    for f in self.alive_fractures:
      f.step()
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

