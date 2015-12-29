# -*- coding: utf-8 -*-

from numpy import pi
from numpy import array
from numpy import arange
from numpy import reshape
from numpy import row_stack
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

  def __init__(self, fractures, fid, start, dx):

    self.i = 0
    self.fractures = fractures
    self.tree = fractures.tree

    self.start = start
    self.inds = [start]
    self.dxs = [dx]
    self.alive = True

    self.fid = fid

  def __find_near_fractures(self, c, h):

    from numpy import column_stack
    from numpy import logical_not
    from operator import itemgetter

    sources = self.fractures.sources

    cx = sources[c,:]
    hx = sources[h,:]
    u = array(
      self.tree.query_ball_point(0.5*(hx+cx),self.fractures.frac_dst),
      'int'
    )

    uc = norm(cx-sources[u,:], axis=1)
    uh = norm(hx-sources[u,:], axis=1)

    ch = norm(cx-hx)
    mm = column_stack([uc,uh]).max(axis=1)
    mask = ch<mm

    a = set(u[logical_not(mask)])
    b = set([c,h])
    relative_neigh_sources = a.difference(b)
    relative_neigh_sources_visited = relative_neigh_sources.intersection(
      self.fractures.visited
    )

    if relative_neigh_sources_visited:
      rnh = [ (r, norm(sources[r,:]-cx)) for r in relative_neigh_sources_visited]
      rnh.sort(key=itemgetter(1))
      return rnh[0][0]

    else:
      return -1

  def step(self, dbg=False):

    self.i += 1
    dbgs = ''

    fractures = self.fractures
    sources = fractures.sources
    frac_dst = fractures.frac_dst
    dt = fractures.frac_dot
    visited = fractures.visited
    stp = fractures.frac_stp

    c = self.inds[-1]
    cx = sources[c,:]
    cdx = self.dxs[-1].reshape((1,2))

    near = self.tree.query_ball_point(cx, frac_dst)

    neardiff = sources[near,:] - cx
    nearnrm = norm(neardiff,axis=1).reshape((-1,1))

    nearnrm[nearnrm<=1e-9] = 1e10
    neardiff /= nearnrm

    mask = (cdx*neardiff).sum(axis=1)>dt

    if mask.sum()<1:
      self.alive = False
      if dbg:
        print(self.fid, 'no nearby sources')
      return False

    masked_diff = neardiff[mask]
    masked_nrm = nearnrm[mask]

    new_dx = (masked_diff/masked_nrm).sum(axis=0).flatten()
    new_dx /= norm(new_dx)

    ind = masked_nrm.argmin()

    if masked_nrm[ind]<stp:

      nonz = mask.nonzero()[0]
      h = near[nonz[ind]]
      
      if h in visited:
        # collision
        dbgs += '{:d}: {:s}'.format(self.fid, 'collision')
        self.alive = False

      else: 
        # no collision

        # this is pretty bad
        collision = self.__find_near_fractures(c,h)
        if collision>-1:
          dbgs += '{:d}: {:s}'.format(self.fid, 'collision (rn)')
          h = collision
          self.alive = False
        else:
          dbgs += '{:d}: {:s}'.format(self.fid, 'no collision')
          self.alive = True
          visited[h] = new_dx

    else:
      # new source
      dbgs += '{:d}: {:s}'.format(self.fid, 'new source')
      new_pos = cx + new_dx*fractures.frac_stp
      h = self.fractures._add_tmp_source(new_pos)
      self.alive = True
      visited[h] = new_dx

    if dbg:
      print(dbgs)

    self.dxs.append(new_dx)
    self.inds.append(h)

    return self.alive

class Fractures(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      source_dst,
      frac_dot,
      frac_dst,
      frac_stp
    ):

    self.i = 0
    self.init_num = init_num
    self.init_rad = init_rad
    self.source_dst = source_dst 
    self.frac_dot = frac_dot
    self.frac_dst = frac_dst
    self.frac_stp = frac_stp

    self.alive_fractures = []
    self.dead_fractures = []

    self.visited = {}

    self.count = 0

    self.tmp_sources = []
    self.__make_sources()

  def blow(self,n, x=array([0.5,0.5])):

    self.tmp_sources = []

    for a in random(size=n)*TWOPI:
      dx = array([cos(a), sin(a)])
      self.__make_fracture(x, dx)

    self._append_tmp_sources()

  def __make_sources(self):

    from scipy.spatial import cKDTree as kdt
    from utils import darts

    sources = darts(self.init_num, 0.5, 0.5, self.init_rad, self.source_dst)
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree

    return len(sources)

  def _add_tmp_source(self, x):

    self.tmp_sources.append(x)
    return len(self.sources)+len(self.tmp_sources)-1

  def _append_tmp_sources(self):

    from scipy.spatial import cKDTree as kdt

    sources = row_stack([self.sources]+self.tmp_sources)
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree
    self.tmp_sources = []

    return len(sources)

  def __make_fracture(self, x, dx):

    _,p = self.tree.query(x,1)
    f = Fracture(self,self.count,p,dx)
    self.count += 1
    res = f.step()
    if res:
      self.alive_fractures.append(f)
    return res

  def make_random_alive_fracture(self, i, angle=0):

    if not self.alive_fractures:
      return False

    dx = self.alive_fractures[i].dxs[-1]
    a = arctan2(dx[1], dx[0])
    x = self.sources[self.alive_fractures[i].inds[-1],:]

    a1 = a + (-1)**randint(2)*HPI + (0.5-random()) * angle
    return self.__make_fracture(x=x, dx=array([cos(a1), sin(a1)]))

  def spawn(self, factor, angle=0.7):

    self.tmp_sources = []

    count = 0
    for i in (random(size=len(self.alive_fractures))<factor).nonzero()[0]:
      count += int(self.make_random_alive_fracture(i, angle))

    self._append_tmp_sources()

    return count

  def step(self):

    self.i += 1

    self.tmp_sources = []

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

    self._append_tmp_sources()

    return len(fracs)>0

  def get_fracture_paths(self):

    paths = []

    for f in self.alive_fractures + self.dead_fractures:
      if len(f.inds)<2:
        continue
      path = row_stack([self.sources[p,:] for p in f.inds])
      paths.append(path)

    return paths

  def get_vertices_and_edges(self):

    vertices = self.sources
    edges = []
    for f in self.alive_fractures + self.dead_fractures:
      if len(f.inds)<2:
        continue
      for i in xrange(len(f.inds)-1):
        edges.append([f.inds[i], f.inds[i+1]])

    return vertices, row_stack(edges)

  def print_stats(self):

    alive = len(self.alive_fractures)
    dead = len(self.dead_fractures)
    print('# {:d} a: {:d} d: {:d} s: {:d}\n'
      .format(self.i, alive, dead, len(self.sources))
    )

