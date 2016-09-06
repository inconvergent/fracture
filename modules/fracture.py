# -*- coding: utf-8 -*-



from numpy import pi
from numpy import array
from numpy import row_stack
from numpy.random import random
from numpy.random import randint
from numpy.linalg import norm
from numpy import cos
from numpy import sin
from numpy import arctan2


TWOPI = pi*2
HPI = pi*0.5


class Fracture(object):

  def __init__(
      self,
      fractures,
      fid,
      start,
      dx,
      frac_spd,
      frac_diminish
    ):

    self.i = 0
    self.fractures = fractures
    self.tree = fractures.tree
    self.frac_spd = frac_spd
    self.frac_diminish = frac_diminish

    self.start = start
    self.inds = [start]
    self.dxs = [dx]
    self.alive = True

    self.fid = fid

  def __relative_neigh_test(self, curr, new):

    from numpy import concatenate
    from numpy import unique
    from scipy.spatial.distance import cdist

    sources = self.fractures.sources
    visited = self.fractures.visited
    tri = self.fractures.tri
    simplices = tri.simplices
    simp = tri.find_simplex(new,bruteforce=True,tol=1e-10)
    neigh = concatenate((tri.neighbors[simp],[simp]))
    vv = set(list(unique(simplices[neigh,:])))

    if curr in vv:
      vv.remove(curr)
    vv = array(list(vv))

    dist = cdist(sources[vv, :], row_stack([new,sources[curr,:]]))
    mas = dist.max(axis=1)

    # curr_new = norm(new-sources[curr,:])
    curr_new = self.fractures.frac_stp

    free = mas<curr_new

    if sum(free)==0:
      return -1
    else:
      col = [k for k in vv[free] if k in visited]
      if col:
        return col[0]
      else:
        return -1

  def step(self, dbg=False):

    self.i += 1
    self.frac_spd *= self.frac_diminish

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
    new_pos = cx + new_dx*stp

    rel = self.__relative_neigh_test(c, new_pos)

    if rel>-1:
      dbgs += '{:d}: {:s}'.format(self.fid, 'collision (rn)')
      h = rel
      self.alive = False
    else:
      # new source
      dbgs += '{:d}: {:s}'.format(self.fid, 'new source')
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
      frac_stp,
      frac_spd=1.0,
      frac_diminish=1.0,
      frac_spawn_diminish=1.0,
      domain='rect'
    ):

    self.i = 0
    self.init_num = init_num
    self.init_rad = init_rad
    self.source_dst = source_dst
    self.frac_dot = frac_dot
    self.frac_dst = frac_dst
    self.frac_stp = frac_stp
    self.frac_spd = frac_spd
    self.frac_diminish = frac_diminish
    self.spawn_diminish = frac_spawn_diminish

    self.alive_fractures = []
    self.dead_fractures = []

    self.visited = {}

    self.count = 0

    self.tmp_sources = []
    self.__make_sources(domain=domain)

  def blow(self,n, x=array([0.5,0.5])):

    self.tmp_sources = []

    for a in random(size=n)*TWOPI:
      dx = array([cos(a), sin(a)])
      self.__make_fracture(x=x, dx=dx)

    self._append_tmp_sources()

  def __make_sources(self, xx=0.5, yy=0.5, rad=None, domain='rect'):

    from scipy.spatial import cKDTree as kdt
    from scipy.spatial import Delaunay as triag
    from iutils.random import darts
    from iutils.random import darts_rect

    if rad is None:
      rad = self.init_rad

    if domain=='circ':
      sources = darts(
        self.init_num,
        xx,
        yy,
        self.init_rad,
        self.source_dst
      )
    elif domain=='rect':
      sources = darts_rect(
        self.init_num,
        xx,
        yy,
        2*rad,
        2*rad,
        self.source_dst
      )
    else:
      raise ValueError('domain must be "rect" or "circ".')
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree
    self.tri = triag(
      self.sources,
      incremental=False,
      qhull_options='QJ Qc'
    )
    self.num_sources = len(self.sources)

    return len(sources)

  def _add_tmp_source(self, x):

    self.tmp_sources.append(x)
    return len(self.sources)+len(self.tmp_sources)-1

  def _append_tmp_sources(self):

    from scipy.spatial import cKDTree as kdt
    from scipy.spatial import Delaunay as triag

    sources = row_stack([self.sources]+self.tmp_sources)
    tree = kdt(sources)
    self.sources = sources
    self.tree = tree
    self.tmp_sources = []
    self.tri = triag(
      self.sources,
      incremental=False,
      qhull_options='QJ Qc'
    )
    self.num_sources = len(self.sources)

    return len(sources)

  def __make_fracture(self, x=None, p=None, dx=None, spd=None):

    if p is None:
      _,p = self.tree.query(x,1)

    if spd is None:
      spd = self.frac_spd

    f = Fracture(
      self,
      self.count,
      p,
      dx,
      spd,
      self.frac_diminish
    )
    self.count += 1
    res = f.step()
    if res:
      self.alive_fractures.append(f)
    return res

  # def spawn_front(self, factor=1.0, angle=0.7):

    # if not self.alive_fractures:
      # return 0

    # self.tmp_sources = []
    # count = 0

    # for i in (random(size=len(self.alive_fractures))<factor).nonzero()[0]:
      # f = self.alive_fractures[i]
      # dx = f.dxs[-1]
      # a = arctan2(dx[1], dx[0]) + (-1)**randint(2)*HPI + (0.5-random()) * angle
      # count += int(self.__make_fracture(p=f.inds[-1], dx=array([cos(a), sin(a)])))

    # self._append_tmp_sources()

    # return count

  def spawn_front(self, factor=1.0, angle=0.7):

    if not self.alive_fractures:
      return 0

    self.tmp_sources = []
    count = 0

    for i,rnd in enumerate(random(size=len(self.alive_fractures))):
      f = self.alive_fractures[i]

      if rnd>f.frac_spd*factor:
        continue

      dx = f.dxs[-1]
      a = arctan2(dx[1], dx[0]) + (-1)**randint(2)*HPI + (0.5-random()) * angle
      count += int(
        self.__make_fracture(
          p=f.inds[-1],
          dx=array([cos(a), sin(a)]),
          spd=f.frac_spd*self.spawn_diminish
        )
      )

    self._append_tmp_sources()

    return count

  def step(self, dbg=False):

    self.i += 1

    self.tmp_sources = []

    fracs = []
    for f in self.alive_fractures:
      f.step(dbg)
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

  def get_vertices_and_paths(self):

    vertices = self.sources
    paths = []
    for f in self.alive_fractures + self.dead_fractures:
      if len(f.inds)<2:
        continue

      paths.append(array(f.inds, 'int'))

    return vertices, paths

  def print_stats(self):

    alive = len(self.alive_fractures)
    dead = len(self.dead_fractures)
    print('# {:d} a: {:d} d: {:d} s: {:d}\n'
      .format(self.i, alive, dead, len(self.sources))
    )

