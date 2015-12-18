# -*- coding: utf-8 -*-

from numpy import pi

TWOPI = pi*2
HPI = pi*0.5

class Fracture(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      init_dst=0.0, 
      crack_dot=0.9,
      crack_dst=0.05
      ):

    self.init_num = init_num
    self.init_rad = init_rad
    self.init_dst = init_dst
    self.crack_dot = crack_dot
    self.crack_dst = crack_dst

    self.cracks = []


    self.__make_sources()
    self.__make_cracks()

  def __make_sources(self):

    from scipy.spatial import cKDTree as kdt
    from utils import darts

    sources = darts(self.init_num, 0.5, 0.5, self.init_rad, self.init_dst)
    tree = kdt(sources)

    self.sources = sources
    self.tree = tree

  def __make_cracks(self):

    from numpy import array

    _,p = self.tree.query([0.5,0.5], 1)
    self.cracks.append([(p, array([0.0,1.0]))])

  def fracture(self):

    from numpy import arctan2
    from numpy import logical_and
    from numpy import arange
    from numpy import array
    from numpy import reshape
    from numpy import abs
    from numpy.linalg import norm

    query = self.tree.query_ball_point
    sources = self.sources
    crack_dst = self.crack_dst
    dl = self.crack_dot

    keep = set()

    for i in range(len(self.cracks)):

      try:

        p,dx = self.cracks[i][-1]
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

        self.cracks[i].append((near[ind], masked_diff[mp,:]))
        print(mp, m, near[ind])
        keep.add(i)

    self.cracks = [c for i,c in enumerate(self.cracks) if i in keep]

