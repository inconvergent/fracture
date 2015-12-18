# -*- coding: utf-8 -*-

class Fracture(object):

  def __init__(
      self, 
      init_num, 
      init_rad, 
      init_dst=0.0, 
      crack_angle=0.3,
      crack_dst=0.1
      ):

    self.init_num = init_num
    self.init_rad = init_rad
    self.init_dst = init_dst
    self.crack_angle = crack_angle
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

    _,p = self.tree.query([0.5,0.5], 1)
    self.cracks.append([(p, 0.0)])

  def crack(self):

    from numpy import arctan2
    from numpy import logical_and

    query = self.tree.query_ball_point
    sources = self.sources
    crack_dst = self.crack_dst
    ca = self.crack_angle*0.5

    for i,c in enumerate(self.cracks):

      if c:

        p,a = c.pop()
        px = sources[p,:]

        near = query(px, crack_dst)
        dd = px - sources[near,:]

        aa = arctan2(dd[:,1], dd[:,0])

        mask = logical_and(
          aa>a-ca,
          aa<a+ca
        )
        print(p,a, near, dd, aa, mask)

    return

