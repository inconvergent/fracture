#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import pi
from numpy import cos
from numpy import sin
from numpy import array


BACK = [1,1,1,1]
FRONT = [0,0,0,0.8]
LIGHT = [0,0,0,0.2]
CYAN = [0,0.5,0.5,0.2]
BLUE = [0,0,1,0.3]


NMAX = 10**6
SIZE = 1000
ONE = 1./SIZE
LINEWIDTH = ONE*1.1

INIT_NUM = 20000
INIT_RAD = 0.45
EDGE = 0.5-INIT_RAD

SOURCE_DST = 2.0*ONE

FRAC_DOT = 0.95
FRAC_DST = 100.*ONE
FRAC_STP = ONE
FRAC_SPD = 1.0

FRAC_DIMINISH = 1.0
FRAC_SPAWN_DIMINISH = 0.85


SPAWN_ANGLE = 2.0
SPAWN_FACTOR = 0.2



def show(render,fractures):

  sources = fractures.sources
  alive_fractures = fractures.alive_fractures
  dead_fractures = fractures.dead_fractures

  def draw_sources():
    for i,s in enumerate(sources):
      if i not in fractures.visited:
        render.circle(*s, r=4*ONE, fill=True)

  def draw_lines(fracs):
    for frac in fracs:
      start = frac.inds[0]
      render.ctx.move_to(*sources[start,:])
      for c in frac.inds[1:]:
        render.ctx.line_to(*sources[c,:])
      render.ctx.stroke()

  render.clear_canvas()

  # render.ctx.set_source_rgba(*LIGHT)
  # draw_sources()

  render.ctx.set_source_rgba(*LIGHT)
  render.set_line_width(3*LINEWIDTH)
  draw_lines(alive_fractures+dead_fractures)

  render.ctx.set_source_rgba(*FRONT)
  render.set_line_width(LINEWIDTH)
  draw_lines(alive_fractures+dead_fractures)

  render.ctx.set_source_rgba(*FRONT)
  render.set_line_width(4*LINEWIDTH)
  render.ctx.move_to(1.0-EDGE, 1.0-EDGE)
  render.ctx.line_to(1.0-EDGE, EDGE)
  render.ctx.line_to(EDGE, EDGE)
  render.ctx.line_to(EDGE, 1.0-EDGE)
  render.ctx.close_path()
  render.ctx.stroke()



  # for f in alive_fractures:
    # for s in sources[f.inds,:]:
      # render.circle(*s, r=2*ONE, fill=False)


def main():

  from render.render import Animate
  from modules.fracture import Fractures

  from dddUtils.ioOBJ import export_2d as export
  from fn import Fn
  fn = Fn(prefix='./res/',postfix='.2obj')

  F = Fractures(
      INIT_NUM,
      INIT_RAD,
      SOURCE_DST,
      FRAC_DOT,
      FRAC_DST,
      FRAC_STP,
      FRAC_SPD,
      FRAC_DIMINISH,
      FRAC_SPAWN_DIMINISH,
      domain = 'rect'
      )

  print(F.sources.shape)

  from numpy.random import random
  for _ in xrange(60):
    a = -pi*0.5
    dx = array([cos(a), sin(a)])
    x = [EDGE + random()*2*INIT_RAD, 1.0-EDGE]
    F.crack(x, dx)

  for _ in xrange(60):
    a = pi*0.5
    dx = array([cos(a), sin(a)])
    x = [EDGE + random()*2*INIT_RAD, EDGE]
    F.crack(x, dx)

  def wrap(render):

    if F.i % 5 == 0:
      show(render,F)
      vertices, paths = F.get_vertices_and_paths()
      # export('fractures', fn.name(), vertices, lines=paths)
      # render.write_to_png('{:04d}.png'.format(F.i))

    F.print_stats()
    res = F.step(dbg=False)

    if not res:
      # vertices, paths = F.get_vertices_and_paths()
      # export('fractures', fn.name(), vertices, lines=paths)
      show(render,F)
      render.write_to_png(fn.name()+'.png')

    return res

  render = Animate(SIZE, BACK, FRONT, wrap)
  render.start()


if __name__ == '__main__':

  main()

