#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import gtk


BACK = [1,1,1,1]
FRONT = [0,0,0,0.8]
LIGHT = [0,0,0,0.2]
CYAN = [0,0.5,0.5,0.2]
BLUE = [0,0,1,0.3]


NMAX = 10**6
SIZE = 1500
ONE = 1./SIZE
LINEWIDTH = ONE*1.1

INIT_NUM = 20000
INIT_RAD = 0.45

SOURCE_DST = 2.0*ONE

FRAC_DOT = 0.99
FRAC_DST = 200.*ONE
FRAC_STP = ONE*3
FRAC_SPD = 1.0

FRAC_DIMINISH = 1.0
FRAC_SPAWN_DIMINISH = 1.0


SPAWN_ANGLE = 2
SPAWN_FACTOR = 0.8



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

  # for f in alive_fractures:
    # for s in sources[f.inds,:]:
      # render.circle(*s, r=2*ONE, fill=False)

def random_uniform_circle(rad, num):

  from numpy.random import random
  from numpy.linalg import norm
  from numpy import array

  while True:
    xy = 0.5-random(size=2)
    if norm(xy)>1.0:
      continue
    r = array([0.5]*2)+xy*rad
    return r



def main():

  from render.render import Animate
  from modules.fracture import Fractures

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
    domain = 'circ'
  )

  # uniform square distribution
  # from numpy.random import random
  # for _ in xrange(3):
    # F.blow(2, random(size=2))

  # uniform circular distribution
  for _ in xrange(5):
    F.blow(3, random_uniform_circle(INIT_RAD, num=1))

  def wrap(render):

    if F.i % 5 == 0:
      show(render,F)
      # render.write_to_png('{:04d}.png'.format(F.i))

    F.print_stats()
    res = F.step(dbg=True)
    n = F.spawn_front(factor=SPAWN_FACTOR, angle=SPAWN_ANGLE)
    print('spawned: {:d}'.format(n))

    # fn = './asdf_{:04d}.png'.format(F.i)
    # render.write_to_png(fn)

    # from dddUtils.ioOBJ import export_2d as export
    # vertices, paths = F.get_vertices_and_paths()
    # fn = './res/export.2obj'.format(F.i)
    # export('fractures', fn, vertices, lines=paths)

    return res

  render = Animate(SIZE, BACK, FRONT, wrap)

  def __write_svg_and_exit(*args):
    gtk.main_quit(*args)
    show(render,F)
    render.write_to_png('./res/on_exit.png')

    from dddUtils.ioOBJ import export_2d as export
    vertices, paths = F.get_vertices_and_paths()
    fn = './res/on_exit.2obj'.format(F.i)
    export('fractures', fn, vertices, lines=paths)

  render.window.connect("destroy", __write_svg_and_exit)

  gtk.main()


if __name__ == '__main__':

  main()

