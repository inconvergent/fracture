#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import gtk

NMAX = 10**6
SIZE = 1500
ONE = 1./SIZE

INIT_NUM = 100000
INIT_RAD = 0.45

SOURCE_DST = 1.5*ONE

FRAC_DOT = 0.90
FRAC_DST = 200.*ONE
FRAC_STP = ONE*10

LINEWIDTH = ONE*1.1

BACK = [1,1,1,1]
FRONT = [0,0,0,0.8]
LIGHT = [0,0,0,0.2]
RED = [1,0,0,0.1]
CYAN = [0,0.5,0.5,0.05]
BLUE = [0,0,1,0.3]



def show(render,fractures):

  sources = fractures.sources
  alive_fractures = fractures.alive_fractures
  dead_fractures = fractures.dead_fractures

  def draw_sources():
    for s in sources:
      render.circle(*s, r=ONE, fill=True)

  def draw_lines(fracs):
    for frac in fracs:
      start = frac.inds[0]
      render.ctx.move_to(*sources[start,:])
      for c in frac.inds[1:]:
        render.ctx.line_to(*sources[c,:])
      render.ctx.stroke()

  render.clear_canvas()

  # render.ctx.set_source_rgba(*RED)
  # draw_sources()

  render.ctx.set_source_rgba(*FRONT)
  render.set_line_width(LINEWIDTH)
  draw_lines(alive_fractures+dead_fractures)


  for f in alive_fractures:
    for s in sources[f.inds,:]:
      render.circle(*s, r=2*ONE, fill=False)



def main():

  from render.render import Animate
  from numpy.random import random
  from modules.fracture import Fractures

  F = Fractures(
    INIT_NUM, 
    INIT_RAD, 
    SOURCE_DST, 
    FRAC_DOT, 
    FRAC_DST, 
    FRAC_STP
  )

  ## init
  for _ in xrange(1):
    F.blow(10, random(size=2))

  def wrap(render):

    show(render,F)
    F.print_stats()
    res = F.step()
    n = F.spawn(factor=0.1, angle=0.7)
    print('spawned: {:d}'.format(n))

    # from modules.utils import export_svg
    # paths = F.get_fracture_paths()
    # fn = './res/asdf_{:05d}.svg'.format(F.i)
    # export_svg(fn, paths, SIZE)

    # fn = './asdf_{:05d}.png'.format(F.i)
    # render.write_to_png(fn)

    return res

  render = Animate(SIZE, BACK, FRONT, wrap)

  def __write_svg_and_exit(*args):                                                                                                                                                                                                                                            
    from modules.utils import export_svg
    export_svg('./res/on_exit.svg', F.get_fracture_paths(), SIZE)
    gtk.main_quit(*args)
  render.window.connect("destroy", __write_svg_and_exit)

  gtk.main()


if __name__ == '__main__':

  main()

