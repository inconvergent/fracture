#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import gtk

NMAX = 10**6
SIZE = 1500
ONE = 1./SIZE

INIT_NUM = 10000
INIT_RAD = 0.45

SOURCE_DST = 1.*ONE

FRAC_DOT = 0.97
FRAC_DST = 200*ONE

LINEWIDTH = ONE*1.1

BACK = [1,1,1,1]
FRONT = [0,0,0,0.8]
LIGHT = [0,0,0,0.2]
RED = [1,0,0,0.1]
CYAN = [0,0.5,0.5,0.2]
BLUE = [0,0,1,0.3]



def show(render, fractures):

  sources = fractures.sources

  def draw_sources():
    for s in sources:
      render.circle(*s, r=3*ONE, fill=True)

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

  # render.ctx.set_source_rgba(*LIGHT)
  # render.set_line_width(LINEWIDTH*5)
  # draw_lines(fractures.alive_fractures + fractures.dead_fractures)

  render.ctx.set_source_rgba(*FRONT)
  render.set_line_width(LINEWIDTH)
  draw_lines(fractures.alive_fractures + fractures.dead_fractures)

  for f in fractures.alive_fractures:
    for s in sources[f.inds,:]:
      render.circle(*s, r=2.*ONE, fill=False)

ADD_SOURCES = True


def step(fractures):

  from modules.utils import export_svg
  from numpy import array
  from numpy.random import randint

  fractures.print_stats()

  res = fractures.step(add_sources=ADD_SOURCES)

  fractures.more_sources(500)

  # hit = array(list(fractures.hit), 'int')
  # num = len(hit)
  # for i in hit[randint(0,num,size=100)]:
    # fractures.make_random_fracture(i)

  if not fractures.alive_fractures:

    return False
  count = 0 
  for i in xrange(len(fractures.alive_fractures)):
    spawned = fractures.make_random_alive_fracture(i, add_sources=ADD_SOURCES)
    if spawned:
      count += 1
  print('spawned: {:d}'.format(count))

  return res


def main():

  from render.render import Animate
  from numpy.random import random
  from modules.fracture import Fractures

  fractures = Fractures(INIT_NUM, INIT_RAD, SOURCE_DST, FRAC_DOT, FRAC_DST)

  ## init
  for _ in xrange(3):
    fractures.blow(5, random(size=2), add_sources = ADD_SOURCES)

  def wrap(render):

    show(render, fractures)
    # paths = fractures.get_fracture_paths()
    # fn = './res/asdf_{:05d}.svg'.format(fractures.i)
    # export_svg(fn, paths, SIZE)
    # fn = './res/asdf_{:05d}.png'.format(fractures.i)
    # render.write_to_png(fn)
    return step(fractures)

  render = Animate(SIZE, BACK, FRONT, wrap)

  def __write_svg_and_exit(*args):                                                                                                                                                                                                                                            
    from modules.utils import export_svg
    export_svg('./res/on_exit.svg', fractures.get_fracture_paths(), SIZE)
    gtk.main_quit(*args)
  render.window.connect("destroy", __write_svg_and_exit)

  gtk.main()


if __name__ == '__main__':

  main()

