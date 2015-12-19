#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import pi

TWOPI = pi*2.

NMAX = 10**6
SIZE = 1500
ONE = 1./SIZE

INIT_NUM = 10000
INIT_RAD = 0.45

SOURCE_DST = 2*ONE

FRAC_DOT = 0.92
FRAC_DST = SOURCE_DST*100

MID = 0.5

LINEWIDTH = ONE*1.1

BACK = [1,1,1,1]
FRONT = [0,0,0,0.9]
RED = [1,0,0,0.3]
CYAN = [0,0.5,0.5,0.05]
BLUE = [0,0,1,0.3]



def show(render,fractures):

  sources = fractures.sources

  def draw_sources():
    render.ctx.set_source_rgba(*RED)
    for s in sources:
      render.circle(*s, r=ONE, fill=True)

  def draw_lines(fracs):

    render.ctx.set_source_rgba(*FRONT)

    for frac in fracs:
      start = frac.inds[0]
      render.ctx.move_to(*sources[start,:])
      for c in frac.inds[1:]:
        render.ctx.line_to(*sources[c,:])
      render.ctx.stroke()

  render.clear_canvas()
  render.ctx.set_source_rgba(*FRONT)
  draw_lines(fractures.alive_fractures + fractures.dead_fractures)


def step(fractures):

  from modules.utils import export_svg

  fractures.print_stats()

  res = fractures.step()

  for _ in xrange(5):
    fractures.make_random_fracture()

  paths = fractures.get_fracture_paths()

  # fn = './res/asdf_{:05d}.svg'.format(i)
  # export_svg(fn, paths, SIZE)

  return res


def main():

  import gtk
  from render.render import Animate
  from numpy.random import random
  from modules.fracture import Fractures

  F = Fractures(INIT_NUM, INIT_RAD, SOURCE_DST, FRAC_DOT, FRAC_DST)

  ## init
  for _ in xrange(1):
    F.blow(1, random(size=2))

  def wrap(render):

    show(render,F)
    return step(F)

  render = Animate(SIZE, BACK, FRONT, wrap)
  render.set_line_width(LINEWIDTH)
  gtk.main()


if __name__ == '__main__':

  main()

