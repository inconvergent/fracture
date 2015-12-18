#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import pi

TWOPI = pi*2.

NMAX = 10**6
SIZE = 1500
ONE = 1./SIZE

INIT_NUM = 20000
INIT_RAD = 0.45
INIT_DST = 1*ONE

CRACK_DOT = 0.90
CRACK_DST = INIT_RAD*4

MID = 0.5

LINEWIDTH = ONE*1.1

BACK = [1,1,1,1]
FRONT = [0,0,0,0.9]
RED = [1,0,0,0.3]
BLUE = [0,0,1,0.3]


i = 0


def show(render,f):

  sources = f.sources

  def lines(fractures):

    for frac in fractures:
      s,_ = frac[0]
      render.ctx.move_to(*sources[s,:].flatten())
      for c,_ in frac[1:]:
        render.ctx.line_to(*sources[c,:].flatten())
      render.ctx.stroke()

  render.clear_canvas()
  # render.ctx.set_source_rgba(*RED)
  # for s in sources:
    # render.circle(*s, r=ONE, fill=True)

  render.ctx.set_source_rgba(*FRONT)
  lines(f.fractures + f.old_fractures)


def step(f):

  global i
  i += 1

  res = f.fracture()

  for _ in xrange(40):
    f.make_fracture_from_old()

  return res


def main():

  import gtk
  from render.render import Animate
  from numpy.random import random
  from modules.fracture import Fracture

  F = Fracture(INIT_NUM, INIT_RAD, INIT_DST, CRACK_DOT, CRACK_DST)

  for _ in xrange(1):
    F.blow(1, random(size=2))

  def wrap(render):

    global i

    show(render,F)
    res = step(F)

    i += 1

    return res

  render = Animate(SIZE, BACK, FRONT, wrap)
  render.set_line_width(LINEWIDTH)
  gtk.main()


if __name__ == '__main__':

  main()

