#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import pi

TWOPI = pi*2.

NMAX = 10**6
SIZE = 1000
ONE = 1./SIZE

INIT_NUM = 4000
INIT_RAD = 0.4
INIT_DST = 5*ONE

MID = 0.5

LINEWIDTH = 5.*ONE

BACK = [1,1,1,1]
FRONT = [0,0,0,5]
RED = [1,0,0,0.3]
BLUE = [0,0,1,0.3]



def show(render,f):

  sources = f.sources

  render.ctx.set_source_rgba(*FRONT)

  render.clear_canvas()

  for s in sources:

    render.circle(*s, r=ONE, fill=True)

  render.ctx.set_source_rgba(*BLUE)

  for cracks in f.cracks:

    s,_ = cracks[0]
    # render.ctx.move_to(*sources[s,:].flatten())
    render.circle(*sources[s,:].flatten(), r=3*ONE, fill=True)

    render.ctx.set_source_rgba(*RED)

    for c,_ in cracks[1:]:
      # render.ctx.line_to(*sources[c,:].flatten())
      render.circle(*sources[c,:].flatten(), r=3*ONE, fill=True)
      # print(*sources[s,:].flatten())

    render.ctx.stroke()


def step(f):

  f.fracture()

  return True


i = 0

def main():

  import gtk
  from render.render import Animate

  from modules.fracture import Fracture

  F = Fracture(INIT_NUM, INIT_RAD, INIT_DST)


  def wrap(render):

    global i

    res = step(F)
    show(render,F)

    i += 1

    return res

  render = Animate(SIZE, BACK, FRONT, wrap)
  render.set_line_width(LINEWIDTH)
  gtk.main()


if __name__ == '__main__':

  main()

