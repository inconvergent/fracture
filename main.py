#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import pi

TWOPI = pi*2.

NMAX = 10**6
SIZE = 800
ONE = 1./SIZE

INIT_NUM = 1000
INIT_RAD = 0.4
INIT_DST = 10*ONE

MID = 0.5

LINEWIDTH = 5.*ONE

BACK = [1,1,1,1]
FRONT = [0,0,0,5]
RED = [1,0,0,0.3]


i = 0

def show(render,f):

  sources = f.sources

  render.ctx.set_source_rgba(*RED)

  for s in sources:

    render.circle(*s, r=ONE, fill=True)

  return

def step(f):

  f.crack()

  return True


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

