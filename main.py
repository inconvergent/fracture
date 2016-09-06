#!/usr/bin/python3
# -*- coding: utf-8 -*-




BACK = [1,1,1,1]
FRONT = [0,0,0,0.8]
LIGHT = [0,0,0,0.2]
CYAN = [0,0.5,0.5,0.2]
BLUE = [0,0,1,0.3]


NMAX = 10**6
SIZE = 1200
ONE = 1./SIZE
LINEWIDTH = ONE*1.1

INIT_NUM = 20000
INIT_RAD = 0.45

SOURCE_DST = 2.0*ONE

FRAC_DOT = 0.85
FRAC_DST = 100.*ONE
FRAC_STP = ONE*2
FRAC_SPD = 1.0

FRAC_DIMINISH = 0.997
FRAC_SPAWN_DIMINISH = 0.9


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

  from iutils.render import Animate
  from modules.fracture import Fractures

  # from dddUtils.ioOBJ import export_2d as export
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

  # uniform square distribution
  from numpy.random import random
  for _ in range(5):
    F.blow(2, random(size=2))

  # uniform circular distribution
  # for _ in xrange(5):
    # F.blow(3, random_uniform_circle(INIT_RAD, num=1))

  def wrap(render):

    if not F.i % 20:
      show(render,F)
      # vertices, paths = F.get_vertices_and_paths()
      # export('fractures', fn.name(), vertices, lines=paths)
      render.write_to_png(fn.name()+'.png')

    F.print_stats()
    res = F.step(dbg=False)
    n = F.spawn_front(factor=SPAWN_FACTOR, angle=SPAWN_ANGLE)
    print('spawned: {:d}'.format(n))

    # fn = './asdf_{:04d}.png'.format(F.i)
    # render.write_to_png(fn)

    # if not res:
      # vertices, paths = F.get_vertices_and_paths()
      # export('fractures', fn.name(), vertices, lines=paths)

    return res

  render = Animate(SIZE, BACK, FRONT, wrap)
  render.start()


if __name__ == '__main__':

  main()

