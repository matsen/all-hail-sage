import csv
import os
from IPython.display import display, Image
from sage.matrix.constructor import Matrix
from sage.misc.temporary_file import tmp_dir, graphics_filename


dirname = tmp_dir()


def nb_show(graphics_object, ext='.png'):
    name = os.path.join(dirname, graphics_filename(ext))
    graphics_object.save(name)
    display(Image(name))


def matrix_of_csv(fname):
    with open(fname, 'rb') as csvfile:
        mreader = csv.reader(csvfile, delimiter=',', quotechar="'")
        return Matrix([[int(x) for x in row] for row in mreader])
