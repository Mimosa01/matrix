from enum import Enum

class Direction(Enum):
  COL = 'col'
  ROW = 'row'


class ReturnedTypeMatrix(Enum):
  LIST = 'list'
  MATRIX = 'vector'