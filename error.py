class MatrixSizeError(ValueError):
  """
    Ошибка размера матрицы
  """
  pass

class MatrixDimensionMismatchError(ValueError):
  """
    Ошибка несоответствия измерений матрицы
  """
  pass

class MatrixMultiplicationError(ValueError):
  """
    Ошибка умножения матриц
  """
  pass

class MatrixShapeMismatchError(ValueError):
  """
    Ошибка несоответствия формы матрицы
  """
  pass

class MatrixNotSquareError(ValueError):
  """
    Ошибка неквадратной матрицы
  """
  pass

class SingularMatrixError(ValueError):
  """
    Ошибка сингулярной матрицы
  """
  pass

class MatrixReshapeError(ValueError):
  """
    Ошибка изменения формы матрицы
  """
  pass

class LUDCompositionError(ValueError):
  """
    Ошибка LU-разложения
  """
  pass

class LinearDependenceError(ValueError):
  """
   Ошибка линейной зависимости
  """
  pass