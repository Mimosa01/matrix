# Реализация класса матриц

# Дополнительно:
# Реализовать методы разложения:
#  сингулярное, гауссово

from typing import List, Union, Tuple
from error import *
from enums import *
from math import sqrt, fabs


class Matrix:
  def __init__(self, *, rows: int = None, columns: int = None, data: List[List[Union[int, float]]] = None) -> None:
    if data is not None:
      self.matrix = data
      self.rows = len(data)
      self.columns = len(data[0]) if data else 0
    elif rows is not None and columns is not None and data is None:
      if rows <= 0 and columns <= 0:
        raise MatrixSizeError("Количество строк и столбцов должно быть положительным числом")
      self.rows = rows
      self.columns = columns
      self.matrix = [[0] * columns for _ in range(rows)]
    else:
      self.matrix = []
      self.rows = 0
      self.columns = 0

  def __add__(self, other: 'Matrix') -> 'Matrix':
    if self.rows != other.rows or self.columns != other.columns:
      raise MatrixDimensionMismatchError("Сложение матриц работает только с одинаковыми по размеру матрицами")
    
    matrix = []

    for i_row in range(self.rows):
      row = []
      for i_column in range(self.columns):
        row.append(self.matrix[i_row][i_column] + other.matrix[i_row][i_column])
      matrix.append(row)

    return Matrix(data=matrix)

  def __mul__(self, other: Union['Matrix', float, int]) -> 'Matrix':
    matrix = []
    if isinstance(other, (float, int)):
      for i_row in range(self.rows):
        row = []
        for i_column in range(self.columns):
          row.append(self.matrix[i_row][i_column] * other)
        matrix.append(row)
    elif isinstance(other, Matrix):
      if self.columns == other.rows:
        for i in range(self.rows):
          row = [sum(self.matrix[i][k] * other.matrix[k][j] for k in range(self.columns)) for j in range(other.columns)]
          matrix.append(row)
        return Matrix(data=matrix)
      else:
        raise MatrixShapeMismatchError("Умножение матриц возможно, только если количество столбцов первой равно количеству строк второй")
    else:
      raise MatrixMultiplicationError("Умножение матрицы возможно только на матрицу или число")
    return Matrix(data=matrix)

  def __sub__(self, other: 'Matrix') -> 'Matrix':
    if self.rows != other.rows or self.columns != other.columns:
      raise MatrixDimensionMismatchError("Вычитание матриц работает только с одинаковыми по размеру матрицами")
    
    matrix = []

    for i_row in range(self.rows):
      row = []
      for i_column in range(self.columns):
        row.append(self.matrix[i_row][i_column] - other.matrix[i_row][i_column])
      matrix.append(row)

    return Matrix(data=matrix)

  def __str__(self):
    return '\n'.join([' '.join(map(str, row)) for row in self.matrix])
  
  @classmethod
  def identity(cls, rows: int, cols: int) -> 'Matrix':
    if rows != cols:
      raise MatrixNotSquareError("Только для квадратной матрицы строиться единичная")
    I = [[0] * rows for _ in range(rows)]
    for i in range(rows):
      I[i][i] = 1
    return Matrix(data=I)
  
  @classmethod
  def _round_matrix(self, *, matrix: 'Matrix', precision: int = 2) -> 'Matrix':
    return Matrix(data=[[round(cell, precision) for cell in row] for row in matrix.matrix])

  def l2_norm(self, *, index: int, vector_type: Direction) -> float:
    if index is None:
      raise IndexError("Не передан индекс строки или колонки")
    if (vector_type == Direction.ROW and index >= self.rows) or \
        (vector_type == Direction.COL and index >= self.columns) or \
        index < 0:
      raise MatrixSizeError("Номер столбца или колонки не должен превышать размер матрицы или быть отрицательным")
    
    vector = self.get_row(row_index=index, returned_type=ReturnedTypeMatrix.LIST) \
      if vector_type == Direction.ROW else \
      self.get_column(col_index=index, returned_type=ReturnedTypeMatrix.LIST)
    
    return sqrt(sum([i ** 2 for i in vector]))
  
  @classmethod
  def dot_product(cls, v1: Union['Matrix', List[Union[float, int]]], v2: Union['Matrix', List[Union[float, int]]]) -> float:
    def extract_vector(v) -> List[Union[float, int]]:
      if isinstance(v, Matrix):
        if v.rows > 1 and v.columns == 1:
          return v.transpose().to_list()
        elif v.columns > 1 and v.rows == 1:
          return v.to_list()
        else:
          raise MatrixShapeMismatchError("Матрицу невозможно преобразовать в вектор")
      elif isinstance(v, List):
        return v
      else:
        raise TypeError("Аргумент должен быть либо матрицей, либо списком")
    
    vector_1 = extract_vector(v1)
    vector_2 = extract_vector(v2)
      
    if len(vector_1) != len(vector_2):
      raise MatrixDimensionMismatchError("Векторы должны быть одинаковой длины")
    if len(vector_1) == 0:
      return 0
    
    return sum(a * b for a, b in zip(vector_1, vector_2))
  
  def insert_row_column(self, *, index: int, vector: List[Union[float, int]], direction: Direction) -> 'Matrix':
    if index is None:
      raise IndexError("Не передан индекс строки или колонки")
    if (direction == Direction.ROW and index >= self.rows) or \
        (direction == Direction.COL and index >= self.columns) or \
        index < 0:
      raise MatrixSizeError("Номер столбца или колонки не должен превышать размер матрицы или быть отрицательным")
    if (direction == Direction.ROW and len(vector) != self.rows) or \
        (direction == Direction.COL and len(vector) != self.columns):
      raise MatrixDimensionMismatchError('Размер вектора не равен размеру матрицы')
    if direction not in (Direction.ROW, Direction.COL):
      raise ValueError(f"Недопустимое направление: {direction}")
    
    matrix = []

    if direction == Direction.ROW:
      matrix = [(row[:] if i_row != index else vector) for i_row, row in enumerate(self.matrix)]
      return Matrix(data=matrix)
    
    matrix = [
        row[:index] + [vector[i]] + row[index+1:]
        for i, row in enumerate(self.matrix)
      ] 
    return Matrix(data=matrix)

  def lu_decomposition(self) -> Tuple['Matrix', 'Matrix', 'Matrix']:
    if self.rows != self.columns:
      raise MatrixNotSquareError("LU-разложение возможно только для квадратных матриц")

    L = Matrix.identity(rows=self.rows, cols=self.columns)
    U = Matrix(rows=self.rows, columns=self.columns)
    P = Matrix.identity(rows=self.rows, cols=self.columns)
    A_copy = Matrix(data=[row[:] for row in self.matrix])

    for i in range(self.rows):
      max_row = i
      for k in range(i + 1, self.rows):
        if abs(A_copy.matrix[k][i]) > abs(A_copy.matrix[max_row][i]):
          max_row = k

      if max_row != i:
        A_copy.matrix[i], A_copy.matrix[max_row] = A_copy.matrix[max_row], A_copy.matrix[i]
        P.matrix[i], P.matrix[max_row] = P.matrix[max_row], P.matrix[i]

      for j in range(i, self.columns):
        U.matrix[i][j] = A_copy.matrix[i][j] - sum(L.matrix[i][k] * U.matrix[k][j] for k in range(i))

      for j in range(i + 1, self.rows):
        if abs(U.matrix[i][i]) < 1e-9:
          raise LUDCompositionError("LU-разложение невозможно, элемент на диагонали равен 0")
        L.matrix[j][i] = (A_copy.matrix[j][i] - sum(L.matrix[j][k] * U.matrix[k][i] for k in range(i))) / U.matrix[i][i]

      for i in range(1, self.rows):
        if P.matrix[i][i] == 0:
          U.matrix[i][i] = -U.matrix[i][i]
    return L, U, P

  def qr_decomposition(self) -> Tuple['Matrix', 'Matrix']:
    if self.rows < self.columns:
      raise MatrixDimensionMismatchError("QR-разложение требует, чтобы число строк было не меньше числа столбцов.")
    
    Q = Matrix(rows=self.rows, columns=self.columns)
    R = Matrix(rows=self.columns, columns=self.columns)
    
    for col in range(self.columns):
      v_col = self.get_column(col_index=col, returned_type=ReturnedTypeMatrix.LIST)
      
      for i in range(col):
        r_ij = sum(Q.matrix[k][i] * v_col[k] for k in range(self.rows))
        R.matrix[i][col] = r_ij
        v_col = [v_col[k] - r_ij * Q.matrix[k][i] for k in range(self.rows)]
      
      norm_v_col = sqrt(sum(x**2 for x in v_col))
      if norm_v_col < 1e-10:
        raise LinearDependenceError("Столбцы матрицы линейно зависимы")

      q_col = Matrix.from_list([[x / norm_v_col for x in v_col]])

      Q = Q.insert_row_column(index=col, vector=q_col.to_list(), direction=Direction.COL)
      R.matrix[col][col] = norm_v_col

    return Q, R

  def normalize_and_sign_correction(self) -> 'Matrix':
    V = Matrix(data=[row[:] for row in self.matrix])
    for i in range(V.rows):
      norm = sqrt(sum([V.matrix[i][j]**2 for j in range(V.columns)]))
      for j in range(V.columns):
        V.matrix[i][j] /= norm
      
      if V.matrix[i][0] < 0:
        for j in range(V.columns):
          V.matrix[i][j] *= -1
    return V

  def spectral_decomposition(self, max_iter: int = 1000, tol: float = 1e-16):
    if self.rows != self.columns:
      raise MatrixNotSquareError("Матрица должна быть квадратной")

    A_copy = Matrix(data=[row[:] for row in self.matrix])
    V = Matrix.identity(rows=self.rows, cols=self.columns)
    
    for _ in range(max_iter):
      Q, R = A_copy.qr_decomposition()
      A_copy = R * Q 
      V = V * Q

      off_diagonal_sum = sum([A_copy.matrix[i][j]**2 for i in range(A_copy.rows) for j in range(A_copy.columns) if i != j])
      
      if off_diagonal_sum < tol:
        break
    # V = V.normalize_and_sign_correction()
    return V, A_copy, V.inverse()

  def svd_decomposition(self):
    if self.rows != self.columns:
      raise MatrixNotSquareError("Матрица должна быть квадратной")

    E = Matrix(rows=self.rows, columns=self.columns)
    A_1 = self * self.transpose()
    A_2 = self.transpose() * self
    V, E_1, _ = A_2.spectral_decomposition()
    U, _, _ = A_1.spectral_decomposition()

    for i in range(self.columns):
      E.matrix[i][i] = sqrt(E_1.matrix[i][i])
    
    return U, E, V.transpose()

  def transpose(self) -> 'Matrix':
    return Matrix(data=list(map(list, zip(*self.matrix))))

  def determinant(self) -> float:
    det = 1
    _, U, P = self.lu_decomposition()
    
    permutation_count = sum(1 for i in range(self.rows) if P.matrix[i][i] != 1)
    det_P = (-1) ** permutation_count
    for i in range(U.rows):
      det *= U.matrix[i][i]
    
    return round(det * det_P, 2) if det * det_P != 0 else 0 

  def inverse(self) -> 'Matrix':
    if self.rows != self.columns:
      raise MatrixNotSquareError("Получение обратной матрицы возможно только из квадратной матрицы")
    if self.determinant() == 0:
      raise SingularMatrixError("Матрица вырождена")

    L, U, P = self.lu_decomposition()
    X, Y = [[0] * self.rows for _ in range(self.rows)], [[0] * self.rows for _ in range(self.rows)]
    I = Matrix.identity(rows=self.rows, cols=self.columns)

    for col in range(self.columns):
      Y[col][0] = I.matrix[col][0]
      for row in range(1, self.rows):
        Y[col][row] = I.matrix[col][row] - sum(L.matrix[row][k] * Y[col][k] for k in range(row))

    for col in range(self.columns-1, -1, -1):
      X[col][self.rows-1] = Y[col][self.rows-1] / U.matrix[self.rows-1][self.rows-1]
      for row in range(self.rows-2, -1, -1):
        X[col][row] = (Y[col][row] - sum(U.matrix[row][k] * X[col][k] for k in range(row + 1, self.rows))) / U.matrix[row][row]

    return Matrix(data=X).transpose() * P


  def get_row(self, *, row_index: int, returned_type: ReturnedTypeMatrix) -> Union['Matrix', List[Union[float, int]]]:
    if not (0 <= row_index < self.columns):
      raise IndexError(f"Индекс столбца {row_index} вне допустимого диапазона [0, {self.columns - 1}]")
    if returned_type == ReturnedTypeMatrix.MATRIX:
      return Matrix(data=self.matrix[row_index])
    elif returned_type == ReturnedTypeMatrix.LIST:
      return self.matrix[row_index]
    else:
      raise ValueError(f"Недопустимый тип возвращаемого значения: {returned_type}")

  def get_column(self, *, col_index: int, returned_type: ReturnedTypeMatrix) -> Union['Matrix', List[Union[float, int]]]:
    if not (0 <= col_index < self.columns):
      raise IndexError(f"Индекс столбца {col_index} вне допустимого диапазона [0, {self.columns - 1}]")
    if returned_type == ReturnedTypeMatrix.MATRIX:
      return Matrix(data=[[i[col_index]] for i in self.matrix])
    elif returned_type == ReturnedTypeMatrix.LIST:
      return [i[col_index] for i in self.matrix]
    else:
      raise ValueError(f"Недопустимый тип возвращаемого значения: {returned_type}")

  def minor(self, row_index: int, col_index: int) -> float:
    matrix = []

    for i_row in range(self.rows):
      if i_row != row_index:
        row = []
        for i_column in range(self.columns):
          if i_column != col_index:
            row.append(self.matrix[i_row][i_column])
        matrix.append(row)
    
    return Matrix(data=matrix).determinant()

  def cofactor(self, row_index: int, col_index: int, precision: int = 2) -> float:
    det = self.minor(row_index, col_index)
    degree = row_index + col_index
    return round((-1)**degree * det, precision)

  def adjugate(self) -> 'Matrix':
    adj = Matrix(rows=self.rows, columns=self.columns)

    for row in range(self.rows):
      for col in range(self.columns):
        adj.matrix[row][col] = self.cofactor(row, col)
    
    return adj.transpose()

  def rank(self) -> int:
    _, U = self.lu_decomposition()
    rank = self.rows
    for row in U.matrix:
      if all(list(map(lambda x: x == 0, row))):
        rank -= 1
    return rank

  def is_square(self) -> bool:
    return self.rows == self.columns

  def is_identity(self) -> bool:
    if self.rows != self.columns:
      return False
    
    for row in range(self.rows):
      for col in range(self.columns):
        if row == col and self.matrix[row][col] != 1:
          return False
        elif row != col and self.matrix[row][col] != 0:
          return False
    return True

  def is_singular(self) -> bool:
    return self.determinant() == 0

  def reshape(self, rows: int, cols: int) -> 'Matrix':
    if (self.rows * self.columns) != (rows * cols):
      raise MatrixReshapeError(f"Невозможно преобразовать в размерность {rows} x {cols}")
    
    lst = [self.matrix[row][col] for row in range(self.rows) for col in range(self.columns)]
    matrix = []
    
    for i in range(rows):
      matrix.append(lst[i * cols:(i + 1) * cols])

    return Matrix(data=matrix)

  def to_list(self) -> Union[List[Union[float, int]], List[List[Union[float, int]]]]:
    if self.rows == 1:
      return self.get_row(row_index=0, returned_type=ReturnedTypeMatrix.LIST)
    return [row[:] for row in self.matrix]

  @classmethod
  def from_list(cls, lst: List[List[int]]) -> 'Matrix':
    return cls(data=lst)


A = Matrix(data=[
  [4,0],
  [3,5]
])
