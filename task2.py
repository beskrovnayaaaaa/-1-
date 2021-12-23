import random
import math
import cmath
import time
import numpy as np
import pandas as pd


def beautiful(matrix):
    """
    Приводит список списков к виду матрицы.
    Например, [[1, 2, 3], [1, 2, 3], [1, 2, 3]] к
    [1, 2, 3]
    [1, 2, 3]
    [1, 2, 3]
    :param matrix: list - введенная матрица
    :return: списки по строкам, матрица
    """
    for i in range(0, len(matrix)):
        print(matrix[i])
    return ''


def el_int(el):
    """
    Проверка на принадлежность к типу int
    :param el: вводится элемент любого типа
    :return: True, если элемент типа int, иначе False
    """
    try:
        int(el)
        return True
    except ValueError:
        return False


def el_float(el):
    """
    Проверка на принадлежность к типу float
    :param el: вводится элемент любого типа
    :return: True, если элемент типа float, иначе False
    """
    try:
        float(el)
        return True
    except ValueError:
        return False


def el_complex(el):
    """
    Проверка на принадлежность к типу complex
    :param el: вводится элемент любого типа
    :return: True, если элемент типа complex, иначе False
    """
    try:
        complex(el)
        return True
    except ValueError:
        return False


def rev_complex(h):
    """
    Функция преобразует комплексное число в нормальный вид, т. е. в вид a + i*b
    Пример: если вы ввели -j + 1, функция преобразует это в 1 - j
    :param h: str — элемент матрицы
    :return: str — преобразованный элемент
    """
    h_rev = ''
    sep = 0
    if h[0] == '+' or h[0] == '-':
        for i in range(1, len(h)):
            if h[i] == '+' or h[i] == '-':
                sep = i
                break
        h_rev = h[sep:len(h)] + h[0:sep]
    else:
        for i in range(0, len(h)):
            if h[i] == '+' or h[i] == '-':
                sep = i
                break
        h_rev = h[sep:len(h)] + '+' + h[0:sep]
    return (h_rev)


def matr(m, n):
    """
    Делает матрицу из введенных пользователем строк матрицы
    :param m: int - кол-во строк матрицы
    :param n: int - кол-во столбцов матрицы
    :return: выводит список списков, где каждый список это строка матрицы
    """
    matrix = []
    print('Введите элементы строки матрицы через пробел:')
    for i in range(0, m):
        a = []
        row = input()
        row = row.split(' ')
        matrix.append(row)
        if len(row) != n:
            print('Некорректное количество элементов в строке матрицы.')
            exit()
        for j in range(0, n):
            el = matrix[i][j]
            k = 0
            while (k == 0):
                if el_int(el):
                    matrix[i][j] = int(el)
                    k = 1
                else:
                    if el_float(el):
                        matrix[i][j] = float(el)
                        k = 1
                    else:
                        if el_complex(el):
                            matrix[i][j] = complex(el)
                            k = 1
                        else:
                            if el_complex(rev_complex(el)):
                                matrix[i][j] = complex(rev_complex(el))
                                k = 1
                            else:
                                el = input('Неверный формат ввода.\
                                            Повторите ввод элемента [{}, {}]:\
                                            '.format(i, j))
    return (matrix)


def transp(matrix):
    """
    Транспонирование матрицы
    :param matrix: list - введенная матрица
    :return: list - транспонированную матрицу
    """
    tr = [[0 for j in range(0, len(matrix))] for i in range(0, len(matrix[0]))]
    for i in range(0, len(matrix[0])):
        for j in range(0, len(matrix)):
            tr[i][j] = matrix[j][i]
    return beautiful(tr)


def mult(a, b):
    """
    Произведение матрицы на матрицу.
    :param a: list - первая матрица
    :param b: list - вторая матрица
    :return: list - итоговая матрица
    """
    s = 0  # сумма
    h = []
    c = []
    for i in range(0, len(a)):
        for j in range(0, len(b[0])):
            for k in range(0, len(a[0])):
                s = s + a[i][k] * b[k][j]
            h.append(s)
            s = 0
        c.append(h)
        h = []
    return c


def num_mult(matrix, num):
    """
    Произведение матриы на число.
    :param matrix: list -  матрица
    :param num: число
    :return: list - матрица
    """
    if el_int(num):
        num = int(num)
    else:
        if el_float(num):
            num = float(num)
        else:
            if el_complex(num):
                num = complex(num)
            else:
                if el_complex(rev_complex(num)):
                    num = complex(rev_complex(num))
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            matrix[i][j] = matrix[i][j] * num
    return matrix


def add(matrix, new_matrix):
    """
    Сложение двух матриц
    :param matrix: list - первая матрица
    :param new_matrix: list - вторая матрица
    :return: list - итоговая матрица
    """
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            matrix[i][j] = matrix[i][j] + new_matrix[i][j]
    return matrix


def sub(matrix, new_matrix):
    """
    Разность матриц
    :param matrix: list - первая матрица
    :param new_matrix: list - матрица, которую вычитают
    :return: итоговая матрица
    """
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            matrix[i][j] = matrix[i][j] - new_matrix[i][j]
    return matrix


def rec_matrix(matrix, n):
    """
    Предоставляет выбор действий с изначальной матрицей
    и выполняет действия.
    :param matrix: list - изначчальная матрица
    :param n: int - кол-во строк
    :return: list - результат действия вашего выбора
    (получившуюся матрицу) и предлагает заново выбор действий
    """
    print('выберите операцию:')
    print('1)Умножение матрицы на число')
    print('2)Умножение матрицы на матрицу')
    print('3)Прибавить матрицу')
    print('4)Вычесть матрицу')
    print('5)Конец ввода')
    choice = int(input())
    if choice == 5:
        return matrix
    else:
        if choice == 2:
            new_matrix = matr(n, n)
            matrix = mult(matrix, new_matrix)
        if choice == 1:
            num = input('Введите число: ')
            matrix = num_mult(matrix, num)
        if choice == 3:
            new_matrix = matr(n, n)
            matrix = add(matrix, new_matrix)
        if choice == 4:
            new_matrix = matr(n, n)
            matrix = sub(matrix, new_matrix)
        print('Результат ')
        print(beautiful(matrix))
        return (rec_matrix(matrix, n))


def det_matrix(matrix):
    """
    Нахождение определителя
    :param matrix: list - матрица
    :return:
    """
    d = 0
    if len(matrix) > 1:
        for i in range(0, len(matrix)):
            new_matrix = [[matrix[g][j] for j in range(0, len(matrix))
                           if j != i] for g in range(1, len(matrix[0]))]
            d = d + (-1) ** (i + 1 + 1) * matrix[0][i] * det_matrix(new_matrix)
        return d
    else:
        return matrix[0][0]


A1 = [[5, 8, 9], [13, 21, 3], [11, 22, 33]]
A2 = [[26, 35, 1 + 1j, 86], [45, 33, 21, 55]]
A3 = [[3, 8, 10], [2, 0, 4], [5, 9, 11]]
A4 = [[1 - 1j, 2], [4, 6 + 3j]]


random_matrix = pd.DataFrame([[int(random.random() * 100)
                               for _ in range(10000)] for _ in range(10000)])

random_matrix.to_csv('random_matrix.csv', header=True, index=False)

random_matrix = pd.read_csv('random_matrix.csv')
spisok = random_matrix.values.tolist()

random_matrix1 = pd.DataFrame([[int(random.random() * 100)
                               for _ in range(10)] for _ in range(10)])

random_matrix1.to_csv('random_matrix1.csv', header=True, index=False)

random_matrix1 = pd.read_csv('random_matrix1.csv')
spisok1 = random_matrix1.values.tolist()


print('Ввод матрицы.')
m_ch = False
while(not m_ch):
    m = input('Введите число строк: ')
    try:
        m = int(m)
        if(m <= 0):
            print("Ошибка. Число строк - \
                положительная величина. Повторите ввод.")
            m_ch = False
        else:
            m_ch = True
    except:
        print("Ошибка. Неправильный формат ввода. Повторите ввод.")
        m_ch = False


n_ch = False
while(not n_ch):
    n = input('Введите число столбцов: ')
    try:
        n = int(n)
        if(n <= 0):
            print("Ошибка. Число столбцов - \
                 положительная величина. Повторите ввод.")
            n_ch = False
        else:
            n_ch = True
    except:
        print("Ошибка. Неправильный формат ввода. Повторите ввод.")
        n_ch = False

matrix = matr(m, n)
print('Введенная матрица: ')
print(beautiful(matrix))
print(beautiful(rec_matrix(matrix, n)))
print('Транспонированная матрица: ')
print(transp(matrix))
print("Определитель матрицы")
det_matrix(A3)
print('Транспонирование')
start = time.perf_counter_ns()
print(np.transpose(A2))
end = time.perf_counter_ns()
print(f"Время выполнения встроенной функции = {end - start}")
start = time.perf_counter_ns()
print(transp(A2))
end = time.perf_counter_ns()
print(f"Время выполнения созданной функции = {end - start}")
print('Транспонирование матрицы 10000*10000')
start = time.perf_counter_ns()
np.transpose(spisok1)
end = time.perf_counter_ns()
print(f"Время выполнения встроенной функции = {end - start}")
start = time.perf_counter_ns()
transp(spisok1)
end = time.perf_counter_ns()
print(f"Время выполнения созданной функции = {end - start}")
print('Умножение матриц')
start = time.perf_counter_ns()
print(np.array(A1).dot(np.array(A3)))
end = time.perf_counter_ns()
print(f"Время выполнения встроенной функции = {end - start}")
start = time.perf_counter_ns()
print(mult(A1, A3))
end = time.perf_counter_ns()
print(f"Время выполнения созданной функции = {end - start}")
print('Вычисление определителя')
start = time.perf_counter_ns()
print(np.linalg.det(np.array(A3)))
end = time.perf_counter_ns()
print(f"Время выполнения встроенной функции = {end - start}")
start = time.perf_counter_ns()
print(det_matrix(A3))
end = time.perf_counter_ns()
print(f"Время выполнения созданной функции = {end - start}")
print('Вычисление определителя')
start = time.perf_counter_ns()
print(np.linalg.det(np.array(spisok1)))
end = time.perf_counter_ns()
print(f"Время выполнения встроенной функции = {end - start}")
start = time.perf_counter_ns()
print(det_matrix(spisok1))
end = time.perf_counter_ns()
print(f"Время выполнения созданной функции = {end - start}")
