import csv
import random
import numpy as np
from numpy import linalg as LA
from fractions import Fraction


def matrix_input(n_f: object, m_f: object) -> object:
    """"
    Функция, для ввода матрицы с клавиатуры

    :param n_f:  int -- количество строк в матрице
    :param m_f:  int -- количество столбцов в матрице

    :return np.array(matrix_f): Матрица

    """
    matrix_1 = []  # матрица, которую мы вводим в функции matrix_input
    for _ in range(n_f):
        line = []  # строка матрицы
        for index_f in range(m_f):
            number = str(input('Введите число '))  # введённое число в матрице
            if 'j' in number:
                plus = number.find('+')  # позиция + в комплексном числе
                minus = number.find('-')  # позиция – в комплексном числе
                ed = number.find('j')  # позиция j в комплексном числе
                if plus != -1:
                    if plus > ed:
                        sign = number[plus:plus + 1]
                        edin = number[:ed + 1]
                        number = number[plus + 1:]
                        number = number + sign + edin
                else:
                    if minus > ed:
                        sign = number[minus:minus + 1]
                        edin = number[:ed + 1]
                        number = number[minus + 1:]
                        number = number + sign + edin
            line.append(number)
        matrix_1.append(line)
    for _ in range(n_f):
        for index_f in range(m_f):
            matrix_1[_][index_f] = type_number(matrix_1[_][index_f])
    return np.array(matrix_1)


def type_number(number):
    """"
    Функция, которая преобразует число в нужный тип

    :param number:  str -- число, которое мы считываем

    :return number: число с нужным типом

    """
    if 'j' in number:
        number = complex(number)
    elif '.' in number:
        number = float(number)
    else:
        number = int(number)
    return number


def fractional_matrix(matrix_f: object) -> object:
    """"
    Функция, которая преобразует матрицу в тип fractions

    :param matrix_f: list -- матрица

    :return np.array(matrix_f): Матрица

    """
    for _ in range(len(matrix_f)):  # _ - индекс строки
        for index in range(len(matrix_f[_])):  # index - индекс столбца
            matrix_f[_][index] = Fraction(str(matrix_f[_][index]))
    return np.array(matrix_f)


def gauss(matrix_f):
    """
    Функция, которая решает СЛАУ методом Гаусса

    :param matrix_f: list -- матрица, соответствующая системе

    :return matrix_f:Решенная матрица

    """
    n_f = len(matrix_f)  # n_f - количество строк в матрице
    for k in range(n_f):  # k - количество итераций
        for _ in range(k, n_f):  # _ - индекс строки
            if _ == k:
                matrix_f[k][len(matrix_f[_]) - 1] = (
                        matrix_f[k][len(matrix_f[_]) - 1] / matrix_f[k][k])
            else:
                matrix_f[_][len(matrix_f[_]) - 1] = (
                    matrix_f[_][len(matrix_f[_]) - 1] -
                    ((matrix_f[k][len(matrix_f[_]) - 1]) * (matrix_f[_][k])))

            for index_f in range(len(matrix_f[_]) - 2, -1, -1):  # index_f -
                # индекс столбца
                if _ == k:
                    matrix_f[_][index_f] = (matrix_f[_][index_f] /
                                            matrix_f[_][_])
                else:
                    matrix_f[_][index_f] = (matrix_f[_][index_f] -
                    matrix_f[k][index_f] * matrix_f[_][k])
    return matrix_f


def jacobi(A, b, eps=1e-3):
    """"
    Функция, которая решает СЛАУ методом Якоби

    :param A: list-- матрица, соответствующая системе, ез вектора свободных членов
    :param b: list -- вектор свободных членов
    :param eps: float -- точность

    :return x1: вектор неизвестных

    """
    x1 = np.zeros(len(A[0]))  # вектор неизвестных
    x = np.ones(len(A[0]))  # вектор из единиц
    d = np.diag(A)  # диагональ матрицы
    u = A - np.diagflat(d)  # матрица без диагонали
    d1 = np.diag(A)  # диагональ матрицы
    u1 = A - np.diagflat(d1)  # матрица без диагонали
    print(f"Диалогальное преобладание A1: {np.sum(abs(d1)) > np.sum(abs(u1))}")

    while max(abs(x1 - x)) > eps:
        x = x1
        x1 = (b - np.dot(u, x)) / d
    return x1


n = int(input('Введите количетво строк матрицы '))
m = n + 1
print(
    "Как будете вводить числа: 1) с клавиатуры,2) заполнение случайными "
    "числами (рандомная генерация),3) считывание "
    "из csv файла")
input_method = int(input())
if input_method == 3:
    print('Введите диапазон, из которого будут генерироваться числа')
    a = int(input('a = '))
    b = int(input('b = '))
    print("Какие числа хотите генерировать?")
    print("1-Целые")
    print("2-Дробные")
    print("3-Комплексные")
    type_of_numbers = int(input("Введите число "))  #
    if type_of_numbers == 1:
        with open('mrx.csv', mode="w", encoding='utf-8') as file:
            file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
            for _ in range(n):
                input_numbers = []
                for index in range(m):
                    input_numbers.append(int(random.uniform(a, b)))
                file_writer.writerow(input_numbers)
    elif type_of_numbers == 2:
        with open("mrx.csv", mode="w", encoding='utf-8') as file:
            file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
            for _ in range(n):
                input_numbers = []
                for index in range(m):
                    input_numbers.append(float(random.uniform(a, b)))
                file_writer.writerow(input_numbers)
    else:
        with open("mrx.csv", mode="w", encoding='utf-8') as file:
            file_writer = csv.writer(file, delimiter=",", lineterminator="\r")
            for _ in range(n):
                input_numbers = []
                for index in range(m):
                    input_numbers.append(int(random.uniform(a, b)))
                file_writer.writerow(input_numbers)

    matrix_string = []
    with open("mrx.csv") as file:
        reader = csv.reader(file, delimiter=',')
        for row in reader:
            matrix_string.append(row)
    if type_of_numbers == 1:
        matrix = np.array(matrix_string, dtype=int)
        print(matrix)
    elif type_of_numbers == 2:
        matrix = np.array(matrix_string, dtype=float)
        print(matrix)
    else:
        matrix_string = np.array(matrix_string, dtype=int)
        matrix = np.array(matrix_string, dtype=complex)
elif input_method == 2:
    matrix = []
    print('Введите диапазон, из которого будут генерироваться числа')
    a = int(input('a = '))
    b = int(input('b = '))
    print("Какие числа хотите генерировать?")
    print("1-Целые")
    print("2-Дробные")
    print("3-Комплексные")
    type_of_numbers = int(input("Введите число"))
    if type_of_numbers == 1:
        matrix = []
        for _ in range(n):
            input_numbers = []
            for ndex in range(m):
                input_numbers.append(int(random.uniform(0, 10)))
            matrix.append(input_numbers)
        matrix = np.array(matrix)
    elif type_of_numbers == 2:
        matrix = []
        for _ in range(n):
            input_numbers = []
            for index in range(m):
                input_numbers.append(float(random.uniform(0, 10)))
            matrix.append(input_numbers)
        matrix = np.array(matrix)
    else:
        matrix = []
        for _ in range(n):
            a1 = []
            for index in range(m):
                a1.append(complex(random.uniform(0, 10)))
            matrix.append(a1)
        matrix = np.array(matrix)

else:
    matrix = []
    print('Введите матрицу системы')
    matrix = matrix_input(n, m)
print("Исходная матрица")
print(matrix)

matrix21 = jacobi(matrix[:, :-1], matrix[:, -1], eps=1e-3)

inverse_matrix = np.linalg.inv(matrix[:, :-1])
matrix_of_coefficients = np.linalg.inv(inverse_matrix)
# cord
a = [0 for i in range(n)]
a_obr = [0 for i in range(n)]
for _ in range(n):
    for j in range(n):
        a_obr[_] += abs(inverse_matrix[_][j])
        a[_] += abs(matrix_of_coefficients[_][j])

norma_A = max(a)
norma_A_obr = max(a_obr)
cord_A = norma_A * norma_A_obr

diagonal = np.diag(matrix[:, :-1])
no_diagonal = matrix[:, :-1] - np.diagflat(diagonal)

if (cord_A <= 100) & (np.sum(abs(diagonal)) > np.sum(abs(no_diagonal))):
    print("СЛАУ")
    print(matrix21)
    print("Обратная матрица")
    print(inverse_matrix)
    print("Прямая матрица коэффициентов")
    print(matrix_of_coefficients)
    print("cord")
    print(cord_A)
    print("Методом простых итераций Якоби")
else:
    matrix2 = gauss(matrix)
    matrix_A = []
    matrix_B = []
    for _ in range(n):
        c = []
        b = []
        for j in range(m - 1):
            c.append(matrix2[_][j])
            b.append(matrix[_][j])
        matrix_A.append(c)
        matrix_B.append(matrix2[_][j + 1])

    matrix_A_obr = LA.inv(matrix_A)
    X = np.dot(matrix_A_obr, matrix_B)

    a = [0 for i in range(n)]
    a_obr = [0 for i in range(n)]
    for _ in range(n):
        for j in range(n):
            a[_] += abs(matrix_A[_][j])
            a_obr[_] += abs(matrix_A_obr[_][j])

    norma_AB = max(a)
    norma_AB_obr = max(a_obr)
    cord_B = norma_AB * norma_AB_obr
    if cord_B <= 100:
        print('Прямая матрица коэффициентов (A)')
        print(matrix_A)

        print('Обратная матрица для матрицы коэффициентов (A-1)')
        print(matrix_A_obr)

        print('Столбец неизвестных – решение СЛАУ (X)')
        print(X)

        print('Cord(A) = ', cord_B)
        print("Прямой алгоритм Гаусса-Жордана")

    else:
        matrix3 = gauss(fractional_matrix(matrix))
        matrix_AA = []
        matrix_BB = []
        for _ in range(n):
            c = []
            b = []
            for j in range(m - 1):
                c.append(matrix3[_][j])
                b.append(matrix[_][j])
            matrix_AA.append(c)
            matrix_BB.append(matrix3[_][j + 1])

        matrix_A_obr1 = LA.inv(matrix_AA)
        Xx = np.dot(matrix_A_obr1, matrix_BB)

        a = [0 for i in range(n)]
        a_obr = [0 for i in range(n)]
        for _ in range(n):
            for j in range(n):
                a[_] += abs(matrix_AA[_][j])
                a_obr[_] += abs(matrix_A_obr1[_][j])

        norma_AC = max(a)
        norma_AC_obr = max(a_obr)
        cord_AC = norma_AC * norma_AC_obr

        print("Прямая матрица коэффициентов (A)")
        print(matrix_AA)

        print("Обратная матрица для матрицы коэффициентов (A-1)")
        print(matrix_A_obr1)

        print("Столбец неизвестных – решение СЛАУ (X)")
        print(Xx)

        print('Cord(A) = ', cord_AC)
        print("Прямой алгоритм Гаусса Жордана, с вычислениями правильных "
              "дробей")
